# ## Reusable neighbor geometry
#
# `neighborgeometry(grid, rings, needs)` builds only what the request implies.
# Distances and bearings are still exposed at the coarsest granularity at which
# they are invariant — one table for a planar rectilinear grid, one table per
# latitude row for a spherical one — but the tables are not built at all unless
# `Distance` or `Bearing` was asked for.
#
# Bearings are measured clockwise from north in the local tangent plane. On a
# sphere, this is the initial great-circle azimuth.

struct NeighborGeometry{G<:AbstractGridSpec,F<:Tuple,P}
    grid::G
    rings::NeighborRings
    fields::F   # The canonical, filtered request
    payload::P
end

# Rectilinear geometry payloads. `NoGeometry{L}` carries the neighbor count so
# that records can be zipped against a same-length tuple of `nothing`s at no
# cost, keeping one code path for "geometry" and "no geometry".
struct NoGeometry{L} end
struct UniformGeometry{T}
    table::T
end
struct RowGeometry{T}
    rows::T
    rowaxis::Int
end

@inline geometryat(g::NoGeometry{L}, I) where {L} = ntuple(Returns(nothing), Val(L))
@inline geometryat(g::UniformGeometry, I) = g.table
@inline geometryat(g::RowGeometry, I) = @inbounds g.rows[I[g.rowaxis]]

geometryat(geom::NeighborGeometry{<:RectilinearGrid}, I) =
    geometryat(geom.payload.geometry, I)

function neighborgeometry(grid::RectilinearGrid, rings::NeighborRings=NeighborRings(),
        needs=DEFAULT_NEEDS)
    request = requestfields(needs)
    logical = logicaloffsets(rings)
    offsets = map(o -> CartesianIndex(storageoffset(grid, o)), logical)
    geometry = _rectgeometry(wantsgeometry(request), grid, logical)
    return NeighborGeometry(grid, rings, request, (; offsets, geometry))
end

_rectgeometry(::Absent, grid::RectilinearGrid, logical) = NoGeometry{length(logical)}()
_rectgeometry(::Requested, grid::RectilinearGrid{<:Planar}, logical) =
    UniformGeometry(geometrytable(grid, logical))
_rectgeometry(::Requested, grid::RectilinearGrid{<:Spherical}, logical) =
    RowGeometry(geometryrows(grid, logical), yaxisnum(grid))

function geometrytable(grid::RectilinearGrid{<:Planar}, logical)
    sx, sy = abs.(grid.spacing)
    return map(logical) do (dx, dy)
        (distance=hypot(dx * sx, dy * sy), bearing=mod(atand(dx * sx, dy * sy), 360.0))
    end
end

function geometryrows(grid::RectilinearGrid{<:Spherical}, logical)
    R = grid.manifold.radius
    sx, sy = abs.(grid.spacing)
    return [
        map(logical) do (dx, dy)
            east = dx * deg2rad(sx) * cosd(lat) * R
            north = dy * deg2rad(sy) * R
            (distance=hypot(east, north), bearing=mod(atand(east, north), 360.0))
        end for lat in grid.lookups[2]
    ]
end

# Cell-grid geometry payloads. `NoCentroids` does no centroid work whatsoever;
# `OnDemandCentroids` computes the from-cell centroid once per visited cell and
# each neighbor centroid per edge, which is O(1) in grid size and therefore safe
# on grids that do not fit in memory. `StoredCentroids` is v2's table, reachable
# only through the explicit `precompute` opt-in below.
struct NoCentroids end
struct OnDemandCentroids{G}
    levelgrid::G
    radius::Float64
end
struct StoredCentroids{G,V}
    levelgrid::G
    radius::Float64
    centroids::V
end

_cellgeometry(::Absent, grid::CellGrid) = NoCentroids()
_cellgeometry(::Requested, grid::CellGrid) =
    OnDemandCentroids(grid.levelgrid, manifold(grid).radius)

function neighborgeometry(grid::CellGrid, rings::NeighborRings=NeighborRings(),
        needs=DEFAULT_NEEDS)
    request = requestfields(needs)
    return NeighborGeometry(grid, rings, request,
        _cellgeometry(wantsgeometry(request), grid))
end

"""
    precompute(geom)

Materialize `geom`'s per-cell geometry, trading O(ncells) residency for fewer
arithmetic operations per edge. Record semantics are identical — this only
changes where the centroids come from. Use it when the grid is known to fit in
memory; the default on-demand provider is what makes out-of-core grids possible.
Rectilinear payloads are already O(rows), so this is a no-op for them.
"""
precompute(geom::NeighborGeometry) = NeighborGeometry(geom.grid, geom.rings, geom.fields,
    _materialize(geom.payload, geom.grid))

_materialize(payload, grid) = payload
# `_materialize(::OnDemandCentroids, ::CellGrid)` — the one that actually
# fetches centroids — is in the DiscreteGlobalGrids extension.

# The from-cell centroid, fetched once per visited cell. The `OnDemandCentroids`
# methods below are the extension's; these two need no cell system at all.
@inline fromcentroid(::NoCentroids, cell, position) = nothing
@inline fromcentroid(p::StoredCentroids, cell, position) = @inbounds p.centroids[position]

# The per-edge distance and bearing. `NoCentroids` never touches DGG.
@inline edgegeometry(::NoCentroids, from, cell, position) = nothing
@inline edgegeometry(p::StoredCentroids, from, cell, position) =
    _edgegeometry(from, (@inbounds p.centroids[position]), p.radius)

# Compute great-circle distance and initial bearing between two unit-vector
# cell centroids.
@inline function _edgegeometry(from, to, R)
    central_angle = acos(clamp(sum(from[k] * to[k] for k in 1:3), -1.0, 1.0))
    lon1, lon2 = atan(from[2], from[1]), atan(to[2], to[1])
    lat1, lat2 = asin(clamp(from[3], -1.0, 1.0)), asin(clamp(to[3], -1.0, 1.0))
    delta_lon = lon2 - lon1
    east = sin(delta_lon) * cos(lat2)
    north = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(delta_lon)
    return (distance=central_angle * R, bearing=mod(rad2deg(atan(east, north)), 360.0))
end
