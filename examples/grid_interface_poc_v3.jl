# # Grid interface proof of concept — v3
#
# v3 implements phases 1–3 of `docs/src/grid-interface-perf-plan.md` on top of
# the design settled in `docs/src/grid-interface-decisions.md`. The grid
# objects, the façade boundary and the three drivers are v2's; what changed is
# how per-neighbor quantities are obtained.
#
# ## What changed from v2
#
# - **Request, don't cache.** A call states which per-neighbor quantities its
#   kernel reads — `needs=(Value(), Distance())` — as a tuple of singleton
#   `NeighborField` types. Because that is a concrete tuple type it specializes
#   through the function barriers, so the driver builds exactly the requested
#   fields with no runtime branching in the hot loop.
# - **Records hold exactly what was requested**, as a NamedTuple in the
#   canonical order `index, value, distance, bearing`. Reading an unrequested
#   field is an immediate field error rather than a silent slow path.
# - **No driver or geometry query allocates O(ncells).** A `CellGrid` with no
#   `Distance`/`Bearing` in the request does no centroid work at all; when
#   geometry is requested, centroids are computed on demand (one per visited
#   cell, one per edge) with zero residency. `precompute(geom)` is the explicit
#   opt-in to v2's whole-grid centroid table.
# - **Value-only kernels reach baseline speed by construction**: a request that
#   is a subset of `{Value}` is mapped onto `DGG.mapneighbors(...;
#   pass=DGG.Values())`, the streaming traversal Geomorphometry's own DGG
#   extension uses.
# - **NaN in, NaN out.** A nodata center cell yields NaN from `steepest_slope`,
#   `flow_direction` and `slope`; nodata neighbors are never selected. Integer
#   rasters are unaffected.
# - **Named algorithms over the drivers**: `topographic_position_index`,
#   `terrain_ruggedness_index` and `roughness`, with Geomorphometry's
#   definitions, on both backends and through the same façade.
# - `threaded` and `order` pass through the `CellGrid` driver to
#   `DGG.mapneighbors`.
# - **The sweep family** (phase 4): `neighbortable` hoists a position-keyed,
#   slot-stable adjacency for queue-driven traversal; `floodsweep` is one
#   priority flood behind a function barrier; `settle`, `flowdirection` and
#   `flowaccumulation` are O(n) passes over its result. This replaces v2's toy
#   sort-descending accumulation and reproduces Geomorphometry's priority-flood
#   D8 exactly, on both backends, from one implementation.

using Rasters
import DiscreteGlobalGrids as DGG
import Stencils
import GeoFormatTypes
using GeoFormatTypes: EPSG
using GeometryOpsCore: Manifold, Planar, Spherical, Geodesic
# The sweep family reuses Geomorphometry's direction encodings, its flow-method
# singletons and its QuickHeaps re-export. The import is qualified because both
# namespaces define `slope`, `roughness`, `flowaccumulation`, `Horn`, ....
import Geomorphometry as GM
using Geomorphometry: D8, DInf, FD8, FlowDirection, FlowDirectionConvention,
    FlowDirectionMethod, LDD, ispit

const AUTHALIC_RADIUS_M = 6_371_007.180918475

# ## Grid objects
#
# `P` maps the logical X and Y dimensions to storage axes. For example,
# `P == (2, 1)` means that X is stored on axis 2 and Y on axis 1. `lookups`
# contains the corresponding X and Y coordinates, regardless of storage order.
#
# The manifold selects planar or spherical geometry methods. The CRS is stored
# as metadata and does not participate in dispatch.

# This interface is independent of `DiscreteGlobalGrids.AbstractGrid`.
abstract type AbstractGridSpec end

struct RectilinearGrid{M<:Manifold,P,L,S,C} <: AbstractGridSpec
    manifold::M
    lookups::L
    spacing::S # Signed spacing in logical (X, Y) order
    crs::C     # Original CRS metadata; not used for dispatch
end

function RectilinearGrid(manifold::Manifold, P::Tuple{Int,Int}, lookups, spacing, crs)
    P in ((1, 2), (2, 1)) || throw(ArgumentError("invalid axis map $P"))
    all(!iszero, spacing) || throw(ArgumentError("cell spacing must be nonzero"))
    manifold isa Geodesic && throw(ArgumentError(
        "Geodesic manifolds are not supported yet; use Spherical (authalic radius)"))
    return RectilinearGrid{typeof(manifold),P,typeof(lookups),typeof(spacing),typeof(crs)}(
        manifold, lookups, spacing, crs)
end

axismap(::RectilinearGrid{M,P}) where {M,P} = P
yaxisnum(grid::RectilinearGrid) = axismap(grid)[2]

function gridsize(grid::RectilinearGrid)
    nx, ny = length(grid.lookups[1]), length(grid.lookups[2])
    return ntuple(axis -> axis == axismap(grid)[1] ? nx : ny, 2)
end

manifold(grid::RectilinearGrid) = grid.manifold

# A cell grid stores the cell collection and its DGG level grid so that DGG can
# provide topology and geometry. Its manifold is always `Spherical{Float64}`:
# the tessellation fixes the topology, but the radius scales distances and
# areas, so a `Spherical` with a custom radius may be supplied. The default is
# the authalic radius assumed by DGG's ISEA projection.

struct CellGrid{C,G} <: AbstractGridSpec
    manifold::Spherical{Float64}
    cells::C
    levelgrid::G
end

_cellmanifold(::Nothing) = Spherical(; radius=AUTHALIC_RADIUS_M)
_cellmanifold(m::Spherical) = Spherical(; radius=Float64(m.radius))
_cellmanifold(m::Manifold) = throw(ArgumentError(
    "a DGG cell grid is spherical; got $(typeof(m))"))

function CellGrid(celldim; manifold=nothing)
    lookup = Rasters.lookup(celldim)
    lookup isa DGG.CellLookup ||
        throw(ArgumentError("Dimension must contain a DGG.CellLookup; found $(typeof(lookup))"))
    cells = lookup.cells
    levelgrid = DGG.levelgrid(DGG.system(cells), DGG.level(cells))
    return CellGrid{typeof(cells),typeof(levelgrid)}(
        _cellmanifold(manifold), cells, levelgrid)
end

manifold(grid::CellGrid) = grid.manifold

# ## Adapters
#
# Matrices use axis 1 for X and axis 2 for Y. Without explicit spacing or a
# manifold, they use unit spacing and planar geometry.

function spatialparts(A::AbstractMatrix; spatialdims=nothing, spacing=nothing, manifold=nothing)
    isnothing(spatialdims) || spatialdims == (1, 2) ||
        throw(ArgumentError("a matrix uses spatialdims=(1, 2)"))
    steps = isnothing(spacing) ? (1.0, 1.0) : spacing
    m = isnothing(manifold) ? Planar() : manifold
    grid = RectilinearGrid(m, (1, 2), axes(A), steps, nothing)
    return A, grid
end

isspatialdim(::Type) = false
isspatialdim(::Type{<:Rasters.XDim}) = true
isspatialdim(::Type{<:Rasters.YDim}) = true
isspatialdim(::Type{<:DGG.Cells}) = true

_astuple(x::Tuple) = x
_astuple(x) = (x,)

# These provisional checks classify a Raster CRS once during grid construction.
# A complete implementation would use the CRS utilities in the Rasters
# extension.
_isgeographiccrs(::Nothing) = false
_isgeographiccrs(crs::EPSG) = GeoFormatTypes.val(crs) == 4326
_isgeographiccrs(crs) = occursin(r"GEOGCS|GEOGCRS"i, string(GeoFormatTypes.val(crs)))

function spatialparts(r::Raster; spatialdims=nothing, spacing=nothing, manifold=nothing)
    found_spatialdims = if isnothing(spatialdims)
        filter(d -> isspatialdim(typeof(d)), Rasters.dims(r))
    else
        Rasters.dims(r, spatialdims)
    end
    selected = _astuple(found_spatialdims)

    if length(selected) == 1 && lookup(only(selected)) isa DGG.CellLookup
        ndims(r) == 1 || throw(ArgumentError("Cells must be the only array axis"))
        isnothing(spacing) || throw(ArgumentError("CellGrid derives geometry from DGG"))
        return r, CellGrid(only(selected); manifold)
    end

    ndims(r) == 2 || throw(ArgumentError("a rectilinear surface must be 2D"))
    length(selected) == 2 ||
        throw(ArgumentError("expected one X and one Y dimension, or one Cells dimension"))

    xdim = only(filter(d -> d isa Rasters.XDim, selected))
    ydim = only(filter(d -> d isa Rasters.YDim, selected))
    P = Rasters.dimnum(r, (Rasters.XDim, Rasters.YDim))
    lookups = (Rasters.lookup(xdim), Rasters.lookup(ydim))
    steps = isnothing(spacing) ? (Float64(step(xdim)), Float64(step(ydim))) : spacing
    rascrs = Rasters.crs(r)
    m = if !isnothing(manifold)
        manifold
    elseif _isgeographiccrs(rascrs)
        Spherical(; radius=AUTHALIC_RADIUS_M)
    else
        Planar()
    end
    grid = RectilinearGrid(m, P, lookups, steps, rascrs)
    return r, grid
end

# ## The neighborhood vocabulary
#
# `NeighborRings(2)` includes both the first and second rings, matching the
# cumulative behavior of `DGG.neighbors(grid, cell, 2)`. The ring count is a
# runtime field so different values do not require separate compilations.

struct NeighborRings
    k::Int
    function NeighborRings(k::Integer=1)
        k >= 1 || throw(ArgumentError("neighbor ring count must be positive"))
        return new(Int(k))
    end
end

# Rectilinear offsets use logical `(dx, dy)` coordinates, where positive X is
# east and positive Y is north. Each square ring starts at north and proceeds
# counterclockwise. The first ring uses named directions; multiple rings use
# an ordered positional stencil.

const NORTH_UP_NEIGHBORS = (
    north=(0, 1), northwest=(-1, 1), west=(-1, 0), southwest=(-1, -1),
    south=(0, -1), southeast=(1, -1), east=(1, 0), northeast=(1, 1),
)

function logicalring(k)
    offsets = NTuple{2,Int}[]
    append!(offsets, ((x, k) for x in 0:-1:-k))
    append!(offsets, ((-k, y) for y in (k - 1):-1:-k))
    append!(offsets, ((x, -k) for x in (-k + 1):k))
    append!(offsets, ((k, y) for y in (-k + 1):k))
    append!(offsets, ((x, k) for x in (k - 1):-1:1))
    return Tuple(offsets)
end

logicaloffsets(rings::NeighborRings) =
    rings.k == 1 ? values(NORTH_UP_NEIGHBORS) :
    Tuple(Iterators.flatten(logicalring(k) for k in 1:rings.k))

function storageoffset(grid::RectilinearGrid, (dx, dy))
    xaxis, yaxis = axismap(grid)
    sx, sy = grid.spacing
    return ntuple(2) do axis
        axis == xaxis ? dx * Int(sign(sx)) : dy * Int(sign(sy))
    end
end

function northupstencil(grid::RectilinearGrid, rings::NeighborRings=NeighborRings())
    if rings.k == 1
        return Stencils.NamedStencil(map(o -> storageoffset(grid, o), NORTH_UP_NEIGHBORS))
    else
        return Stencils.Positional(map(o -> storageoffset(grid, o), logicaloffsets(rings)))
    end
end

# These methods expose stable cell identifiers and convert them to storage
# indices without requiring algorithms to know the grid representation.

cellindices(ras, ::RectilinearGrid) = CartesianIndices(axes(ras))

function cellindices(ras, grid::CellGrid)
    return (DGG.SubsetPositionedCell(grid.cells[p], p) for p in eachindex(grid.cells))
end

storagekey(::RectilinearGrid, I) = I
storagekey(::CellGrid, cell) = DGG.cellposition(cell)

# ## The request: which per-neighbor quantities does this kernel read?
#
# Every per-neighbor quantity is named by a singleton type under a shared
# supertype, and a call requests a tuple of them. The tuple's *type* carries
# the request, so it specializes through `Stencils.mapstencil` and
# `DGG.mapneighbors` — both of which are function barriers — and the driver
# never branches on the request at run time.

abstract type NeighborField end

"""Position of the neighbor: a `CartesianIndex` or a positioned DGG cell."""
struct Index <: NeighborField end
"""The neighbor's value in the array being mapped."""
struct Value <: NeighborField end
"""Distance from the center cell to the neighbor, in the manifold's units."""
struct Distance <: NeighborField end
"""Tangent-plane bearing to the neighbor, degrees clockwise from north."""
struct Bearing <: NeighborField end

# Records always list their fields in this order, whatever order the request
# used, so a kernel written against `(Value(), Distance())` also accepts records
# built from `(Distance(), Value())`.
const CANONICAL_FIELDS = (Index(), Value(), Distance(), Bearing())

# `Index` is deliberately *not* in the default: it is the field that forces the
# `CellGrid` driver off DGG's streaming value pass.
const DEFAULT_NEEDS = (Value(),)

# Request membership is answered with singleton types rather than `Bool`s, so
# every branch below is a dispatch that resolves at compile time even when
# constant propagation does not fire.
struct Requested end
struct Absent end

@inline _asked(::Type{F}, ::Tuple{}) where {F<:NeighborField} = Absent()
@inline _asked(::Type{F}, needs::Tuple) where {F<:NeighborField} =
    _askedhead(F, first(needs), Base.tail(needs))
@inline _askedhead(::Type{F}, ::F, rest::Tuple) where {F<:NeighborField} = Requested()
@inline _askedhead(::Type{F}, ::NeighborField, rest::Tuple) where {F<:NeighborField} =
    _asked(F, rest)

@inline _either(::Absent, ::Absent) = Absent()
@inline _either(::Requested, ::Any) = Requested()
@inline _either(::Absent, ::Requested) = Requested()

# `Distance` and `Bearing` share one geometry provider, so they are gated
# together.
@inline wantsgeometry(fields::Tuple) =
    _either(_asked(Distance, fields), _asked(Bearing, fields))

@inline _keepfield(::Requested, f, rest::Tuple) = (f, rest...)
@inline _keepfield(::Absent, f, rest::Tuple) = rest

@inline _canonical(needs::Tuple, ::Tuple{}) = ()
@inline _canonical(needs::Tuple, fields::Tuple) =
    _keepfield(_asked(typeof(first(fields)), needs), first(fields),
        _canonical(needs, Base.tail(fields)))

function requestfields(needs)
    needs isa Tuple || throw(ArgumentError(
        "needs must be a tuple of NeighborField singletons, got $(typeof(needs))"))
    all(n -> n isa NeighborField, needs) || throw(ArgumentError(
        "needs must be a tuple of NeighborField singletons, " *
        "e.g. (Value(), Distance()); got $needs"))
    return _canonical(needs, CANONICAL_FIELDS)
end

# One record field per requested quantity. `geo` carries the pair the geometry
# provider produced, or `nothing` when no geometry was requested — in which case
# the `Distance`/`Bearing` methods are unreachable by construction.
@inline _recordfield(::Index, index, value, geo) = (; index)
@inline _recordfield(::Value, index, value, geo) = (; value)
@inline _recordfield(::Distance, index, value, geo) = (; distance=geo.distance)
@inline _recordfield(::Bearing, index, value, geo) = (; bearing=geo.bearing)

@inline neighborrecord(::Tuple{}, index, value, geo) = NamedTuple()
@inline neighborrecord(fields::Tuple, index, value, geo) =
    merge(_recordfield(first(fields), index, value, geo),
        neighborrecord(Base.tail(fields), index, value, geo))

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
_materialize(p::OnDemandCentroids, grid::CellGrid) = StoredCentroids(p.levelgrid, p.radius,
    [DGG.cell_centroid(p.levelgrid, c) for c in grid.cells])

# The from-cell centroid, fetched once per visited cell.
@inline fromcentroid(::NoCentroids, cell, position) = nothing
@inline fromcentroid(p::OnDemandCentroids, cell, position) =
    DGG.cell_centroid(p.levelgrid, cell)
@inline fromcentroid(p::StoredCentroids, cell, position) = @inbounds p.centroids[position]

# The per-edge distance and bearing. `NoCentroids` never touches DGG.
@inline edgegeometry(::NoCentroids, from, cell, position) = nothing
@inline edgegeometry(p::OnDemandCentroids, from, cell, position) =
    _edgegeometry(from, DGG.cell_centroid(p.levelgrid, cell), p.radius)
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

# ## Neighbor iteration for traversal algorithms
#
# Queue-driven algorithms such as priority flood, accumulation, and HAND
# hoist the geometry once and reuse it while visiting cells:
#
#     geom = neighborgeometry(grid, NeighborRings(1), (Index(), Value(), Distance()))
#     visit!(out, geom, ras, queue) # a function barrier; see `_accumulate!`
#
# Both grid types yield records carrying exactly the requested fields.

@inline _neighborvalue(::Requested, ras, key) = ras[key]
@inline _neighborvalue(::Absent, ras, key) = nothing

function eachneighbor(geom::NeighborGeometry{<:RectilinearGrid}, ras, I::CartesianIndex)
    valid = CartesianIndices(ras)
    table = geometryat(geom, I)
    request = geom.fields
    wantvalue = _asked(Value, request)
    return (
        neighborrecord(request, I + o, _neighborvalue(wantvalue, ras, I + o), g)
        for (o, g) in zip(geom.payload.offsets, table) if (I + o) in valid
    )
end

function eachneighbor(geom::NeighborGeometry{<:CellGrid}, ras, cell)
    cells = geom.grid.cells
    payload = geom.payload
    request = geom.fields
    wantvalue = _asked(Value, request)
    p0 = DGG.cellposition(cell)
    from = fromcentroid(payload, cells[p0], p0)
    positions = DGG.neighbors(cells, p0, geom.rings.k)
    return (
        begin
            c = cells[p]
            nb = DGG.SubsetPositionedCell(c, p)
            neighborrecord(request, nb, _neighborvalue(wantvalue, ras, nb),
                edgegeometry(payload, from, c, p))
        end for p in positions
    )
end

# ## Cell areas
#
# `cellarea(grid)` returns cell areas in storage order. The return type reflects
# whether area is constant for the whole grid, constant within each latitude
# row, or different for every cell. The two lazy array types below provide the
# behavior needed here without depending on `FillArrays.Fill` in this example.

struct ConstantMatrix{T} <: AbstractMatrix{T}
    value::T
    dims::NTuple{2,Int}
end
Base.size(A::ConstantMatrix) = A.dims
Base.getindex(A::ConstantMatrix, ::Int, ::Int) = A.value
# Traversals hold their work arrays as flat vectors over storage positions, so
# these must answer a linear index directly rather than through a Cartesian
# conversion they do not need.
Base.IndexStyle(::Type{<:ConstantMatrix}) = Base.IndexLinear()
Base.getindex(A::ConstantMatrix, ::Int) = A.value
# Summing a constant grid does not require iterating over its cells.
Base.sum(A::ConstantMatrix) = A.value * prod(A.dims)

struct RowConstantMatrix{T,V<:AbstractVector{T}} <: AbstractMatrix{T}
    rows::V          # One value for each index along `rowaxis`
    rowaxis::Int     # Storage axis along which area varies
    dims::NTuple{2,Int}
end
Base.size(A::RowConstantMatrix) = A.dims
Base.getindex(A::RowConstantMatrix, i::Int, j::Int) = A.rows[A.rowaxis == 1 ? i : j]
Base.IndexStyle(::Type{<:RowConstantMatrix}) = Base.IndexLinear()
@inline function Base.getindex(A::RowConstantMatrix, p::Int)
    d1 = A.dims[1]
    return @inbounds A.rows[A.rowaxis == 1 ? mod1(p, d1) : (p - 1) ÷ d1 + 1]
end

cellarea(grid::RectilinearGrid{<:Planar}) =
    ConstantMatrix(abs(prod(grid.spacing)), gridsize(grid))

function cellarea(grid::RectilinearGrid{<:Spherical})
    R = grid.manifold.radius
    sx, sy = abs.(grid.spacing)
    areas = [
        R^2 * deg2rad(sx) * abs(sind(lat + sy / 2) - sind(lat - sy / 2))
        for lat in grid.lookups[2]
    ] # Exact spherical band area for each latitude row
    return RowConstantMatrix(areas, yaxisnum(grid), gridsize(grid))
end

_cellareaof(levelgrid, c, R) = DGG.cell_area(levelgrid, c) * R^2
# Z7 cells provide a closed-form area in square meters on the authalic sphere;
# rescale it when the grid's radius differs.
_cellareaof(levelgrid, c::DGG.Z7Cell, R) =
    DGG.IGeo7.cell_area(DGG.rawid(c)) * (R / AUTHALIC_RADIUS_M)^2

# `cellarea(grid)` on a cell grid is lazy. An IGeo7 area is closed form — one
# table read and a pentagon test — so materializing a vector buys arithmetic
# nobody was going to repeat and costs 8 bytes per cell of residency (123 MiB at
# level 13, and 625 ms to fill). Random access is a `CellVector` window lookup
# rather than O(1), so bulk consumers walk `p` ascending; `cellarea(grid, cell)`
# is the pointwise verb.
struct CellAreas{C,G} <: AbstractVector{Float64}
    cells::C
    levelgrid::G
    radius::Float64
end
Base.size(A::CellAreas) = (length(A.cells),)
Base.IndexStyle(::Type{<:CellAreas}) = Base.IndexLinear()
Base.@propagate_inbounds Base.getindex(A::CellAreas, p::Int) =
    _cellareaof(A.levelgrid, A.cells[p], A.radius)

cellarea(grid::CellGrid) = CellAreas(grid.cells, grid.levelgrid, manifold(grid).radius)

# Pointwise methods give traversal algorithms a common way to request the area
# of the current cell.
cellarea(grid::RectilinearGrid, I) = cellarea(grid)[I]
cellarea(grid::CellGrid, cell) =
    _cellareaof(grid.levelgrid, DGG.cellid(cell), manifold(grid).radius)

# ## `mapneighbors`: kernels that consume neighbor records
#
# The stencil boundary uses an isbits padding value, but neighbors whose indices
# fall outside the Raster are removed before the callback runs. The padding
# value therefore cannot reach the callback or widen its output type. The
# callback remains responsible for NaN or `missing` values inside the Raster.
#
# `Stencils.mapstencil` may run callbacks concurrently, so callbacks must not
# mutate shared state. It always threads through KernelAbstractions and offers
# no opt-out, which is why the rectilinear method below accepts no `threaded`
# keyword: passing one is an error rather than a silently ignored request.
#
# The three-argument form takes a `NeighborGeometry` directly, so a caller can
# hoist it, `precompute` it, and reuse it across calls.

mapneighbors(f, ras, grid::AbstractGridSpec, rings::NeighborRings=NeighborRings();
    needs=DEFAULT_NEEDS, kw...) =
    mapneighbors(f, ras, neighborgeometry(grid, rings, needs); kw...)

function mapneighbors(f, ras, geom::NeighborGeometry{<:RectilinearGrid})
    stencil = northupstencil(geom.grid, geom.rings)
    return _rectmapneighbors(f, ras, cellindices(ras, geom.grid), stencil,
        geom.fields, geom.payload.geometry)
end

# The stencil, stored offsets, and geometry table use the same logical offset
# order, so each zipped value, index, and geometry record describes the same
# neighbor. The `in valid` filter relies on `Stencils.indices` leaving
# out-of-bounds indices unclamped under `Remove`.
function _rectmapneighbors(f::F, ras, keys, stencil, request, geometry) where {F}
    valid = CartesianIndices(ras)
    # The value comes out of the stencil hood whether or not it was requested,
    # so there is nothing to gate here — only the record's field list.
    return Stencils.mapstencil(stencil, ras, keys;
        boundary=Stencils.Remove(zero(eltype(ras)))) do hood, I
        table = geometryat(geometry, I)
        neighbors = (
            neighborrecord(request, CartesianIndex(si), value, g)
            for (value, si, g) in
                zip(Stencils.neighbors(hood), Stencils.indices(stencil, I), table)
            if CartesianIndex(si) in valid
        )
        f(I, Stencils.center(hood), neighbors)
    end
end

# DGG performs cell traversal and infers the output type. A request that is a
# subset of `{Value}` rides DGG's streaming `Values()` pass — the same code path
# Geomorphometry's own DGG extension uses — so value-only kernels reach baseline
# speed by construction. Anything else needs the cell handles, and takes the
# `Neighbors()` pass.
function mapneighbors(f, ras::AbstractVector, geom::NeighborGeometry{<:CellGrid};
        threaded=true, order=DGG.StorageOrder())
    # DGG's traversal always hands the callback the one-ring, so a wider request
    # ignores it and re-queries the topology per cell.
    geom.rings.k == 1 || return _cellmapwiderings(f, ras, geom, threaded, order)
    return _cellmapneighbors(_streamable(geom.fields), f, ras, geom, threaded, order)
end

# Only a request with no `Index`, `Distance` or `Bearing` can stream values.
@inline _streamable(request::Tuple) = _streamable(_asked(Index, request),
    _asked(Distance, request), _asked(Bearing, request))
@inline _streamable(::Absent, ::Absent, ::Absent) = Requested()
@inline _streamable(::Any, ::Any, ::Any) = Absent()

function _cellmapneighbors(::Requested, f::F, ras, geom, threaded, order) where {F}
    request = geom.fields
    return DGG.mapneighbors(ras; pass=DGG.Values(), threaded, order) do cell, value, values
        f(cell, value,
            (neighborrecord(request, nothing, v, nothing) for v in values))
    end
end

# DGG already clipped the one-ring, so build records from the handles it hands
# us. Values are fetched per neighbor only when `Value` was requested.
function _cellmapneighbors(::Absent, f::F, ras, geom, threaded, order) where {F}
    return DGG.mapneighbors(ras; pass=DGG.Neighbors(), threaded, order) do cell, nbrs
        f(cell, ras[cell], _cellrecords(geom, ras, cell, nbrs))
    end
end

function _cellmapwiderings(f::F, ras, geom, threaded, order) where {F}
    return DGG.mapneighbors(ras; pass=DGG.Neighbors(), threaded, order) do cell, _
        f(cell, ras[cell], eachneighbor(geom, ras, cell))
    end
end

@inline function _cellrecords(geom, ras, cell, nbrs)
    payload = geom.payload
    request = geom.fields
    wantvalue = _asked(Value, request)
    from = fromcentroid(payload, DGG.cellid(cell), DGG.cellposition(cell))
    return (
        neighborrecord(request, nb, _neighborvalue(wantvalue, ras, nb),
            edgegeometry(payload, from, DGG.cellid(nb), DGG.cellposition(nb)))
        for nb in nbrs
    )
end

# ## `mapwindow`: named north-up windows for rectilinear grids
#
# Window kernels access neighbors by geographic direction, independent of
# storage order. For example, `w.northeast` always refers to the cell northeast
# of the center. This interface is not defined for hexagonal cell grids.

function mapwindow(f, ras, grid::RectilinearGrid)
    stencil = northupstencil(grid, NeighborRings(1))
    interior = CartesianIndices(map(r -> (first(r) + 1):(last(r) - 1), axes(ras)))
    return Stencils.mapstencil(stencil, ras, cellindices(ras, grid);
        boundary=Stencils.Remove(zero(eltype(ras)))) do hood, I
        I in interior || return NaN
        f(hood, metricspacing(grid, I))
    end
end

mapwindow(f, ras, ::CellGrid) = throw(ArgumentError(
    "windowed drivers require a rectilinear grid; use a neighbor-fit method (e.g. PlaneFit())"))

metricspacing(grid::RectilinearGrid{<:Planar}, I) = abs.(grid.spacing)

function metricspacing(grid::RectilinearGrid{<:Spherical}, I)
    R = grid.manifold.radius
    sx, sy = abs.(grid.spacing)
    lat = grid.lookups[2][I[yaxisnum(grid)]]
    return (deg2rad(sx) * cosd(lat) * R, deg2rad(sy) * R)
end

# ## Algorithms
#
# Each algorithm declares the fields it reads. Nothing else is computed for it.
#
# Nodata follows Geomorphometry's convention: NaN in, NaN out at that cell. The
# test is on the value, not on the grid, so integer rasters — which have no NaN
# — keep the identical code path with the predicate folded away.

@inline _isnodata(x::AbstractFloat) = isnan(x)
@inline _isnodata(x) = false

const STEEPEST_NEEDS = (Value(), Distance())
# `flow_direction` reports the steepest downhill bearing, so it reads one more
# field than `steepest_slope` and no index.
const DIRECTION_NEEDS = (Value(), Distance(), Bearing())
const PLANEFIT_NEEDS = (Value(), Distance(), Bearing())
# The sweep family's D8 rule is "lowest neighbor", not "steepest gradient", so it
# reads no geometry whatsoever — `neighbortable` serves it instead of the record
# API. This is what a multi-direction method would have to request, and why
# `DInf()` and `FD8()` are refused there rather than silently downgraded.
const MULTIDIRECTION_NEEDS = (Index(), Value(), Distance(), Bearing())
# Local statistics read values only, which is what puts them on DGG's streaming
# pass.
const LOCAL_NEEDS = (Value(),)

# A nodata center has no defined gradient; nodata neighbors are never selected.
function steepest(value, neighbors)
    _isnodata(value) && return (NaN, NaN)
    best_gradient = 0.0
    best_bearing = NaN
    for neighbor in neighbors
        _isnodata(neighbor.value) && continue
        gradient = (value - neighbor.value) / neighbor.distance
        if gradient > best_gradient
            best_gradient = gradient
            best_bearing = neighbor.bearing
        end
    end
    return best_gradient, best_bearing
end

# `steepest_slope` never reads a bearing, so it requests one field fewer than
# `flow_direction` and the driver skips the bearing arithmetic entirely.
function _steepest_slope(_, value, neighbors)
    _isnodata(value) && return NaN
    best_gradient = 0.0
    for neighbor in neighbors
        _isnodata(neighbor.value) && continue
        best_gradient = max(best_gradient, (value - neighbor.value) / neighbor.distance)
    end
    return atand(best_gradient)
end

function _flow_direction(_, value, neighbors)
    _, bearing = steepest(value, neighbors)
    return bearing # Clockwise from north; NaN indicates a pit, edge outlet or nodata
end

steepest_slope(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_steepest_slope, ras, grid, NeighborRings(); needs=STEEPEST_NEEDS, kw...)
flow_direction(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_flow_direction, ras, grid, NeighborRings(); needs=DIRECTION_NEEDS, kw...)

# Slope dispatches to an estimator chosen for the grid type. Horn uses a named
# rectilinear window. PlaneFit converts each neighbor's distance and bearing to
# local east and north coordinates, so it works with either grid type.

struct Horn end
struct PlaneFit end

slope(ras, grid::AbstractGridSpec; method=defaultmethod(slope, grid), kw...) =
    _slope(method, ras, grid; kw...)

defaultmethod(::typeof(slope), ::RectilinearGrid) = Horn()
defaultmethod(::typeof(slope), ::CellGrid) = PlaneFit()

function _horn_slope(w, (sx, sy))
    # Horn reads only the ring, so the center's nodata state must be tested
    # explicitly to keep the NaN-in/NaN-out convention.
    _isnodata(Stencils.center(w)) && return NaN
    dzdx = ((w.northeast + 2w.east + w.southeast) -
            (w.northwest + 2w.west + w.southwest)) / (8sx)
    dzdy = ((w.northeast + 2w.north + w.northwest) -
            (w.southeast + 2w.south + w.southwest)) / (8sy)
    return atand(hypot(dzdx, dzdy))
end

function _planefit_slope(_, value, neighbors)
    _isnodata(value) && return NaN
    see = sen = snn = sze = szn = 0.0
    for n in neighbors
        _isnodata(n.value) && continue
        east = n.distance * sind(n.bearing)
        north = n.distance * cosd(n.bearing)
        dz = n.value - value
        see += east * east; sen += east * north; snn += north * north
        sze += dz * east; szn += dz * north
    end
    det = see * snn - sen * sen
    iszero(det) && return NaN
    ddeast = (sze * snn - szn * sen) / det
    ddnorth = (szn * see - sze * sen) / det
    return atand(hypot(ddeast, ddnorth))
end

_slope(::Horn, ras, grid::RectilinearGrid) = mapwindow(_horn_slope, ras, grid)
_slope(::Horn, ras, grid::CellGrid) = throw(ArgumentError(
    "Horn() requires a rectilinear grid; use PlaneFit()"))
_slope(::PlaneFit, ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_planefit_slope, ras, grid, NeighborRings(); needs=PLANEFIT_NEEDS, kw...)

# ### Local statistics
#
# These match Geomorphometry's definitions (`src/relative.jl`) and read values
# only, so on a cell grid they run on DGG's streaming pass. Their nodata
# behavior is Geomorphometry's: NaN propagates through the arithmetic rather
# than being skipped, so a NaN neighbor makes the whole neighborhood NaN.

# TPI: the center minus the mean of its neighbors.
function _tpi_kernel(_, value, neighbors)
    total = 0.0
    count = 0
    for n in neighbors
        total += n.value
        count += 1
    end
    return oftype(float(value), value - total / count)
end

topographic_position_index(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_tpi_kernel, ras, grid, NeighborRings(); needs=LOCAL_NEEDS, kw...)
const TPI = topographic_position_index

# Roughness: the largest absolute difference between the center and a neighbor.
function _roughness_kernel(_, value, neighbors)
    o = zero(value)
    for n in neighbors
        o = max(o, abs(n.value - value))
    end
    return o
end

roughness(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_roughness_kernel, ras, grid, NeighborRings(); needs=LOCAL_NEEDS, kw...)

# TRI: the root of the summed squared differences to the neighbors. `normalize`
# divides by the neighbor count first; `squared=false` gives the mean absolute
# difference instead.
struct RuggednessKernel
    normalize::Bool
    squared::Bool
end

function (k::RuggednessKernel)(_, value, neighbors)
    total = 0.0
    count = 0
    for n in neighbors
        difference = abs(n.value - value)
        total += k.squared ? difference^2 : difference
        count += 1
    end
    k.normalize && (total /= count)
    return oftype(float(value), k.squared ? sqrt(total) : total)
end

function terrain_ruggedness_index(ras, grid::AbstractGridSpec;
        normalize=false, squared=true, kw...)
    if !normalize && !squared
        @warn "TRI: normalize=false and squared=false is not recommended."
    end
    return mapneighbors(RuggednessKernel(normalize, squared), ras, grid, NeighborRings();
        needs=LOCAL_NEEDS, kw...)
end
const TRI = terrain_ruggedness_index

# ## The sweep family
#
# A priority flood is random access in a data-dependent order, which is the one
# access pattern the record API is not built for: `neighborgeometry` serves
# whole-grid sweeps in storage order. So the sweep family gets its own hoisted
# structure — the traversal analogue of `neighborgeometry` — exposing the one
# thing the flood needs and the records do not carry: a *slot-stable,
# position-keyed* adjacency.
#
# Everything below is keyed by storage position (`1:ncells`) on both backends,
# never by cell id. The rectilinear and cell paths differ only in
# `neighbortable`, `_fillcellareas!` and the direction codec; the flood, the
# settling pass and the accumulation pass are one implementation each.

_data(A::AbstractArray) = A
_data(r::Raster) = parent(r)

# ### `neighbortable`: the traversal primitive
#
# `slots(tbl, p)` yields `(k, q)`: `k` is the ring slot — stable for the cell,
# not for the region — and `q` is the neighbor's storage position, or `0` when
# that ring member is outside the domain. `nslots(tbl, p)` is the complete
# degree. Slot order *is* `logicaloffsets(rings)` on a rectilinear grid and the
# complete ring order on a cell grid: slot ids, logical offsets and storage
# offsets are three parallel orderings of one list, so anyone reordering one
# reorders all three.
#
# Conceptually this is `needs = (Index(),)` specialized for random access.
# `Value` is served by indexing the value vector directly; `Distance` and
# `Bearing` are never requested, because the D8 rule is "lowest neighbor", not
# "steepest gradient" — distance does not enter it.

function neighbortable end

# On a rectilinear grid the table is arithmetic, so it holds no memory at all.
struct RectNeighborTable{K,C,L}
    cart::C
    lin::L
    offsets::NTuple{K,CartesianIndex{2}} # Storage offsets, in logical slot order
    logical::NTuple{K,NTuple{2,Int}}     # Logical (dx, dy): +x east, +y north
end

function neighbortable(grid::RectilinearGrid, rings::NeighborRings=NeighborRings(1))
    sz = gridsize(grid)
    logical = logicaloffsets(rings)
    offsets = map(o -> CartesianIndex(storageoffset(grid, o)), logical)
    return RectNeighborTable(CartesianIndices(sz), LinearIndices(sz), offsets, logical)
end

Base.length(t::RectNeighborTable) = length(t.lin)
nslots(::RectNeighborTable{K}, ::Int) where {K} = K

struct RectSlots{K,C,L}
    table::RectNeighborTable{K,C,L}
    I::CartesianIndex{2}
end
Base.length(::RectSlots{K}) where {K} = K
Base.eltype(::Type{<:RectSlots}) = Tuple{Int,Int}
@inline function Base.iterate(s::RectSlots{K}, k::Int=1) where {K}
    k > K && return nothing
    t = s.table
    J = s.I + @inbounds t.offsets[k]
    return ((k, J in t.cart ? @inbounds(t.lin[J]) : 0), k + 1)
end
@inline slots(t::RectNeighborTable, p::Int) = RectSlots(t, @inbounds t.cart[p])

# A cell grid takes one CSR table from DGG, built once per public call.
# `halo = :mark` rather than the clipped default buys two things: a row carries
# a `0` exactly where the cell has a complete-ring neighbor outside the region,
# which is the seed set for free, and complete-width rows preserve *slot
# identity*, so slot `k` is ring member `k` of the complete grid — which is what
# makes the direction codec a scan of a row already in cache.
struct CellNeighborTable{T}
    table::T
end

function neighbortable(grid::CellGrid, rings::NeighborRings=NeighborRings(1))
    rings.k == 1 ||
        throw(ArgumentError("the sweep family needs a one-ring neighbor table"))
    return CellNeighborTable(DGG.adjacency(grid.cells; halo=:mark, threaded=true))
end

Base.length(t::CellNeighborTable) = length(t.table)
@inline nslots(t::CellNeighborTable, p::Int) = length(@inbounds t.table[p])
@inline slots(t::CellNeighborTable, p::Int) = enumerate(@inbounds t.table[p])

# A cell is on the domain boundary when some ring member is outside the domain.
# On a rectilinear grid that is exactly the array frame; on a region of a cell
# grid it is exactly the rim. On a *complete* level grid it is empty — a sphere
# has no edges — which is the case `floodsweep` must not answer silently.
@inline function isboundary(tbl, p::Int)
    for (_, q) in slots(tbl, p)
        q == 0 && return true
    end
    return false
end

function _boundaryseeds(tbl, closedv)
    seeds = Int[]
    @inbounds for p in 1:length(tbl)
        (closedv[p] || !isboundary(tbl, p)) && continue
        push!(seeds, p)
    end
    return seeds
end

# ### The flood
#
# `order` is the pop sequence and `down[p]` the position that first reached `p`.
# On the surface the flood implies — the depression-filled one — `down[p]` is a
# downstream neighbor, which is Barnes' "flow directions from Priority-Flood":
# routing out of depressions is a property of the pop order, not of a separate
# filling pass.
struct FloodSweep{P<:Integer}
    order::Vector{P} # Pop sequence, truncated to the number of cells visited
    down::Vector{P}  # Downstream position; 0 at an outlet or an unvisited cell
    nseeds::Int
end

_closedmask(::Nothing, n) = falses(n)
function _closedmask(closed, n)
    mask = collect(Bool, vec(_data(closed)))
    length(mask) == n ||
        throw(ArgumentError("`closed` has $(length(mask)) entries for $n cells"))
    return mask
end

"""
    floodsweep(ras, grid; closed=nothing, seeds=nothing, table=neighbortable(grid))

Priority-flood traversal from the domain boundary inward, taking the lowest
queued value each step. Returns a `FloodSweep`.

`closed` cells are never enqueued, never visited and never receive or emit flow.
`seeds` overrides the boundary seed set with real hydrological outlets. `table`
is threaded through so a composite operation builds the adjacency once.

Serial by nature: the queue order is the algorithm. Allocates `order`, `down`, a
`closed` scratch and the queue, and nothing else that is O(ncells).
"""
function floodsweep(ras, grid::AbstractGridSpec; closed=nothing, seeds=nothing,
        table=neighbortable(grid, NeighborRings(1)))
    z = vec(_data(ras))
    n = length(z)
    closedv = _closedmask(closed, n)
    # Positions index the grid, so `Int32` carries any grid below 2.1e9 cells at
    # half the residency. The branch is answered once, above two barriers.
    return n <= typemax(Int32) ? _floodsweep(Int32, z, closedv, table, seeds) :
        _floodsweep(Int64, z, closedv, table, seeds)
end

function _floodsweep(::Type{P}, z, closedv, table, seeds) where {P}
    seedv = seeds === nothing ? _boundaryseeds(table, closedv) : collect(Int, seeds)
    isempty(seedv) && (seedv = _noboundaryseeds(z, closedv))
    order = Vector{P}(undef, length(z) - count(closedv))
    down = zeros(P, length(z))
    nvisited, nseeded = _floodsweep!(order, down, closedv, z, table, seedv)
    resize!(order, nvisited) # Unreachable cells leave a tail nobody may read
    return FloodSweep(order, down, nseeded)
end

# THE FUNCTION BARRIER. `table`'s concrete type is runtime-determined, so this
# loop cannot live at the hoist site: inlining a traversal there cost 35x time
# and ~170x allocations in the v2 measurements.
function _floodsweep!(order::Vector{P}, down::Vector{P}, closedv, z, tbl, seeds) where {P}
    # Array-backed, keyed by position: no hash and no `Dict` per enqueue.
    open = GM.FastPriorityQueue{eltype(z)}(length(z))
    nseeded = 0
    @inbounds for s in seeds
        p = Int(s)
        closedv[p] && continue # A seed the caller also closed is not a seed
        closedv[p] = true
        nseeded += 1
        GM.enqueue!(open, p, z[p])
    end
    nord = length(order)
    i = 0
    @inbounds while !isempty(open)
        p = GM.dequeue!(open)
        i += 1
        i <= nord || _sweepoverrun(i, nord)
        order[i] = P(p)
        for (_, q) in slots(tbl, p)
            (q == 0 || closedv[q]) && continue
            closedv[q] = true
            down[q] = P(p)
            GM.enqueue!(open, q, z[q])
        end
    end
    return i, nseeded
end

@noinline _sweepoverrun(i, n) = error(
    "the flood visited $i cells but `order` holds $n: the seed set and `closed` " *
    "disagree about how many cells are reachable")

# A complete level grid has no boundary, so there is no outlet to seed and the
# honest answers are "seed the lowest cell" or "refuse". Seeding the global
# minimum keeps the call total and conserves area; what it must not do is return
# the silent `acc == cellarea` answer an empty queue produces.
function _noboundaryseeds(z, closedv)
    isempty(z) && return Int[]
    minima = _globalminima(z, closedv)
    if isempty(minima)
        @warn "no seed cells: every cell is closed or non-finite. Nothing is visited."
    else
        @warn "no domain boundary: seeding the global minimum. A complete level grid " *
              "has no outlet — pass `seeds=` for a real one." nseeds = length(minima)
    end
    return minima
end

function _globalminima(z::AbstractVector{T}, closedv) where {T}
    found = false
    best = zero(T)
    @inbounds for p in eachindex(z)
        (closedv[p] || !isfinite(z[p])) && continue
        if !found || z[p] < best
            best = z[p]
            found = true
        end
    end
    found || return Int[]
    return [p for p in eachindex(z) if !closedv[p] && z[p] == best]
end

# ### `settle`: the depression-filled surface

"""
    settle(ras, grid; closed=nothing, sweep=floodsweep(ras, grid; closed))

Depression-filled ("settled") elevations: `settled[p] = max(z[p], settled[down[p]])`,
with outlets keeping their own elevation. This is the minimax (spill-point) fill,
so it is value-identical to Geomorphometry's `filldepressions` — flats are filled
to exactly the spill elevation, with no epsilon.

Non-finite input stays non-finite and does not propagate upstream; a cell the
flood never visited is `NaN`. Integer input is promoted to `float(eltype)`.
Serial: it walks the pop order forward.
"""
function settle(ras, grid::AbstractGridSpec; closed=nothing,
        sweep=floodsweep(ras, grid; closed))
    z = _floatvalues(vec(_data(ras)))
    out = similar(z)
    _settlepass!(out, z, sweep.order, sweep.down)
    return reshape(out, size(_data(ras)))
end

_floatvalues(z::AbstractVector{<:AbstractFloat}) = z
_floatvalues(z::AbstractVector) = float.(z) # So that "unvisited" is expressible

function _settlepass!(out, z, order, down)
    fill!(out, NaN)
    @inbounds for p in order # Forward: a parent always pops before its children
        zp = z[p]
        d = down[p]
        out[p] = (d != 0 && isfinite(zp) && isfinite(out[d])) ? max(zp, out[d]) : zp
    end
    return out
end

# ### Direction output
#
# `down` stays the internal representation; a direction code is produced in one
# final pass that is embarrassingly parallel. Rectilinear and IGeo7 grids get
# Geomorphometry's `LDD` numpad, which makes the output bit-comparable with it.

"""
Ring-slot flow directions: the value is the neighbor's slot in the cell's
complete ring, and `0` is a pit. Cell systems with no relative-cell arithmetic
use this, because an LDD numpad code has no meaning on an arbitrary ring.
"""
struct RingSlot <: FlowDirectionConvention end

GM._arrow(::Type{RingSlot}, d::Integer) =
    iszero(d) ? '·' : Char(0x2460 + (Int(d) - 1)) # ①②③…
GM.ispit(d::FlowDirection{RingSlot}) = iszero(d.value)

_directiontype(::RectilinearGrid) = FlowDirection{LDD,UInt8}
_directiontype(grid::CellGrid) = _celldirectiontype(eltype(grid.cells))
_celldirectiontype(::Type{<:DGG.Z7Cell}) = FlowDirection{LDD,UInt8}
_celldirectiontype(::Type) = FlowDirection{RingSlot,UInt8}

"""
    flowdirection(ras, grid; closed=nothing, table=…, sweep=…)

The downstream direction of each cell on the settled surface, as
`FlowDirection{LDD}` on rectilinear and IGeo7 grids and `FlowDirection{RingSlot}`
on other cell grids. Outlets, and cells the flood never visited, are pits.

A cell whose value is nodata is *not* a pit: the flood visits it last and gives
it a parent, matching Geomorphometry. Pass `closed` to exclude it instead.
"""
function flowdirection(ras, grid::AbstractGridSpec; closed=nothing,
        table=neighbortable(grid, NeighborRings(1)),
        sweep=floodsweep(ras, grid; closed, table))
    out = Array{_directiontype(grid)}(undef, size(_data(ras))) # Rule A: plain Array
    _encodedirections!(vec(out), sweep.down, grid, table)
    return out
end

# The rectilinear slot already *is* the logical direction, so the numpad code is
# read straight off `logicaloffsets`. No `cellsize` sign logic is needed and none
# appears: storage order was handled once, in `storageoffset`, and both storage
# orders are therefore correct by construction rather than by a compensating flip.
_lddcode((dx, dy)) = -1 <= dx <= 1 && -1 <= dy <= 1 ?
    GM._ldd_ci2dir[CartesianIndex(dx, dy)] :
    throw(ArgumentError("the LDD numpad encodes one ring; got offset ($dx, $dy)"))

_encodedirections!(out, down, ::RectilinearGrid, t::RectNeighborTable{K}) where {K} =
    _rectencode!(out, down, t, ntuple(k -> _lddcode(t.logical[k]), Val(K)), _lddcode((0, 0)))

function _rectencode!(out, down, t, codes::NTuple{K,UInt8}, pit) where {K}
    Threads.@threads for r in _chunkranges(length(out))
        @inbounds for p in r
            d = Int(down[p])
            code = pit
            if d != 0
                for (k, q) in slots(t, p)
                    q == d && (code = codes[k]; break)
                end
            end
            out[p] = FlowDirection{LDD,UInt8}(code)
        end
    end
    return out
end

# The `:mark` row is complete-width, so slot `k` is ring member `k` of the
# complete grid — and on IGeo7 that is exactly `DGG.directioncode`, checked
# below. The codec is therefore a scan of a CSR row already in cache: no ring
# rebuild (Geomorphometry's `_ringslot`), no `CellVector` window lookup and no
# relative-cell arithmetic.
#
# Table from `ext/GeomorphometryDiscreteGlobalGridsExt.jl`'s `IGEO7_TO_LDD`:
# ·, E, NE, NW, W, SW, SE. The numpad's N(8) and S(2) are unused on a hex ring.
const IGEO7_TO_LDD = (0x05, 0x06, 0x09, 0x07, 0x04, 0x01, 0x03)

_encodedirections!(out, down, grid::CellGrid, t::CellNeighborTable) =
    _cellencode!(out, down, t, _celldirectiontype(eltype(grid.cells)))

function _cellencode!(out, down, t, ::Type{FlowDirection{LDD,UInt8}})
    Threads.@threads for r in _chunkranges(length(out))
        @inbounds for p in r
            out[p] = FlowDirection{LDD,UInt8}(IGEO7_TO_LDD[_downslot(t, p, Int(down[p])) + 1])
        end
    end
    return out
end

function _cellencode!(out, down, t, ::Type{FlowDirection{RingSlot,UInt8}})
    Threads.@threads for r in _chunkranges(length(out))
        @inbounds for p in r
            out[p] = FlowDirection{RingSlot,UInt8}(_downslot(t, p, Int(down[p])))
        end
    end
    return out
end

@inline function _downslot(tbl, p::Int, d::Int)
    d == 0 && return 0
    for (k, q) in slots(tbl, p)
        q == d && return k
    end
    return 0
end

"""
    downstreamposition(table, p, direction::FlowDirection{RingSlot}) -> Int

Decode a ring-slot direction back to a storage position, or `0` at a pit. O(1)
off the row the encoder read, which is what makes `RingSlot` cheap in both
directions on a grid with no relative-cell arithmetic.
"""
@inline downstreamposition(tbl::CellNeighborTable, p::Int, d::FlowDirection{RingSlot}) =
    ispit(d) ? 0 : (@inbounds tbl.table[p][Int(d)])

# ### Accumulation

"""
    flowaccumulation(ras, grid; method=D8(), closed=nothing) -> (acc, directions)

`acc[c]` is the total upstream area draining through `c`, in the grid's area
units, including `c`'s own cell area. `directions` is the matching
[`flowdirection`](@ref) output.

`acc` is `Float64`: at level 13 the largest accumulation is ~2.7e9 m² built from
16.2M additions of ~526 m² each, and a `Float32` mantissa resolves only ~161 m²
there — a third of a cell lost per addition near the outlet.

The flood and the downhill pass are serial; the cell-area fill and the direction
encoding are threaded over contiguous chunks.
"""
flowaccumulation(ras, grid::AbstractGridSpec; method=D8(), closed=nothing) =
    _flowaccumulation(method, ras, grid, closed)

function _flowaccumulation(::D8, ras, grid::AbstractGridSpec, closed)
    table = neighbortable(grid, NeighborRings(1)) # Built once, threaded through
    sweep = floodsweep(ras, grid; closed, table)
    acc = Vector{Float64}(undef, length(_data(ras)))
    _fillcellareas!(acc, grid)
    _accumulatedown!(acc, sweep.order, sweep.down)
    dirs = flowdirection(ras, grid; sweep, table)
    return reshape(acc, size(_data(ras))), dirs
end

_flowaccumulation(method::FlowDirectionMethod, ras, grid, closed) = throw(ArgumentError(
    "flowaccumulation with $(nameof(typeof(method))) splits flow between several " *
    "neighbors, so it needs per-neighbor bearings and distances — " *
    "needs = (Index(), Value(), Distance(), Bearing()) — and a ragged sweep result. " *
    "The sweep family carries one downstream position per cell. Use method = D8()."))

# Reverse pop order is a topological order: a child always pops after its parent,
# so reversing guarantees a cell's own inflow is complete before it drains.
function _accumulatedown!(acc, order, down)
    @inbounds for j in reverse(eachindex(order))
        p = order[j]
        d = down[p]
        d == 0 && continue
        acc[d] += acc[p]
    end
    return acc
end

# `cellarea` supplies the only unit-carrying input to the whole family. Nothing
# else in it needs geometry: no distances, no bearings, no manifold dispatch.
_fillcellareas!(acc, grid::RectilinearGrid) = _copychunks!(acc, cellarea(grid))
_fillcellareas!(acc, grid::CellGrid) =
    _cellareachunks!(acc, grid.cells, grid.levelgrid, manifold(grid).radius)

function _copychunks!(acc, areas)
    Threads.@threads for r in _chunkranges(length(acc))
        @inbounds for p in r
            acc[p] = areas[p]
        end
    end
    return acc
end

# Each chunk walks `p` ascending, which is what a `CellVector` window cursor wants.
function _cellareachunks!(acc, cells, levelgrid, radius)
    Threads.@threads for r in _chunkranges(length(acc))
        @inbounds for p in r
            acc[p] = _cellareaof(levelgrid, cells[p], radius)
        end
    end
    return acc
end

function _chunkranges(n::Int, nchunks::Int=Threads.nthreads())
    n <= 0 && return UnitRange{Int}[]
    width = cld(n, max(1, min(nchunks, n)))
    return [i:min(i + width - 1, n) for i in 1:width:n]
end

# ## The façade boundary
#
# This signature defines all grid-construction keywords. Other keywords are
# forwarded to the algorithm method, where Julia reports unsupported keywords
# normally. `spatialdims`, `spacing`, and `manifold` are reserved by the public
# API.
#
# `needs` is deliberately *not* a façade keyword: every named algorithm knows
# its own request. Users writing their own kernels call `mapneighbors` or
# `eachneighbor` and pass `needs` explicitly.

splitspatial(; spatialdims=nothing, spacing=nothing, manifold=nothing, kwargs...) =
    (; spatialdims, spacing, manifold), kwargs

# Algorithms return data without rebuilding input metadata. These methods
# preserve plain-array outputs and restore Raster metadata, including the output
# name and floating-point missing-value convention.

rebuildoutput(input, grid, data; name) = data

# One `rebuild`, not two. An intermediate Raster whose `missingval` does not
# match its new eltype makes Rasters cook one up, and cooking one up is
# `typemin`/`typemax` on the eltype — which `FlowDirection` does not define.
# This is the single place the missingval machinery is reached, which is what
# structurally disarms that crash rather than working around it per call site.
function rebuildoutput(input::Raster, grid, data; name)
    mv = eltype(data) <: AbstractFloat ? eltype(data)(NaN) : nothing
    return data isa Raster ? Rasters.rebuild(data; name, missingval=mv) :
        Rasters.rebuild(input; data, name, missingval=mv)
end

# An algorithm may return several products; the façade names them positionally.
rebuildoutput(input, grid, data::Tuple, names::Tuple) =
    map((d, nm) -> rebuildoutput(input, grid, d; name=nm), data, names)
rebuildoutput(input, grid, data, names::Tuple) =
    rebuildoutput(input, grid, data; name=first(names))

const OUTPUTNAMES = Dict(
    :flowaccumulation => (:flowaccumulation, :flowdirection),
    :flowdirection => (:flowdirection,),
    :settle => (:settle,),
)

for f in (:steepest_slope, :flow_direction, :slope, :settle, :flowdirection,
          :flowaccumulation, :topographic_position_index,
          :terrain_ruggedness_index, :roughness)
    names = get(OUTPUTNAMES, f, (f,))
    @eval function $f(input; kwargs...)
        spatial, rest = splitspatial(; kwargs...)
        ras, grid = spatialparts(input; spatial...)
        data = $f(ras, grid; rest...)
        return rebuildoutput(input, grid, data, $(QuoteNode(names)))
    end
end

# ## Behavior checks
#
# Store the same plane in X,Y and Y,X axis order. The reversed Y coordinates and
# unequal X and Y spacing expose calculations that incorrectly depend on
# storage-axis order.

xs = 0.0:2.0:8.0
ys = 20.0:-5.0:5.0
surface(x, y) = 100.0 - 3.0x

data_xy = [surface(x, y) for x in xs, y in ys]
raster_xy = Raster(data_xy, (X(xs), Y(ys)))
raster_yx = Raster(permutedims(data_xy), (Y(ys), X(xs)))

_, grid_xy = spatialparts(raster_xy)
_, grid_yx = spatialparts(raster_yx)
@assert axismap(grid_xy) == (1, 2)
@assert axismap(grid_yx) == (2, 1)
# A Raster without a CRS defaults to planar geometry.
@assert grid_xy.manifold isa Planar
@assert collect(cellindices(raster_xy, grid_xy)) == collect(CartesianIndices(raster_xy))

# ### The request API
#
# The canonical field order is independent of the order the request used, and
# a record carries exactly the requested fields — no more.
@assert requestfields((Value(),)) === (Value(),)
@assert requestfields((Bearing(), Value(), Index())) === (Index(), Value(), Bearing())
@assert requestfields((Value(), Value())) === (Value(),)
@assert requestfields(()) === ()
bad_needs = try requestfields((:value,)); nothing catch e; e end
@assert bad_needs isa ArgumentError
bad_needs2 = try requestfields(Value()); nothing catch e; e end
@assert bad_needs2 isa ArgumentError

record = neighborrecord(requestfields((Distance(), Value())), CartesianIndex(1, 1),
    3.0, (distance=2.0, bearing=90.0))
@assert record === (value=3.0, distance=2.0)
@assert keys(record) == (:value, :distance)
# (a) Reading a field that was not requested is a loud failure: Julia's own
# NamedTuple field error, naming the field and listing what is available.
missing_field = try record.bearing; nothing catch e; e end
@assert missing_field isa Exception
@assert occursin("bearing", sprint(showerror, missing_field))
missing_index = try record.index; nothing catch e; e end
@assert missing_index isa Exception
@assert occursin("index", sprint(showerror, missing_index))
# The same failure reaches a kernel: PlaneFit's request has no `index`.
kernel_field_error = try
    mapneighbors(raster_xy, grid_xy; needs=PLANEFIT_NEEDS) do _, value, neighbors
        sum(n -> n.index[1], neighbors)
    end
    nothing
catch e
    e
end
@assert kernel_field_error isa Exception

# Named stencil fields refer to the same geographic directions in both storage
# layouts.
tagged(x, y) = 100.0x + y
tagged_xy = Raster([tagged(x, y) for x in xs, y in ys], (X(xs), Y(ys)))
tagged_yx = Raster(permutedims(parent(tagged_xy)), (Y(ys), X(xs)))
stencil_xy, stencil_yx = northupstencil(grid_xy), northupstencil(grid_yx)
hood_xy = Stencils.stencil(Stencils.StencilArray(tagged_xy, stencil_xy), (3, 2))
hood_yx = Stencils.stencil(Stencils.StencilArray(tagged_yx, stencil_yx), (2, 3))
@assert all(name -> getproperty(hood_xy, name) == getproperty(hood_yx, name),
    keys(NORTH_UP_NEIGHBORS))
@assert hood_xy.north == tagged(xs[3], ys[1])
@assert hood_xy.east == tagged(xs[4], ys[2])

# Increasing the ring count expands the neighborhood without changing the
# callback. Index filtering removes out-of-bounds neighbors.
neighborcount(_, _, neighbors) = count(_ -> true, neighbors)
counts1_xy = mapneighbors(neighborcount, raster_xy, grid_xy)
@assert eltype(counts1_xy) == Int # Boundary handling does not add `missing`
@assert parent(counts1_xy)[1, 1] == 3
counts2_xy = mapneighbors(neighborcount, raster_xy, grid_xy, NeighborRings(2))
counts2_yx = mapneighbors(neighborcount, raster_yx, grid_yx, NeighborRings(2))
@assert all(parent(counts2_xy) .>= parent(counts1_xy))
@assert maximum(parent(counts2_xy)) > maximum(parent(counts1_xy))
@assert parent(counts2_xy) == permutedims(parent(counts2_yx))

# A rectilinear geometry table is only built when distance or bearing is asked
# for; the default request leaves the payload a zero-size singleton.
@assert neighborgeometry(grid_xy).payload.geometry isa NoGeometry
let withgeom = neighborgeometry(grid_xy, NeighborRings(1), (Value(), Distance()))
    @assert withgeom.payload.geometry isa UniformGeometry
end

slope_xy = steepest_slope(raster_xy)
slope_yx = steepest_slope(raster_yx; spatialdims=(X, Y))
direction_xy = flow_direction(raster_xy)
direction_yx = flow_direction(raster_yx)

@assert slope_xy isa Raster
# Boundary handling preserves the inferred output type.
@assert eltype(slope_xy) == Float64
# The public method restores output metadata.
@assert Rasters.name(slope_xy) == :steepest_slope
@assert parent(slope_xy) ≈ permutedims(parent(slope_yx))
@assert all(isequal.(parent(direction_xy), permutedims(parent(direction_yx))))
@assert all(d -> isnan(d) || d ≈ 90.0, parent(direction_xy))

# Plain matrices use the same algorithms through an inferred rectilinear grid.
# An isbits padding value also supports integer input arrays.
matrix_slope = steepest_slope(data_xy; spacing=(2.0, -5.0))
@assert matrix_slope isa Matrix{Float64} # Plain arrays produce plain-array outputs
@assert matrix_slope ≈ parent(slope_xy)
int_slope = steepest_slope(Int.(data_xy); spacing=(2.0, -5.0))
@assert int_slope ≈ matrix_slope

# Grid-construction keywords are consumed by `spatialparts`. Unknown algorithm
# keywords raise a `MethodError`, and unsupported manifolds raise an
# `ArgumentError` during grid construction.
bad_kwarg = try steepest_slope(raster_xy; bogus=1); nothing catch e; e end
@assert bad_kwarg isa MethodError
geodesic_err = try spatialparts(raster_xy; manifold=Geodesic()); nothing catch e; e end
@assert geodesic_err isa ArgumentError
# `Stencils.mapstencil` always threads, so the rectilinear driver has no
# `threaded` knob and says so by refusing the keyword.
no_threading = try steepest_slope(raster_xy; threaded=false); nothing catch e; e end
@assert no_threading isa MethodError

# Horn and PlaneFit recover the same slope for an exact plane through their
# respective window and neighbor-record interfaces.
horn = slope(raster_xy) # Rectilinear grids use Horn by default
@assert all(v -> v ≈ atand(3.0), parent(horn)[2:end-1, 2:end-1])
# A complete Horn window is unavailable at the border.
@assert all(isnan, parent(horn)[1, :])
pf = slope(raster_xy; method=PlaneFit())
@assert all(v -> v ≈ atand(3.0), parent(pf)) # Plane fitting uses available edge neighbors
pf_yx = slope(raster_yx; spatialdims=(X, Y), method=PlaneFit())
@assert parent(pf) ≈ permutedims(parent(pf_yx))

# Every cell in this planar grid has the same area.
ca = cellarea(grid_xy)
@assert ca[3, 2] == 10.0
@assert sum(ca) == 10.0 * length(raster_xy)
@assert cellarea(grid_xy, CartesianIndex(1, 1)) == 10.0

# Flow accumulation is the priority-flood sweep, and it returns two named
# products. All area reaches a cell with no downstream neighbor, which is an
# outlet; `missingval = nothing` on the direction raster is what keeps Rasters
# from cooking one up out of `typemax(FlowDirection)`.
acc, dirs = flowaccumulation(raster_xy)
@assert Rasters.name(acc) == :flowaccumulation
@assert Rasters.name(dirs) == :flowdirection
@assert eltype(dirs) == FlowDirection{LDD,UInt8}
@assert isnothing(Rasters.missingval(dirs))
@assert isnan(Rasters.missingval(acc))
sweep_xy = floodsweep(raster_xy, grid_xy)
@assert sum(parent(acc)[sweep_xy.down .== 0]) ≈ sum(cellarea(grid_xy))
@assert all(parent(acc) .>= 10.0) # Every cell carries at least its own area
@assert count(ispit, parent(dirs)) == sweep_xy.nseeds

# ### (e) Local statistics against hand-computed values
#
# `bump` is flat except for one raised cell, so every neighborhood statistic has
# an obvious closed form.
bump = zeros(5, 5)
bump[3, 3] = 10.0
bump_ras = Raster(bump, (X(1.0:5.0), Y(1.0:5.0)))
tpi_bump = topographic_position_index(bump_ras)
tri_bump = terrain_ruggedness_index(bump_ras)
rough_bump = roughness(bump_ras)
@assert Rasters.name(tpi_bump) == :topographic_position_index
# The peak sits 10 above the mean of its eight zero neighbors.
@assert parent(tpi_bump)[3, 3] ≈ 10.0
# A cell orthogonally adjacent to the peak has 8 neighbors, one of them the peak.
@assert parent(tpi_bump)[3, 2] ≈ -10.0 / 8
# A corner cell has only 3 neighbors, none of them the peak.
@assert parent(tpi_bump)[1, 1] ≈ 0.0
@assert parent(tri_bump)[3, 3] ≈ sqrt(8 * 10.0^2)
@assert parent(tri_bump)[3, 2] ≈ 10.0
@assert parent(tri_bump)[1, 1] ≈ 0.0
@assert parent(terrain_ruggedness_index(bump_ras; normalize=true))[3, 3] ≈ sqrt(8 * 100 / 8)
@assert parent(terrain_ruggedness_index(bump_ras; normalize=true, squared=false))[3, 3] ≈ 10.0
@assert parent(rough_bump)[3, 3] ≈ 10.0
@assert parent(rough_bump)[3, 2] ≈ 10.0
@assert parent(rough_bump)[1, 1] ≈ 0.0
# Storage order does not change a local statistic.
bump_yx = Raster(permutedims(bump), (Y(1.0:5.0), X(1.0:5.0)))
@assert parent(topographic_position_index(bump_yx)) ≈ permutedims(parent(tpi_bump))
@assert parent(roughness(bump_yx)) ≈ permutedims(parent(rough_bump))

# ### (d) Nodata semantics
#
# A NaN center yields NaN; a NaN neighbor is never selected; every other cell is
# untouched.
nan_data = copy(data_xy)
nan_data[3, 2] = NaN
nan_ras = Raster(nan_data, (X(xs), Y(ys)))
nan_slope = steepest_slope(nan_ras)
nan_dir = flow_direction(nan_ras)
nan_pf = slope(nan_ras; method=PlaneFit())
nan_horn = slope(nan_ras)
@assert isnan(parent(nan_slope)[3, 2])
@assert isnan(parent(nan_dir)[3, 2])
@assert isnan(parent(nan_pf)[3, 2])
@assert isnan(parent(nan_horn)[3, 2])
# A nodata neighbor is skipped rather than selected, and never poisons the
# result. The cell west of the hole would have flowed straight east into it; it
# now reports the steepest *available* descent, which is the southeast diagonal.
@assert isfinite(parent(nan_slope)[2, 2])
@assert parent(nan_slope)[2, 2] ≈ atand(6.0 / hypot(2.0, 5.0))
@assert parent(nan_dir)[2, 2] ≈ mod(atand(2.0, -5.0), 360.0)
# The NaN cell's eastern neighbor would have flowed east anyway; the point is
# that it never picks the NaN cell, whose gradient comparison is skipped.
@assert parent(nan_slope)[4, 2] ≈ parent(slope_xy)[4, 2]
@assert parent(nan_dir)[4, 2] ≈ 90.0
# Cells that do not touch the hole are bit-identical to the NaN-free run.
for column in (1, 5)
    @assert all(isequal.(parent(nan_slope)[column, :], parent(slope_xy)[column, :]))
    @assert all(isequal.(parent(nan_dir)[column, :], parent(direction_xy)[column, :]))
    @assert all(isequal.(parent(nan_pf)[column, :], parent(pf)[column, :]))
end
# Integer rasters have no NaN and keep working through the same kernels.
@assert steepest_slope(Int.(data_xy); spacing=(2.0, -5.0)) ≈ matrix_slope
@assert eltype(roughness(Int.(data_xy))) == Int

# An EPSG:4326 Raster uses spherical geometry. East-west distances and cell
# areas vary by latitude row.
lons = 0.0:1.0:5.0
lats = 60.0:-1.0:55.0
geo_data = [1000.0 - 10.0 * lon for lon in lons, lat in lats]
geo_raster = Raster(geo_data, (X(lons), Y(lats)); crs=EPSG(4326))
_, grid_geo = spatialparts(geo_raster)
@assert grid_geo.manifold isa Spherical

geom_geo = neighborgeometry(grid_geo, NeighborRings(1), (Value(), Distance(), Bearing()))
@assert geom_geo.payload.geometry isa RowGeometry
east = findfirst(==(:east), collect(keys(NORTH_UP_NEIGHBORS)))
t60 = geometryat(geom_geo, CartesianIndex(3, 1))
t55 = geometryat(geom_geo, CartesianIndex(3, 6))
@assert t60[east].bearing == 90.0
@assert t60[east].distance ≈ deg2rad(1.0) * cosd(60.0) * AUTHALIC_RADIUS_M
@assert t55[east].distance > t60[east].distance # Longitude spacing widens toward the equator

slope_geo = steepest_slope(geo_raster)
@assert parent(slope_geo)[3, 1] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(60.0) * AUTHALIC_RADIUS_M))
@assert parent(slope_geo)[3, 6] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(55.0) * AUTHALIC_RADIUS_M))
horn_geo = slope(geo_raster)
@assert parent(horn_geo)[3, 2] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(59.0) * AUTHALIC_RADIUS_M))

ca_geo = cellarea(grid_geo)
@assert ca_geo[1, 1] ≈ AUTHALIC_RADIUS_M^2 * deg2rad(1.0) * (sind(60.5) - sind(59.5))
@assert ca_geo[1, 6] > ca_geo[1, 1]

# A DGG Raster uses the same public algorithms through `CellGrid`. Elevation is
# proportional to each centroid's Z coordinate, so most cells have a downhill
# neighbor.
level = 2
dgg = DGG.levelgrid(DGG.IGeo7System(), level)
cell_lookup = DGG.CellLookup(dgg)
cell_values = [10_000.0 * DGG.cell_centroid(dgg, cell)[3] for cell in cell_lookup]
cell_raster = Raster(cell_values, (DGG.Cells(cell_lookup),))

_, cell_grid = spatialparts(cell_raster; spatialdims=DGG.Cells)
@assert cell_grid isa CellGrid
@assert cell_grid.cells === cell_lookup.cells
@assert manifold(cell_grid) isa Spherical
@assert manifold(cell_grid).radius == AUTHALIC_RADIUS_M

# A user-provided radius rescales every derived quantity; on a unit sphere the
# cell areas sum to the full sphere's solid angle. Non-spherical manifolds are
# rejected.
_, unit_grid = spatialparts(cell_raster; manifold=Spherical(; radius=1.0))
@assert manifold(unit_grid).radius == 1.0
@assert sum(cellarea(unit_grid)) ≈ 4π
planar_cells_err = try spatialparts(cell_raster; manifold=Planar()); nothing catch e; e end
@assert planar_cells_err isa ArgumentError

cell_index = first(cellindices(cell_raster, cell_grid))
@assert cell_index isa DGG.SubsetPositionedCell
@assert cell_raster[cell_index] == first(cell_values)

cell_counts2 = mapneighbors(neighborcount, cell_raster, cell_grid, NeighborRings(2))
@assert parent(cell_counts2)[1] == length(DGG.neighbors(cell_lookup, 1, 2))

cell_areas = cellarea(cell_grid)
@assert sum(cell_areas) ≈ 4π * AUTHALIC_RADIUS_M^2 # IGeo7 cells cover the sphere
@assert cellarea(cell_grid, cell_index) == cell_areas[1]

cell_slope = steepest_slope(cell_raster)
cell_direction = flow_direction(cell_raster)
@assert cell_slope isa Raster
@assert all(>=(0.0), parent(cell_slope))
@assert any(isfinite, parent(cell_direction))
@assert all(d -> isnan(d) || 0.0 <= d < 360.0, parent(cell_direction))

pf_cells = slope(cell_raster) # Cell grids use PlaneFit by default
@assert all(isfinite, parent(pf_cells))
horn_err = try slope(cell_raster; method=Horn()); nothing catch e; e end
@assert horn_err isa ArgumentError

# The same local statistics run on the cell backend, through DGG's streaming
# value pass, and reproduce the same hand-checkable relationships.
cell_tpi = topographic_position_index(cell_raster)
cell_tri = terrain_ruggedness_index(cell_raster)
cell_rough = roughness(cell_raster)
@assert cell_tpi isa Raster && length(cell_tpi) == length(cell_raster)
@assert Rasters.name(cell_rough) == :roughness
@assert all(>=(0.0), parent(cell_rough))
@assert all(isfinite, parent(cell_tpi))
# TRI with no normalization is the root of the summed squares of the same
# differences roughness takes the maximum of, so it is never the smaller one.
@assert all(parent(cell_tri) .>= parent(cell_rough) .- 1e-9)
# `threaded` and `order` reach DGG.mapneighbors; the answer must not depend on
# either.
@assert parent(topographic_position_index(cell_raster, cell_grid; threaded=false)) ==
        parent(cell_tpi)
@assert parent(mapneighbors(_tpi_kernel, cell_raster, cell_grid, NeighborRings();
    needs=LOCAL_NEEDS, order=reverse(1:length(cell_raster)))) == parent(cell_tpi)

# ### (b) A value-only cell request builds no centroid table
#
# The payload for the default request is a zero-size singleton: there is no
# centroid storage to be found, at any grid size. When geometry *is* requested,
# the on-demand provider holds only the level grid and the radius.
ncells = length(cell_raster)
value_geom = neighborgeometry(cell_grid, NeighborRings(1), (Value(),))
ondemand_geom = neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
@assert value_geom.payload isa NoCentroids
@assert sizeof(value_geom.payload) == 0
@assert ondemand_geom.payload isa OnDemandCentroids
@assert _streamable(value_geom.fields) isa Requested   # -> DGG pass=Values()
@assert _streamable(ondemand_geom.fields) isa Absent   # -> DGG pass=Neighbors()
neighborgeometry(cell_grid, NeighborRings(1), (Value(),))            # warm up
neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
value_geom_bytes = @allocated neighborgeometry(cell_grid, NeighborRings(1), (Value(),))
ondemand_bytes = @allocated neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
# A materialized centroid table is 24 bytes per cell (three Float64 direction
# cosines). Both on-demand paths stay orders of magnitude below that.
@assert value_geom_bytes < 24 * ncells / 8
@assert ondemand_bytes < 24 * ncells / 8
precompute(ondemand_geom) # warm up
precomputed_bytes = @allocated precompute(ondemand_geom)
@assert precomputed_bytes >= 24 * ncells # The opt-in really does materialize

# ### (c) `precompute` is a performance choice, not a semantic one
cell_geom_steep = neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
cell_geom_fit = neighborgeometry(cell_grid, NeighborRings(1), PLANEFIT_NEEDS)
@assert precompute(cell_geom_steep).payload isa StoredCentroids
@assert precompute(cell_geom_steep).fields === cell_geom_steep.fields
@assert all(isequal.(
    parent(mapneighbors(_steepest_slope, cell_raster, cell_geom_steep)),
    parent(mapneighbors(_steepest_slope, cell_raster, precompute(cell_geom_steep)))))
@assert all(isequal.(
    parent(mapneighbors(_planefit_slope, cell_raster, cell_geom_fit)),
    parent(mapneighbors(_planefit_slope, cell_raster, precompute(cell_geom_fit)))))
# `precompute` on a rectilinear grid is a no-op: its tables are already O(rows).
@assert precompute(geom_geo).payload.geometry === geom_geo.payload.geometry

# A NaN cell propagates on the cell backend too.
nan_cell_values = copy(cell_values)
nan_cell_values[5] = NaN
nan_cell_raster = Raster(nan_cell_values, (DGG.Cells(cell_lookup),))
@assert isnan(parent(steepest_slope(nan_cell_raster))[5])
@assert isnan(parent(flow_direction(nan_cell_raster))[5])
@assert isnan(parent(slope(nan_cell_raster))[5])
# A cell that is not a neighbor of the NaN cell is untouched.
untouched = findfirst(p -> p != 5 && !(5 in DGG.neighbors(cell_lookup, p, 1)),
    1:length(cell_values))
@assert parent(steepest_slope(nan_cell_raster))[untouched] ==
        parent(cell_slope)[untouched]

# `eachneighbor` honors the request identically on both backends.
rect_records = collect(eachneighbor(
    neighborgeometry(grid_xy, NeighborRings(1), (Index(), Value())),
    raster_xy, CartesianIndex(1, 1)))
@assert length(rect_records) == 3
@assert all(r -> keys(r) == (:index, :value), rect_records)
cell_records = collect(eachneighbor(value_geom, cell_raster, cell_index))
@assert all(r -> keys(r) == (:value,), cell_records)
@assert length(cell_records) == length(DGG.neighbors(cell_lookup, 1, 1))
cell_records_full = collect(eachneighbor(
    neighborgeometry(cell_grid, NeighborRings(1), PLANEFIT_NEEDS), cell_raster, cell_index))
@assert all(r -> keys(r) == (:value, :distance, :bearing), cell_records_full)
@assert all(r -> r.distance > 0 && 0 <= r.bearing < 360, cell_records_full)

# ### (f) The sweep family
#
# Geomorphometry is the reference. Where the answer is uniquely determined — the
# settled surface always, the flood tree whenever no two cells tie — the
# assertion is `==`; where it is not, the assertion is the invariant, never a
# loosened tolerance.

# A logger that collects rather than prints, so the deliberate no-outlet warning
# is asserted instead of leaked into the output.
struct CollectedLogs <: Base.CoreLogging.AbstractLogger
    messages::Vector{String}
end
Base.CoreLogging.min_enabled_level(::CollectedLogs) = Base.CoreLogging.Debug
Base.CoreLogging.shouldlog(::CollectedLogs, args...) = true
Base.CoreLogging.catch_exceptions(::CollectedLogs) = false
Base.CoreLogging.handle_message(logger::CollectedLogs, level, message, _module, group,
    id, file, line; kwargs...) = push!(logger.messages, string(message))
collectlogs(f) =
    (l = CollectedLogs(String[]); (Base.CoreLogging.with_logger(f, l), l))

# One flood per fixture, so `down`, `acc` and the directions all describe the
# same tree — and so the tests exercise the same `table`-threading a composite
# operation uses.
function sweepparts(ras, grid; closed=nothing)
    table = neighbortable(grid)
    sweep = floodsweep(ras, grid; closed, table)
    acc = Vector{Float64}(undef, length(_data(ras)))
    _fillcellareas!(acc, grid)
    _accumulatedown!(acc, sweep.order, sweep.down)
    dirs = flowdirection(ras, grid; sweep, table)
    return sweep, reshape(acc, size(_data(ras))), dirs
end

# Walk `down` to its root, refusing to loop forever.
function downroot(down, p)
    steps = 0
    while down[p] != 0
        p = Int(down[p])
        (steps += 1) <= length(down) || error("`down` contains a cycle")
    end
    return p
end

# The linear `getindex` the traversal work arrays use must agree with the
# Cartesian one the record API uses.
@assert all(I -> ca[LinearIndices(size(ca))[I]] == ca[I[1], I[2]],
    CartesianIndices(size(ca)))
@assert all(I -> ca_geo[LinearIndices(size(ca_geo))[I]] == ca_geo[I[1], I[2]],
    CartesianIndices(size(ca_geo)))
@assert Base.IndexStyle(typeof(ca)) isa Base.IndexLinear
@assert Base.IndexStyle(typeof(ca_geo)) isa Base.IndexLinear
# A cell grid's areas are lazy: no residency, whatever the grid size.
@assert cellarea(cell_grid) isa CellAreas
@assert sizeof(cellarea(cell_grid)) < 8 * length(cell_raster)

# #### Fixture 1 — a bowl with a flat and a pit
#
# Two rectilinear DEMs, each run in `(X, Y)` and `(Y, X)` storage with unequal,
# sign-flipped spacing. `pit_z` has a uniform rim, so the whole interior lies
# below it and the fill raises everything; `bowl_z` has a graded rim, so only the
# enclosed cells rise — the two-cell flat and the single-cell pit spill over the
# ring of 9s.
pit_z = Float64[
    10 10 10 10 10 10
    10 8 7 7 6 10
    10 7 3 3 5 10
    10 6 3 1 4 10
    10 5 4 4 4 10
    10 10 10 10 10 10]
bowl_z = Float64[
    1 2 3 4 5 6
    2 9 9 9 9 7
    3 9 4 4 9 8
    4 9 4 1 9 9
    5 9 9 9 9 10
    6 7 8 9 10 11]
sweep_xs, sweep_ys = 0.0:2.0:10.0, 25.0:-5.0:0.0
asxy(z) = Raster(z, (X(sweep_xs), Y(sweep_ys)))
asyx(z) = Raster(permutedims(z), (Y(sweep_ys), X(sweep_xs)))

for z in (pit_z, bowl_z)
    # (2) Nested depressions: `settle` is the minimax fill, so it equals
    # `filldepressions` exactly. On Float64 there is no tolerance to hide in.
    @assert parent(settle(asxy(z))) == GM.filldepressions(z)
    # (8) Both storage orders. Settled elevation is a geographic quantity, so the
    # two layouts are transposes of each other.
    @assert parent(settle(asxy(z))) == permutedims(parent(settle(asyx(z))))
    @assert all(parent(settle(asxy(z))) .>= z)
    # (7) Integer elevations take the same path, promoted to float on the way in.
    @assert parent(settle(asxy(Int.(z)))) == GM.filldepressions(Int.(z))
    @assert parent(settle(asxy(UInt8.(z)))) == GM.filldepressions(UInt8.(z))
end
@assert parent(settle(asxy(pit_z)))[4, 4] == 10.0 # A uniform rim fills the whole bowl
@assert parent(settle(asxy(bowl_z)))[4, 4] == 9.0 # The pit spills over the 9-ring

function checkrectsweep(z)
    _, g = spatialparts(asxy(z))
    sweep, area, direction = sweepparts(asxy(z), g)
    @assert length(sweep.order) == length(z) # Every cell reachable, every slot used
    # Conservation: all area reaches a cell with no downstream neighbor.
    @assert sum(vec(area)[sweep.down .== 0]) ≈ sum(cellarea(g))
    @assert all(area .>= 10.0)
    # (3) No cycles: every `down` chain terminates at an outlet.
    @assert all(p -> sweep.down[downroot(sweep.down, p)] == 0, eachindex(vec(z)))
    # Pits are exactly the seeds; no interior pit survives the fill.
    @assert count(ispit, direction) == sweep.nseeds ==
            count(p -> isboundary(neighbortable(g), p), eachindex(vec(z)))
    return nothing
end
checkrectsweep(pit_z)
checkrectsweep(bowl_z)

# #### Fixture 2 — a tie-free surface, exact parity with Geomorphometry
#
# Ties make the flood tree non-unique, and the two priority queues break them
# differently (Geomorphometry's is `Dict`-backed; this one is array-backed). With
# every value distinct the tree is unique, so the direction rasters must agree
# cell for cell. `splitmix01` stands in for a seeded RNG, which is not a
# dependency here.
function splitmix01(k::Integer)
    h = UInt64(k) * 0x9E3779B97F4A7C15
    h ⊻= h >> 29
    h *= 0xBF58476D1CE4E5B9
    h ⊻= h >> 32
    h *= 0x94D049BB133111EB
    h ⊻= h >> 31
    return Float64(h >> 11) * (1 / 9007199254740992.0)
end
tf_n = 200
tf_z = [100.0 * splitmix01(i + tf_n * (j - 1)) + 1e-9 * (i + tf_n * (j - 1))
        for i in 1:tf_n, j in 1:tf_n]
@assert length(unique(tf_z)) == length(tf_z) # The premise of every `==` below
tf_ras = Raster(tf_z, (X(range(0.0; step=2.0, length=tf_n)),
    Y(range(1000.0; step=-5.0, length=tf_n))))
_, tf_grid = spatialparts(tf_ras)
tf_sweep, tf_acc, tf_dirs = sweepparts(tf_ras, tf_grid)
tf_settled = settle(tf_ras, tf_grid; sweep=tf_sweep)
gm_tf_acc, gm_tf_dirs = GM.flowaccumulation(tf_z; method=D8(), cellsize=(2.0, -5.0))

@assert tf_settled == GM.filldepressions(tf_z)
# The flood tree itself, compared through the LDD codec. No `_orient`, no
# `cellsize` sign logic — `storageoffset` already absorbed storage order.
@assert Int.(tf_dirs) == Int.(gm_tf_dirs)
# Geomorphometry accumulates in Float32 over the identical addition sequence, so
# the difference is its rounding, not a different answer.
@assert maximum(abs.(tf_acc .- Float64.(gm_tf_acc))) / maximum(tf_acc) < 1e-5
@assert all(isapprox.(Float32.(tf_acc), gm_tf_acc; rtol=1e-5))
@assert sum(vec(tf_acc)[tf_sweep.down .== 0]) ≈ sum(cellarea(tf_grid))
# Only the algorithm's own output and work arrays: `order` and `down` at 4 bytes
# each, `closed` at a bit, `acc` at 8, the directions at 1, and the queue. The
# rectilinear neighbor table is arithmetic and holds no memory at all.
@assert sizeof(neighbortable(tf_grid)) < 1024
flowaccumulation(tf_ras, tf_grid) # Warm up, then measure the steady state
@assert (@allocated flowaccumulation(tf_ras, tf_grid)) < 128 * length(tf_z)

# (8) Storage order, on the fixture where accumulation is determined: identical
# LDD codes, not permuted ones, because the codes are geographic.
tf_yx = Raster(permutedims(tf_z), (Y(range(1000.0; step=-5.0, length=tf_n)),
    X(range(0.0; step=2.0, length=tf_n))))
_, tf_grid_yx = spatialparts(tf_yx)
tf_sweep_yx, tf_acc_yx, tf_dirs_yx = sweepparts(tf_yx, tf_grid_yx)
@assert tf_acc == permutedims(tf_acc_yx)
@assert Int.(tf_dirs) == permutedims(Int.(tf_dirs_yx))
@assert tf_settled == permutedims(settle(tf_yx, tf_grid_yx; sweep=tf_sweep_yx))

# (1) A single-cell pit: filled to its spill elevation, and routed rather than
# left a sink. Its direction points uphill on the raw DEM, which is what correct
# depression routing looks like.
tf_table = neighbortable(tf_grid)
tf_pits = findall(eachindex(tf_z)) do p
    !isboundary(tf_table, p) && all(((k, q),) -> tf_z[q] > tf_z[p], slots(tf_table, p))
end
@assert !isempty(tf_pits)
let p = first(tf_pits)
    @assert tf_settled[p] > tf_z[p]
    @assert !ispit(tf_dirs[p])
    @assert tf_settled[p] == max(tf_z[p], tf_settled[tf_sweep.down[p]])
end

# #### (3) A flat draining to two outlets
#
# A plateau at 5 inside a rim at 9, with two low notches in the rim. Which flat
# cell goes to which notch is heap-order dependent, so the assertions are the
# partition and the conservation, not the assignment.
flat_z = fill(5.0, 5, 7)
flat_z[1, :] .= 9.0
flat_z[5, :] .= 9.0
flat_z[:, 1] .= 9.0
flat_z[:, 7] .= 9.0
flat_z[1, 3] = 1.0
flat_z[5, 5] = 2.0
_, flat_grid = spatialparts(flat_z; spacing=(2.0, -5.0))
flat_sweep, flat_acc, flat_dirs = sweepparts(flat_z, flat_grid)
flat_lin = LinearIndices(flat_z)
flat_interior = vec(flat_lin[2:4, 2:6])
@assert all(p -> flat_sweep.down[downroot(flat_sweep.down, p)] == 0, eachindex(vec(flat_z)))
flat_roots = unique(downroot(flat_sweep.down, p) for p in flat_interior)
@assert Set(flat_roots) == Set((flat_lin[1, 3], flat_lin[5, 5])) # The two notches, both used
@assert flat_acc[1, 3] + flat_acc[5, 5] ≈ (length(flat_interior) + 2) * 10.0
@assert sum(vec(flat_acc)[flat_sweep.down .== 0]) ≈ sum(cellarea(flat_grid))
@assert all(parent(settle(flat_z, flat_grid; sweep=flat_sweep))[flat_interior] .== 5.0)

# #### (4, 6) Nodata
#
# By default NaN participates, exactly as Geomorphometry does: NaN sorts last in
# the queue, so those cells are visited at the end, get a parent, and inject their
# own area into the accumulation. `closed = isnan.(z)` is the one-liner that
# excludes them instead.
nan_z = Float64[20.0 - i - 0.3j for i in 1:8, j in 1:8]
nan_z[3:5, 3:5] .= NaN # A 3x3 hole, so its center's whole neighborhood is nodata
_, nan_grid = spatialparts(nan_z; spacing=(2.0, -5.0))
nan_table = neighbortable(nan_grid)
nan_holes = findall(isnan, vec(nan_z))
@assert all(((k, q),) -> q == 0 || isnan(nan_z[q]), slots(nan_table, LinearIndices(nan_z)[4, 4]))

nan_sweep, nan_acc, nan_dirs = sweepparts(nan_z, nan_grid)
nan_settled = settle(nan_z, nan_grid; sweep=nan_sweep)
@assert length(nan_sweep.order) == length(nan_z)
@assert all(p -> nan_sweep.down[p] != 0, nan_holes) # Visited, and given a parent
@assert all(p -> isnan(nan_settled[p]), nan_holes)  # But with no settled elevation
@assert all(p -> isfinite(nan_settled[p]), setdiff(eachindex(vec(nan_z)), nan_holes))
@assert sum(vec(nan_acc)[nan_sweep.down .== 0]) ≈ sum(cellarea(nan_grid))

closed_sweep, closed_acc, closed_dirs = sweepparts(nan_z, nan_grid; closed=isnan.(nan_z))
@assert length(closed_sweep.order) == length(nan_z) - length(nan_holes)
@assert all(p -> closed_sweep.down[p] == 0, nan_holes)
@assert all(p -> closed_acc[p] == 10.0, nan_holes) # Closed cells emit and receive nothing
@assert all(p -> ispit(closed_dirs[p]), nan_holes)
let reached = (vec(closed_sweep.down) .== 0) .& .!isnan.(vec(nan_z))
    @assert sum(vec(closed_acc)[reached]) ≈ sum(cellarea(nan_grid)) - length(nan_holes) * 10.0
end

# (6) A raster with nothing to seed says so instead of looping or silently
# returning `acc == cellarea`.
allnan_z = fill(NaN, 4, 4)
_, allnan_grid = spatialparts(allnan_z)
allnan_sweep, allnan_logs = collectlogs() do
    floodsweep(allnan_z, allnan_grid; closed=isnan.(allnan_z))
end
@assert any(m -> occursin("every cell is closed", m), allnan_logs.messages)
@assert isempty(allnan_sweep.order) && all(==(0), allnan_sweep.down)

# #### (11) `closed` disconnecting a component
#
# A ring of closed cells encloses the bowl's four inner cells, so fewer cells are
# reachable than `order` was sized for. The truncation is load-bearing: without
# it the untruncated tail of ones would inflate one arbitrary cell.
_, bowl_grid = spatialparts(bowl_z; spacing=(2.0, -5.0))
bowl_ring = falses(6, 6)
bowl_ring[2:5, 2] .= true
bowl_ring[2:5, 5] .= true
bowl_ring[2, 2:5] .= true
bowl_ring[5, 2:5] .= true
cut_sweep, cut_acc, cut_dirs = sweepparts(bowl_z, bowl_grid; closed=bowl_ring)
cut_inner = vec(LinearIndices((6, 6))[3:4, 3:4])
@assert length(cut_sweep.order) == length(bowl_z) - count(bowl_ring) - length(cut_inner)
@assert length(cut_sweep.order) < length(bowl_z) - count(bowl_ring)
@assert all(p -> cut_sweep.down[p] == 0, cut_inner)
@assert all(p -> cut_acc[p] == 10.0, cut_inner) # Never touched by the `ones` tail
let reached = (vec(cut_sweep.down) .== 0) .& .!vec(bowl_ring)
    reached[cut_inner] .= false
    @assert sum(vec(cut_acc)[reached]) ≈ length(cut_sweep.order) * 10.0
end

# #### (10) Degenerate grids
one_acc, one_dirs = flowaccumulation(fill(5.0, 1, 1))
@assert one_acc == fill(1.0, 1, 1) # Its own outlet, carrying its own area
@assert ispit(only(one_dirs))
@assert only(floodsweep(fill(5.0, 1, 1), last(spatialparts(fill(5.0, 1, 1)))).order) == 1
empty_acc, empty_dirs = flowaccumulation(zeros(0, 0))
@assert isempty(empty_acc) && isempty(empty_dirs)

# #### (12) Multi-direction methods refuse clearly
dinf_error = try flowaccumulation(asxy(bowl_z); method=DInf()); nothing catch e; e end
@assert dinf_error isa ArgumentError
@assert occursin("Bearing", sprint(showerror, dinf_error))
@assert occursin("D8", sprint(showerror, dinf_error))
fd8_error = try flowaccumulation(asxy(bowl_z); method=FD8()); nothing catch e; e end
@assert fd8_error isa ArgumentError

# #### (13) Round-trip through the façade
matrix_acc, matrix_dirs = flowaccumulation(bowl_z; spacing=(2.0, -5.0))
@assert matrix_acc isa Matrix{Float64}
@assert matrix_dirs isa Matrix{FlowDirection{LDD,UInt8}}
raster_acc, raster_dirs = flowaccumulation(asxy(bowl_z))
@assert raster_acc isa Raster && raster_dirs isa Raster
@assert (Rasters.name(raster_acc), Rasters.name(raster_dirs)) ==
        (:flowaccumulation, :flowdirection)
@assert isnan(Rasters.missingval(raster_acc))
@assert isnothing(Rasters.missingval(raster_dirs))
@assert Rasters.name(settle(asxy(bowl_z))) == :settle
@assert Rasters.name(flowdirection(asxy(bowl_z))) == :flowdirection
# Algorithm keywords reach the inner method through the façade unchanged.
@assert isnan(parent(settle(Raster(nan_z, (X(1.0:8.0), Y(8.0:-1.0:1.0)));
    closed=isnan.(nan_z)))[4, 4])
# Rule A: below the façade every output is a plain array. Rasters' missingval
# machinery is reached exactly once, in `rebuildoutput`.
inner_acc, inner_dirs = flowaccumulation(asxy(bowl_z), last(spatialparts(asxy(bowl_z))))
@assert inner_acc isa Matrix{Float64} && !(inner_acc isa Raster)
@assert inner_dirs isa Matrix{FlowDirection{LDD,UInt8}} && !(inner_dirs isa Raster)
# A geographic raster keeps its CRS, and its per-row cell areas still conserve.
geo_acc, geo_dirs = flowaccumulation(geo_raster)
@assert Rasters.crs(geo_acc) == Rasters.crs(geo_raster)
@assert sum(vec(parent(geo_acc))[floodsweep(geo_raster, grid_geo).down .== 0]) ≈
        sum(cellarea(grid_geo))

# #### Fixture 3 — cell grids
#
# `cell_raster` is a *complete* level-2 grid, and a sphere has no edges: the
# boundary seed set is empty, so the no-outlet fallback must fire rather than let
# an empty queue return `acc == cellarea` everywhere.
sphere_table = neighbortable(cell_grid)
@assert !any(p -> isboundary(sphere_table, p), 1:length(sphere_table))
@assert isempty(_boundaryseeds(sphere_table, falses(length(cell_raster))))
sphere_sweep, sphere_logs = collectlogs() do
    floodsweep(cell_raster, cell_grid; table=sphere_table)
end
@assert any(m -> occursin("no domain boundary", m), sphere_logs.messages)
@assert length(sphere_sweep.order) == length(cell_raster)
sphere_acc = let a = Vector{Float64}(undef, length(cell_raster))
    _fillcellareas!(a, cell_grid)
    _accumulatedown!(a, sphere_sweep.order, sphere_sweep.down)
end
@assert sum(sphere_acc[sphere_sweep.down .== 0]) ≈ sum(cellarea(cell_grid)) ≈
        4π * AUTHALIC_RADIUS_M^2
# A one-ring table is the contract; a wider request is refused, not clipped.
wide_table_error = try neighbortable(cell_grid, NeighborRings(2)); nothing catch e; e end
@assert wide_table_error isa ArgumentError
# The adjacency table is built once per public call and threaded through the
# flood and the direction encoding — a second build would double this.
collectlogs() do
    neighbortable(cell_grid)
    flowaccumulation(cell_raster, cell_grid) # Warm up before measuring
    @assert (@allocated flowaccumulation(cell_raster, cell_grid)) <
            2 * (@allocated neighbortable(cell_grid))
end

# A *region* of a cell grid has a rim, so Geomorphometry can run on it too. The
# values carry a per-position ramp to make the flood tree unique.
sub_region = DGG.subtree(DGG.IGeo7System(),
    DGG.cellindex(DGG.levelgrid(DGG.IGeo7System(), 0), 5), 3)
sub_lookup = DGG.CellLookup(sub_region)
sub_cells = sub_lookup.cells
sub_vals = [10_000.0 * DGG.cell_centroid(sub_region, c)[3] + 1e-6 * p
            for (p, c) in enumerate(sub_cells)]
@assert length(unique(sub_vals)) == length(sub_vals)
sub_raster = Raster(sub_vals, (DGG.Cells(sub_lookup),))
_, sub_grid = spatialparts(sub_raster; spatialdims=DGG.Cells)
sub_table = neighbortable(sub_grid)

# The seed set read off the zero slots is exactly DGG's own rim, and exactly the
# rim Geomorphometry rebuilds per cell from `neighborcount(complete, cv[p])`.
@assert _boundaryseeds(sub_table, falses(length(sub_cells))) ==
        sort(collect(DGG.border(sub_cells)))
# The invariant the IGeo7 codec rests on: a complete-width row preserves slot
# identity, and on IGeo7 slot `k` is direction code `k`. That is why no ring is
# rebuilt and no relative cell is constructed.
@assert all(1:length(sub_cells)) do p
    all(((k, q),) -> q == 0 || DGG.directioncode(sub_cells[q] - sub_cells[p]) == k,
        slots(sub_table, p))
end

# `nslots` is the complete degree, which is what makes a slot code portable
# between regions: a clipped row would renumber it.
@assert nslots(tf_table, 1) == 8 == length(collect(slots(tf_table, 1)))
@assert all(p -> nslots(sub_table, p) == length(collect(slots(sub_table, p))),
    1:length(sub_cells))
@assert Set(nslots(sub_table, p) for p in 1:length(sub_cells)) ⊆ Set((5, 6))

sub_sweep, sub_acc, sub_dirs = sweepparts(sub_raster, sub_grid)
gm_sub_acc, gm_sub_dirs = GM.flowaccumulation(sub_raster; method=D8())
@assert eltype(sub_dirs) == FlowDirection{LDD,UInt8}
@assert Int.(sub_dirs) == Int.(parent(gm_sub_dirs))
@assert maximum(abs.(sub_acc .- Float64.(parent(gm_sub_acc)))) / maximum(sub_acc) < 1e-5
@assert sum(sub_acc[sub_sweep.down .== 0]) ≈ sum(cellarea(sub_grid))
@assert count(ispit, sub_dirs) == sub_sweep.nseeds == length(collect(DGG.border(sub_cells)))

# Geomorphometry has no `filldepressions` for cell grids, so `settle` is checked
# against its defining invariants instead.
sub_settled = settle(sub_raster, sub_grid; sweep=sub_sweep)
@assert all(sub_settled .>= sub_vals)
@assert all(eachindex(sub_vals)) do p
    d = sub_sweep.down[p]
    d == 0 ? sub_settled[p] == sub_vals[p] :
    sub_settled[p] == max(sub_vals[p], sub_settled[d])
end

# The generic seam. A cell system with no relative-cell arithmetic gets the ring
# slot itself, which decodes back to a position in O(1) off the same row — no
# ring rebuild and no `findfirst`. It is exercised here through the IGeo7 table
# because slot `k` is what the IGeo7 codec is built on in the first place.
slot_dirs = Vector{FlowDirection{RingSlot,UInt8}}(undef, length(sub_cells))
_cellencode!(slot_dirs, sub_sweep.down, sub_table, FlowDirection{RingSlot,UInt8})
@assert all(p -> downstreamposition(sub_table, p, slot_dirs[p]) == sub_sweep.down[p],
    eachindex(sub_cells))
@assert count(ispit, slot_dirs) == sub_sweep.nseeds
@assert Int.(sub_dirs) == [IGEO7_TO_LDD[Int(d) + 1] for d in slot_dirs]
@assert sprint(show, first(filter(ispit, slot_dirs))) == "·"
@assert sprint(show, first(filter(!ispit, slot_dirs))) != "·"

# (5) The tutorial's nodata overhang, in miniature: nodata cells sitting on the
# region rim, excluded with `closed`. This is the configuration in which
# Geomorphometry seeds rim cells it was told were closed and then writes past the
# end of `order`.
sub_nan = copy(sub_vals)
sub_nan[first(collect(DGG.border(sub_cells)), 10)] .= NaN
sub_nan_raster = Raster(sub_nan, (DGG.Cells(sub_lookup),))
sub_nan_sweep, sub_nan_acc, _ = sweepparts(sub_nan_raster, sub_grid; closed=isnan.(sub_nan))
@assert length(sub_nan_sweep.order) == length(sub_nan) - count(isnan, sub_nan)
@assert all(p -> sub_nan_sweep.down[p] == 0, findall(isnan, sub_nan))

println("PoC v3 passed: requested neighbor fields, on-demand cell geometry, " *
        "NaN semantics, named local statistics, and the sweep family — one " *
        "priority flood behind `neighbortable`, with settle, D8 accumulation " *
        "and LDD directions matching Geomorphometry on both backends.")
