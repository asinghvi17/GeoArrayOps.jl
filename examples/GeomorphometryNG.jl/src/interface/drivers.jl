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

# `eachneighbor(::NeighborGeometry{<:CellGrid}, ...)` is in the
# DiscreteGlobalGrids extension.

# Only a request with no `Index`, `Distance` or `Bearing` can stream values.
# This is a question about the request tuple alone, so it is answered here even
# though only the cell driver in the DiscreteGlobalGrids extension asks it.
@inline _streamable(request::Tuple) = _streamable(_asked(Index, request),
    _asked(Distance, request), _asked(Bearing, request))
@inline _streamable(::Absent, ::Absent, ::Absent) = Requested()
@inline _streamable(::Any, ::Any, ::Any) = Absent()

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

# ## `mapwindow`: named north-up windows for rectilinear grids
#
# Window kernels access neighbors by geographic direction, independent of
# storage order. For example, `w.northeast` always refers to the cell northeast
# of the center. This interface is not defined for hexagonal cell grids.

function mapwindow(f, ras, grid::RectilinearGrid; radius::Int=1)
    stencil = northupwindow(grid, radius)
    interior = CartesianIndices(
        map(r -> (first(r) + radius):(last(r) - radius), axes(ras)))
    # The kernel is handed the *metric* spacing of the window it actually reads,
    # so a dilated window reports `radius * spacing`. Geomorphometry's
    # `scaled8nb(radius)` dilates the offsets but keeps dividing by the undilated
    # `cellsize`, which makes its `radius > 1` derivatives off by `radius` (and
    # by `radius^2` for the second derivatives).
    return Stencils.mapstencil(stencil, ras, cellindices(ras, grid);
        boundary=Stencils.Remove(zero(eltype(ras)))) do hood, I
        I in interior || return NaN
        f(hood, metricspacing(grid, I) .* radius)
    end
end

mapwindow(f, ras, ::CellGrid; kw...) = throw(ArgumentError(
    "windowed drivers require a rectilinear grid; use a neighbor-fit method (e.g. PlaneFit())"))

metricspacing(grid::RectilinearGrid{<:Planar}, I) = abs.(grid.spacing)

function metricspacing(grid::RectilinearGrid{<:Spherical}, I)
    R = grid.manifold.radius
    sx, sy = abs.(grid.spacing)
    lat = grid.lookups[2][I[yaxisnum(grid)]]
    return (deg2rad(sx) * cosd(lat) * R, deg2rad(sy) * R)
end
