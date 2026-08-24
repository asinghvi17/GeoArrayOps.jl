"""
    GeomorphometryNGDiscreteGlobalGridsExt

The cell-grid backend: every method of GeomorphometryNG that needs a discrete
global grid to answer.

The split is by *what a method calls*, not by what it dispatches on. All of the
cell types — `CellGrid`, the centroid payloads, `CellAreas`, `CellNeighborTable`
— are plain storage and live in the package, as do the generic traversal and
direction-encoding passes written against them. What is here is the set of
methods that reach into DGG for topology (`neighbors`, `adjacency`,
`mapneighbors`), geometry (`cell_centroid`), area (`cell_area`) or cell identity
(`cellid`, `localindex`, `SubsetIndexedCell`).
"""
module GeomorphometryNGDiscreteGlobalGridsExt

using Rasters
import DiscreteGlobalGrids as DGG
using Geomorphometry: FlowDirection, LDD

import GeomorphometryNG as NG
# Extended here. Every one of these already has a rectilinear or grid-agnostic
# method in the package, so these are additions, not definitions.
import GeomorphometryNG: cellarea, cellindices, cellkey, eachneighbor, edgegeometry,
    fromcentroid, isspatialdim, mapneighbors, neighbortable, storagekey,
    storageposition, _cellareaof, _celldirectiontype, _iscelllookup, _materialize
# Used, not extended.
using GeomorphometryNG: AUTHALIC_RADIUS_M, Absent, CellGrid, CellNeighborTable,
    Index, NeighborGeometry, NeighborRings, NoCentroids, OnDemandCentroids,
    Requested, StoredCentroids, Value, manifold, neighborrecord, _asked,
    _cellmanifold, _edgegeometry, _neighborvalue, _streamable

# ## Grid construction
#
# A `Cells` dimension is a spatial dimension, and a `CellLookup` is what makes
# `spatialparts` take the cell branch. Neither question is answerable without
# DGG loaded, and `false` is the right answer when it is not.

isspatialdim(::Type{<:DGG.Cells}) = true
_iscelllookup(::DGG.CellLookup) = true

function NG.CellGrid(celldim; manifold=nothing)
    lookup = Rasters.lookup(celldim)
    lookup isa DGG.CellLookup ||
        throw(ArgumentError("Dimension must contain a DGG.CellLookup; found $(typeof(lookup))"))
    cells = lookup.cells
    levelgrid = DGG.levelgrid(DGG.system(cells), DGG.level(cells))
    return CellGrid{typeof(cells),typeof(levelgrid)}(
        _cellmanifold(manifold), cells, levelgrid)
end

# ## Cell identity
#
# Stable cell handles, and their conversion to storage positions.

function cellindices(ras, grid::CellGrid)
    return (DGG.SubsetIndexedCell(grid.cells[p], p) for p in eachindex(grid.cells))
end

storagekey(::CellGrid, cell) = DGG.localindex(cell)

# The random-access half of the bijection, which a queue-driven traversal needs
# and `cellindices` — a generator here — cannot serve.
cellkey(ras, grid::CellGrid, p::Int) =
    DGG.SubsetIndexedCell((@inbounds grid.cells[p]), p)
storageposition(::CellGrid, cell) = DGG.localindex(cell)

# ## Geometry
#
# `precompute`'s materialization, and the two centroid fetches the on-demand
# provider makes — one per visited cell, one per edge.

_materialize(p::OnDemandCentroids, grid::CellGrid) = StoredCentroids(p.levelgrid, p.radius,
    [DGG.cell_centroid(p.levelgrid, c) for c in grid.cells])

@inline fromcentroid(p::OnDemandCentroids, cell, position) =
    DGG.cell_centroid(p.levelgrid, cell)

@inline edgegeometry(p::OnDemandCentroids, from, cell, position) =
    _edgegeometry(from, DGG.cell_centroid(p.levelgrid, cell), p.radius)

# ## Neighbor iteration
#
# Queue-driven algorithms hoist the geometry once and reuse it while visiting
# cells; the records carry exactly the requested fields.

function eachneighbor(geom::NeighborGeometry{<:CellGrid}, ras, cell)
    cells = geom.grid.cells
    payload = geom.payload
    request = geom.fields
    wantvalue = _asked(Value, request)
    p0 = DGG.localindex(cell)
    from = fromcentroid(payload, cells[p0], p0)
    positions = DGG.neighbors(cells, p0, geom.rings.k)
    return (
        begin
            c = cells[p]
            nb = DGG.SubsetIndexedCell(c, p)
            neighborrecord(request, nb, _neighborvalue(wantvalue, ras, nb),
                edgegeometry(payload, from, c, p))
        end for p in positions
    )
end

# ## The cell driver

# DGG performs cell traversal and infers the output type. A request that is a
# subset of `{Value}` rides DGG's streaming `Values()` pass — the same code path
# Geomorphometry's own DGG extension uses — so value-only kernels reach baseline
# speed by construction. Anything wider states its fields in one `needs=`
# request and reads them off the membership clip the sweep already made.
function mapneighbors(f, ras::AbstractVector, geom::NeighborGeometry{<:CellGrid};
        threaded=true, order=DGG.StorageOrder())
    # DGG's traversal always hands the callback the one-ring, so a wider request
    # ignores it and re-queries the topology per cell.
    geom.rings.k == 1 || return _cellmapwiderings(f, ras, geom, threaded, order)
    return _cellmapneighbors(_streamable(geom.fields), f, ras, geom, threaded, order)
end

function _cellmapneighbors(::Requested, f::F, ras, geom, threaded, order) where {F}
    request = geom.fields
    return DGG.mapneighbors(ras; pass=DGG.Values(), threaded, order) do cell, value, values
        f(cell, value,
            (neighborrecord(request, nothing, v, nothing) for v in values))
    end
end

# ### The field-request route
#
# A request that reads an index or a geometry states those fields to DGG rather
# than fetching them per neighbor from inside the callback: one `needs=` tuple,
# one membership clip, and the callback reaches back into nothing. The rings
# arrive *field-major* — `rings[j]` is need `j` for every neighbor, slot `i` of
# every ring naming the same neighbor in `DGG.neighbors` order — so the
# neighbor-major records the kernels read are a `zip` on this side.
#
# `DGG.Index(DGG.Local())` is in every request, because the callback's first
# argument is the visited cell's own handle and the request tuple is the only
# channel that names the visited cell's index. Its *ring* is spare work when the
# record carries no `index`; `DGG.Cell()` joins it only when the record does.
function _cellmapneighbors(::Absent, f::F, ras, geom, threaded, order) where {F}
    return _cellmapfields(geom.payload, _asked(Index, geom.fields), f, ras, geom,
        threaded, order)
end

# Distance and bearing are the kernel's own arithmetic over two streamed
# centroids: DGG answers on the unit sphere, and the radius stays here.
function _cellmapfields(payload::OnDemandCentroids, ::Absent, f::F, ras, geom,
        threaded, order) where {F}
    request = geom.fields
    cells = geom.grid.cells
    R = payload.radius
    needs = (DGG.Index(DGG.Local()), DGG.Value(ras), DGG.Centroid())
    return DGG.mapneighbors(ras; needs, threaded, order) do center, rings
        position, value, from = center
        _, values, centroids = rings
        f(DGG.SubsetIndexedCell((@inbounds cells[position]), position), value,
            (neighborrecord(request, nothing, v, _edgegeometry(from, to, R))
             for (v, to) in zip(values, centroids)))
    end
end

# The record names its neighbor, so the cell ids ride the same clip and the
# handles are rebuilt from two rings instead of one array read per neighbor.
function _cellmapfields(payload::OnDemandCentroids, ::Requested, f::F, ras, geom,
        threaded, order) where {F}
    request = geom.fields
    R = payload.radius
    needs = (DGG.Cell(), DGG.Index(DGG.Local()), DGG.Value(ras), DGG.Centroid())
    return DGG.mapneighbors(ras; needs, threaded, order) do center, rings
        cell, position, value, from = center
        cellring, positions, values, centroids = rings
        f(DGG.SubsetIndexedCell(cell, position), value,
            (neighborrecord(request, DGG.SubsetIndexedCell(c, p), v,
                 _edgegeometry(from, to, R))
             for (c, p, v, to) in zip(cellring, positions, values, centroids)))
    end
end

# No geometry was requested, so no centroid is: this is the `Index`-only shape,
# the one `_streamable` turns away for its index alone.
function _cellmapfields(::NoCentroids, ::Requested, f::F, ras, geom,
        threaded, order) where {F}
    request = geom.fields
    needs = (DGG.Cell(), DGG.Index(DGG.Local()), DGG.Value(ras))
    return DGG.mapneighbors(ras; needs, threaded, order) do center, rings
        cell, position, value = center
        cellring, positions, values = rings
        f(DGG.SubsetIndexedCell(cell, position), value,
            (neighborrecord(request, DGG.SubsetIndexedCell(c, p), v, nothing)
             for (c, p, v) in zip(cellring, positions, values)))
    end
end

# `precompute` already paid O(ncells) to store the centroids, so a `Centroid()`
# request would recompute exactly what it bought. The stored payload keeps the
# handle route and reads its table.
function _cellmapfields(::StoredCentroids, ::Any, f::F, ras, geom,
        threaded, order) where {F}
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
    from = fromcentroid(payload, DGG.cellid(cell), DGG.localindex(cell))
    return (
        neighborrecord(request, nb, _neighborvalue(wantvalue, ras, nb),
            edgegeometry(payload, from, DGG.cellid(nb), DGG.localindex(nb)))
        for nb in nbrs
    )
end

# ## Cell areas

_cellareaof(levelgrid, c, R) = DGG.cell_area(levelgrid, c) * R^2
# Z7 cells provide a closed-form area in square meters on the authalic sphere;
# rescale it when the grid's radius differs.
_cellareaof(levelgrid, c::DGG.Z7Cell, R) =
    DGG.IGeo7.cell_area(DGG.rawid(c)) * (R / AUTHALIC_RADIUS_M)^2

cellarea(grid::CellGrid, cell) =
    _cellareaof(grid.levelgrid, DGG.cellid(cell), manifold(grid).radius)

# ## The traversal primitive

function neighbortable(grid::CellGrid, rings::NeighborRings=NeighborRings(1))
    rings.k == 1 ||
        throw(ArgumentError("the sweep family needs a one-ring neighbor table"))
    return CellNeighborTable(DGG.adjacency(grid.cells; halo=:mark, threaded=true))
end

# ## Direction encoding
#
# Z7 cells have relative-cell arithmetic, so they get Geomorphometry's LDD
# numpad rather than the generic ring slot.

_celldirectiontype(::Type{<:DGG.Z7Cell}) = FlowDirection{LDD,UInt8}

end # module GeomorphometryNGDiscreteGlobalGridsExt
