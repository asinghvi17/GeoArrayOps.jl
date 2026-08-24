# ## `neighbortable`: the traversal primitive
#
# A priority flood is random access in a data-dependent order, which is the one
# access pattern the record API is not built for: `neighborgeometry` serves
# whole-grid sweeps in storage order. So the sweep family gets its own hoisted
# structure — the traversal analogue of `neighborgeometry` — exposing the one
# thing the flood needs and the records do not carry: a *slot-stable,
# position-keyed* adjacency. That makes it grid-interface machinery rather than
# an algorithm, which is why it lives here beside `neighborgeometry`; the
# algorithms that consume it are `algorithms/hydrology.jl` and
# `algorithms/flowdir.jl`.
#
# Everything below is keyed by storage position (`1:ncells`) on both backends,
# never by cell id.
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

# `neighbortable(::CellGrid, ::NeighborRings)`, which asks DGG for the CSR
# adjacency, is in the DiscreteGlobalGrids extension.

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

# ## Per-slot edge geometry: the traversal table's missing half
#
# `neighbortable` is position-keyed and carries no geometry; `neighborgeometry`
# carries geometry but streams records in storage order, keyed by the array's
# own index type. A queue-driven algorithm that needs *both* — cost distance,
# and the multi-direction flow methods, which split flow by bearing — falls
# between them.
#
# `slotgeometry` closes that gap without a third representation: it reuses the
# geometry providers `neighborgeometry` already builds, and exposes them keyed
# by `(p, k, q)`, which is the traversal table's key. It caches nothing that
# `neighborgeometry` would not: a rectilinear grid gets the same one-table (or
# one-table-per-row) payload, and a cell grid the same on-demand centroids.
#
# The two-step `cellorigin` / `edgeat` shape is what keeps the cell backend from
# refetching the from-centroid per edge: hoist once per cell, read per slot.

struct SlotGeometry{G<:AbstractGridSpec,P,C}
    grid::G
    payload::P
    cart::C
end

slotgeometry(grid::RectilinearGrid) = SlotGeometry(grid,
    _rectgeometry(Requested(), grid, logicaloffsets(NeighborRings(1))),
    CartesianIndices(gridsize(grid)))

slotgeometry(grid::CellGrid) =
    SlotGeometry(grid, _cellgeometry(Requested(), grid), nothing)

# On a rectilinear grid the "origin" is the whole slot table for the cell's row,
# so an edge read is a tuple index.
@inline cellorigin(sg::SlotGeometry{<:RectilinearGrid}, p::Int) =
    geometryat(sg.payload, @inbounds sg.cart[p])
@inline edgeat(::SlotGeometry{<:RectilinearGrid}, table, ::Int, k::Int, ::Int) =
    @inbounds table[k]

@inline cellorigin(sg::SlotGeometry{<:CellGrid}, p::Int) =
    fromcentroid(sg.payload, (@inbounds sg.grid.cells[p]), p)
@inline edgeat(sg::SlotGeometry{<:CellGrid}, from, ::Int, ::Int, q::Int) =
    edgegeometry(sg.payload, from, (@inbounds sg.grid.cells[q]), q)
