# ## Direction output
#
# The direction half of the sweep family: the conventions, the two codecs and
# the public `flowdirection`. The flood that produces `down` is in
# `hydrology.jl`.
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
_celldirectiontype(::Type) = FlowDirection{RingSlot,UInt8}
# The DiscreteGlobalGrids extension adds the `Z7Cell` method that selects `LDD`.

"""
    flowdirection(ras, grid; closed=nothing, table=…, sweep=…)

The downstream direction of each cell on the settled surface, as
`FlowDirection{LDD}` on rectilinear and IGeo7 grids and `FlowDirection{RingSlot}`
on other cell grids. Outlets, and cells the flood never visited, are pits.

A cell whose value is nodata is *not* a pit: the flood visits it last and gives
it a parent, matching Geomorphometry. Pass `closed` to exclude it instead.
"""
function flowdirection(ras, grid::AbstractGridSpec; method=D8(), closed=nothing,
        table=neighbortable(grid, NeighborRings(1)),
        sweep=floodsweep(ras, grid; closed, table))
    return _flowdirection(method, ras, grid, table, sweep)
end

function _flowdirection(::D8, ras, grid::AbstractGridSpec, table, sweep)
    out = Array{_directiontype(grid)}(undef, size(_data(ras))) # Rule A: plain Array
    _encodedirections!(vec(out), sweep.down, grid, table)
    return out
end

_flowdirection(method::Union{DInf,FD8}, ras, grid::AbstractGridSpec, table, sweep) =
    _multidirections(ras, grid, table,
        flowpartition(method, ras, grid; table, sweep))

_flowdirection(method, ras, grid, table, sweep) = throw(ArgumentError(
    "flowdirection has no rule for $(nameof(typeof(method))). The supported " *
    "methods are D8(), DInf() and FD8()."))

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

# ## Multi-direction output
#
# A cell that splits its flow has several downstream neighbors, so its direction
# is a *set*. Geomorphometry spells that as `D8D`, a bitmask over the eight
# numpad directions, and a rectilinear grid gets exactly that — bit-comparable
# with Geomorphometry. A hexagonal ring has no numpad, so it gets the bitmask
# over ring slots, which is the set-valued analogue of [`RingSlot`](@ref).

"""
Bitmask over ring slots: bit `k - 1` set means the cell drains to ring slot `k`.
Zero is a pit. This is [`RingSlot`](@ref) generalized to split flow, and the
cell-grid counterpart of Geomorphometry's `D8D`.
"""
struct RingMask <: FlowDirectionConvention end

GM.ismulti(::Type{RingMask}) = true
GM.ispit(d::FlowDirection{RingMask}) = iszero(d.value)
GM._arrow(::Type{RingMask}, d::Integer) =
    iszero(d) ? '·' :
    isone(count_ones(d)) ? Char(0x2460 + trailing_zeros(d)) : '✳'

_multidirectiontype(::RectilinearGrid) = FlowDirection{D8D,UInt8}
_multidirectiontype(::CellGrid) = FlowDirection{RingMask,UInt8}

function _multidirections(ras, grid::AbstractGridSpec, table, part)
    out = Array{_multidirectiontype(grid)}(undef, size(_data(ras)))
    _encodemulti!(vec(out), part, grid, table)
    return out
end

# One ring is what the numpad encodes, and `logicaloffsets` is that ring, so the
# D8D bit is read straight off the slot's logical offset — the same argument
# `_lddcode` makes, in the other convention.
_d8dcode((dx, dy)) = -1 <= dx <= 1 && -1 <= dy <= 1 ?
    GM._d8_ci2dir[CartesianIndex(dx, dy)] :
    throw(ArgumentError("the D8D bitmask encodes one ring; got offset ($dx, $dy)"))

_encodemulti!(out, part, ::RectilinearGrid, t::RectNeighborTable{K}) where {K} =
    _rectencodemulti!(out, part, t, ntuple(k -> _d8dcode(t.logical[k]), Val(K)))

function _rectencodemulti!(out, part, t, codes::NTuple{K,UInt8}) where {K}
    Threads.@threads for r in _chunkranges(length(out))
        @inbounds for p in r
            code = 0x00
            for e in partitionrange(part, p)
                q = Int(part.targets[e])
                for (k, qq) in slots(t, p)
                    qq == q && (code |= codes[k]; break)
                end
            end
            out[p] = FlowDirection{D8D,UInt8}(code)
        end
    end
    return out
end

function _encodemulti!(out, part, ::CellGrid, t::CellNeighborTable)
    Threads.@threads for r in _chunkranges(length(out))
        @inbounds for p in r
            code = 0x00
            for e in partitionrange(part, p)
                k = _downslot(t, p, Int(part.targets[e]))
                k == 0 || (code |= (0x01 << (k - 1)))
            end
            out[p] = FlowDirection{RingMask,UInt8}(code)
        end
    end
    return out
end
