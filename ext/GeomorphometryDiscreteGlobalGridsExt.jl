module GeomorphometryDiscreteGlobalGridsExt

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import Rasters
using Rasters: Raster

const IGeo7Raster{T,D} =
    Raster{T,1,D} where {T,D<:Tuple{<:DGG.Cells{<:DGG.CellLookup{DGG.Z7Cell}}}}
const RelativeIndex = DGG.RelativeZ7Cell
const Cell = DGG.Z7Cell

_lookup(r::IGeo7Raster) = Rasters.lookup(r, DGG.Cells)
_cells(r::IGeo7Raster) = parent(_lookup(r))

# `CellVector` is lazy and its `in` method searches compressed position windows
# rather than scanning cell ids.
Base.eachindex(r::IGeo7Raster) = _cells(r)

function _position(r::IGeo7Raster, c::Cell)
    p = DGG.cellposition(_lookup(r), c)
    # Tuple-wrapped: `showerror` iterates the index field, and a `Cell` is not iterable.
    isnothing(p) && throw(BoundsError(r, (c,)))
    return p
end

Base.getindex(r::IGeo7Raster, c::Cell) = parent(r)[_position(r, c)]
Base.setindex!(r::IGeo7Raster, value, c::Cell) =
    setindex!(parent(r), value, _position(r, c))
Base.checkbounds(::Type{Bool}, r::IGeo7Raster, c::Cell) =
    DGG.cellposition(_lookup(r), c) !== nothing

GM.neighbors(r::IGeo7Raster, c::Cell) = DGG.neighbors(_lookup(r), c)

# --- the neighborhood forms -------------------------------------------------

GM.neighbors(r::IGeo7Raster) = DGG.neighbors(_cells(r))

function GM.mapneighbors(f::F, r::IGeo7Raster; order = nothing,
        threaded = true) where {F}
    out = DGG.mapneighbors(f, _cells(r), parent(r);
        order = order === nothing ? DGG.StorageOrder() : order, threaded)
    return Rasters.rebuild(r; data = out)
end

# The rim: cells whose clipped ring is shorter than their complete-level
# degree. One threaded sweep instead of two resolved neighbor sets per cell.
function GM.outlets(r::IGeo7Raster)
    cells = _cells(r)
    complete = DGG.levelgrid(DGG.system(cells), DGG.level(cells))
    rim = DGG.mapneighbors(
        (c, nbrs) -> length(nbrs) < DGG.neighborcount(complete, DGG.cellid(c)),
        cells; threaded = true)
    return Cell[cells[p] for p in eachindex(rim) if rim[p]]
end

# --- the priority-flood queue phase in position space -----------------------

# Settle every cell from the rim inward in ascending elevation — the queue
# phase `flowaccumulation!` and `height_above_nearest_drainage` share — with
# positions as the only currency: the queue is keyed on Int positions with a
# dense locator (no hashing), each cell's clipped ring is one row of a
# `HaloTable` built once, and the neighbor that settles a cell is recorded as
# a downstream POSITION in `down`. `order` receives the settle sequence;
# positions closed on entry are skipped exactly as the cell-indexed loop
# skips them, and every rim cell is seeded regardless of that mask.
function _settle!(order::Vector{Int64}, down::Vector{Int}, closedv, zv, cv,
        table)
    n = length(cv)
    complete = DGG.levelgrid(DGG.system(cv), DGG.level(cv))
    open = GM.FastPriorityQueue{eltype(zv)}(n)
    @inbounds for p in 1:n
        if length(table[p]) < DGG.neighborcount(complete, cv[p])
            GM.enqueue!(open, p, zv[p])
            closedv[p] = true
        end
    end
    i = 1
    @inbounds while !isempty(open)
        p = GM.dequeue!(open)
        order[i] = p
        i += 1
        for q in table[p]
            closedv[q] && continue
            closedv[q] = true
            down[q] = p
            GM.enqueue!(open, q, zv[q])
        end
    end
    return nothing
end

# The relative-cell directions the settle implies, materialized once per cell
# from the downstream positions; zero where nothing settled the cell.
function _directions(dem::IGeo7Raster, cv, down)
    zrel = first(cv) - first(cv)
    dir = similar(dem, typeof(zrel))
    dirv = parent(dir)
    @inbounds for p in eachindex(down)
        dirv[p] = down[p] == 0 ? zrel : cv[down[p]] - cv[p]
    end
    return dir
end

function GM.flowaccumulation(dem::IGeo7Raster,
        closed = fill!(similar(dem, Bool), false);
        method = GM.DInf(), cellsize = GM.cellsize(dem))
    cv = _cells(dem)
    acc = similar(dem, Float32)
    accv = parent(acc)
    @inbounds for p in eachindex(accv)
        accv[p] = DGG.IGeo7.cell_area(DGG.rawid(cv[p]))
    end
    return GM.flowaccumulation!(dem, acc, copy(closed); method, cellsize)
end

function GM.flowaccumulation!(dem::IGeo7Raster, acc::AbstractArray{<:Real},
        closed = fill!(similar(dem, Bool), false);
        method = GM.DInf(), cellsize = GM.cellsize(dem))
    cv = _cells(dem)
    closedv = parent(closed)
    order = ones(Int64, length(cv) - count(closedv))
    down = zeros(Int, length(cv))
    table = DGG.HaloTable(cv)
    _settle!(order, down, closedv, parent(dem), cv, table)
    return _postsettle(method, dem, acc, order, down, cv, cellsize)
end

# D8 needs no cell arithmetic after the settle: accumulation follows the
# downstream positions and the persisted directions — cell-relative, exactly
# as the generic path writes them — are encoded once per cell.
function _postsettle(::GM.D8, dem, acc, order, down, cv, cellsize)
    accv = parent(acc)
    @inbounds for j in reverse(eachindex(order))
        p = order[j]
        d = down[p]
        d == 0 && continue
        accv[d] += accv[p]
    end
    zrel = first(cv) - first(cv)
    output = similar(dem, GM.FlowDirection{GM.LDD, UInt8})
    outv = parent(output)
    @inbounds for p in eachindex(down)
        rel = down[p] == 0 ? zrel : cv[down[p]] - cv[p]
        outv[p] = GM.FlowDirection{GM.LDD}(rel)
    end
    return acc, output
end

# The DInf and FD8 accumulate phases read the DEM's geometry around each
# cell, so they keep the generic implementation; the settle's product is
# materialized into the direction raster they consume.
function _postsettle(method::GM.FlowDirectionMethod, dem, acc, order, down,
        cv, cellsize)
    dir = _directions(dem, cv, down)
    dirs = GM._accumulate!(method, acc, order, dir, cv, dem, cellsize)
    return acc, dirs
end

function GM.height_above_nearest_drainage(dem::IGeo7Raster;
        method = GM.D8(), cellsize = GM.cellsize(dem), threshold = 100)
    cv = _cells(dem)
    n = length(cv)
    output = zero(dem)
    acc = similar(dem, Float32)
    accv = parent(acc)
    @inbounds for p in 1:n
        accv[p] = DGG.IGeo7.cell_area(DGG.rawid(cv[p]))
    end
    closed = fill!(similar(dem, Bool), false)
    order = ones(Int64, n)
    down = zeros(Int, n)
    table = DGG.HaloTable(cv)
    _settle!(order, down, parent(closed), parent(dem), cv, table)
    dir = _directions(dem, cv, down)
    flowdirs = GM._accumulate!(method, acc, order, dir, cv, dem, cellsize)
    stream_mask = acc .>= threshold
    GM._burn_streams!(stream_mask, order, dir, cv)
    GM._fill_flow_depressions!(acc, order, flowdirs, cv, dem, cellsize)
    GM._hand!(output, order, flowdirs, cv, acc, stream_mask, cellsize)
    return output
end

function GM.cellarea(r::IGeo7Raster, c::Cell; cellsize=nothing)
    _position(r, c)
    return DGG.IGeo7.cell_area(DGG.rawid(c))
end

function GM.celldistance(
    r::IGeo7Raster,
    from::Cell,
    to::Cell;
    cellsize=nothing,
)
    _position(r, from)
    _position(r, to)
    grid = DGG.levelgrid(DGG.system(_cells(r)), DGG.level(from))
    a = DGG.cell_centroid(grid, from)
    b = DGG.cell_centroid(grid, to)
    angle = acos(clamp(a[1] * b[1] + a[2] * b[2] + a[3] * b[3], -1.0, 1.0))
    return angle * DGG.IGeo7.R_AUTHALIC
end

function GM.cellbearing(r::IGeo7Raster, from::Cell, to::Cell)
    _position(r, from)
    _position(r, to)
    from == to && return 0.0
    grid = DGG.levelgrid(DGG.system(_cells(r)), DGG.level(from))
    a = DGG.cell_centroid(grid, from)
    b = DGG.cell_centroid(grid, to)
    lon1, lon2 = atan(a[2], a[1]), atan(b[2], b[1])
    lat1, lat2 = asin(clamp(a[3], -1.0, 1.0)), asin(clamp(b[3], -1.0, 1.0))
    delta_lon = lon2 - lon1
    east = sin(delta_lon) * cos(lat2)
    north = cos(lat1) * sin(lat2) -
            sin(lat1) * cos(lat2) * cos(delta_lon)
    return mod(rad2deg(atan(east, north)), 360.0)
end

const IGEO7_TO_LDD = (0x05, 0x06, 0x09, 0x07, 0x04, 0x01, 0x03)
const IGEO7_TO_D8D = (0x00, 0x01, 0x80, 0x20, 0x10, 0x08, 0x02)

GM.FlowDirection{GM.LDD}(d::RelativeIndex) =
    GM.FlowDirection{GM.LDD}(IGEO7_TO_LDD[DGG.directioncode(d)+1])
GM.FlowDirection{GM.D8D}(d::RelativeIndex) =
    GM.FlowDirection{GM.D8D}(IGEO7_TO_D8D[DGG.directioncode(d)+1])

function _inverse_code(codes, value)
    index = findfirst(==(UInt8(value)), codes)
    return isnothing(index) ? 0xff : UInt8(index - 1)
end

const LDD_TO_IGEO7 = ntuple(value -> _inverse_code(IGEO7_TO_LDD, value), 9)
const D8D_BIT_TO_IGEO7 =
    ntuple(bit -> _inverse_code(IGEO7_TO_D8D, 1 << (bit - 1)), 8)

@inline function _directioncode(direction::GM.FlowDirection{GM.LDD})
    value = Int(direction)
    1 <= value <= length(LDD_TO_IGEO7) ||
        throw(ArgumentError("invalid LDD direction $value"))
    code = @inbounds LDD_TO_IGEO7[value]
    code != 0xff ||
        throw(ArgumentError("LDD direction $value has no IGeo7 equivalent"))
    return code
end

@inline function _directioncode(direction::GM.FlowDirection{GM.D8D})
    iszero(direction) && return 0x00
    value = Int(direction)
    ispow2(value) ||
        throw(ArgumentError("D8D direction must be decomposed before conversion"))
    bit = trailing_zeros(value) + 1
    bit <= length(D8D_BIT_TO_IGEO7) ||
        throw(ArgumentError("D8D direction $value has no IGeo7 equivalent"))
    code = @inbounds D8D_BIT_TO_IGEO7[bit]
    code != 0xff ||
        throw(ArgumentError("D8D direction $value has no IGeo7 equivalent"))
    return code
end

function GM.decompose(
    ::Type{RelativeIndex},
    direction::GM.FlowDirection,
    center::Cell,
)
    return Iterators.map(GM.decompose(direction)) do component
        code = _directioncode(component)
        iszero(code) && return RelativeIndex(center, 0)
        adjacent = DGG.neighbors(
            DGG.levelgrid(DGG.IGeo7System(), DGG.level(center)),
            center,
        )
        code <= length(adjacent) ||
            throw(ArgumentError("direction $code does not exist at pentagon $center"))
        return adjacent[code] - center
    end
end

end # module
