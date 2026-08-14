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
    p === nothing && throw(BoundsError(r, c))
    return p
end

Base.getindex(r::IGeo7Raster, c::Cell) = parent(r)[_position(r, c)]
Base.setindex!(r::IGeo7Raster, value, c::Cell) =
    setindex!(parent(r), value, _position(r, c))
Base.checkbounds(::Type{Bool}, r::IGeo7Raster, c::Cell) =
    DGG.cellposition(_lookup(r), c) !== nothing

GM.neighbors(r::IGeo7Raster, c::Cell) = DGG.neighbors(_lookup(r), c)

function GM.outlets(r::IGeo7Raster)
    cells = _cells(r)
    complete = DGG.levelgrid(DGG.system(cells), DGG.level(cells))
    return Cell[
        c for c in cells
        if length(DGG.neighbors(cells, c)) < length(DGG.neighbors(complete, c))
    ]
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

const LDD_CODES = (0x05, 0x06, 0x09, 0x07, 0x04, 0x01, 0x03)
const D8D_CODES = (0x00, 0x01, 0x80, 0x20, 0x10, 0x08, 0x02)

GM.FlowDirection{GM.LDD}(d::RelativeIndex) =
    GM.FlowDirection{GM.LDD}(LDD_CODES[DGG.directioncode(d)+1])
GM.FlowDirection{GM.D8D}(d::RelativeIndex) =
    GM.FlowDirection{GM.D8D}(D8D_CODES[DGG.directioncode(d)+1])

const LDD_TO_IGEO7 = (0x05, 0xff, 0x06, 0x04, 0x00, 0x01, 0x03, 0xff, 0x02)
const D8D_BIT_TO_IGEO7 = (0x01, 0x06, 0xff, 0x05, 0x04, 0x03, 0xff, 0x02)

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
        neighbors = DGG.neighbors(
            DGG.levelgrid(DGG.IGeo7System(), DGG.level(center)),
            center,
        )
        code <= length(neighbors) ||
            throw(ArgumentError("direction $code does not exist at pentagon $center"))
        return neighbors[code] - center
    end
end

end # module
