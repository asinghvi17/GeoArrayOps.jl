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

# The per-cell area of a DGG cell. The DiscreteGlobalGrids extension supplies
# the methods; the lazy vector below and the traversal's area fill are written
# against this one name.
function _cellareaof end

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
# `cellarea(::CellGrid, cell)` is in the DiscreteGlobalGrids extension.
