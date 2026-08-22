# # A tiny grid interface for Geomorphometry
#
# This proof of concept keeps DimensionalData vocabulary at the boundary.
# Numerical algorithms receive only `(ras, grid)` and ask the grid to map
# neighborhoods.  It is deliberately small: eight-neighbor steepest descent,
# a slope angle, and a compass flow direction.

using Rasters

# ## Grid objects
#
# `P == (xaxis, yaxis)`.  The dimensions and lookups are retained for adapters
# and rebuilding, but the algorithms below never inspect them.

abstract type AbstractGrid end

struct RectilinearGrid{P,D,L,S} <: AbstractGrid
    dims::D
    lookups::L
    spacing::S # signed semantic (x, y) spacing
end

function RectilinearGrid{P}(dims, lookups, spacing) where {P}
    P in ((1, 2), (2, 1)) || throw(ArgumentError("invalid axis map $P"))
    return RectilinearGrid{P,typeof(dims),typeof(lookups),typeof(spacing)}(
        dims, lookups, spacing,
    )
end

# A real DGG extension would keep a `Cells` dimension and `CellLookup` here.
# For the PoC, adjacency and cell centres stand in for that backend.  There is
# intentionally no axis parameter: cells must be the only axis.

struct CellGrid{D,L,A,C} <: AbstractGrid
    dims::D
    lookup::L
    adjacency::A
    centres::C
end

function CellGrid(dims, lookup, adjacency, centres)
    length(dims) == 1 || throw(ArgumentError("Cells must be the only axis"))
    length(lookup) == length(adjacency) == length(centres) ||
        throw(DimensionMismatch("cell topology fields must have equal lengths"))
    return CellGrid{typeof(dims),typeof(lookup),typeof(adjacency),typeof(centres)}(
        dims, lookup, adjacency, centres,
    )
end

# ## Adapters
#
# A matrix has the package's existing convention: axis 1 is X and axis 2 is Y.

function spatialparts(A::AbstractMatrix; spatialdims=nothing, spacing=nothing)
    isnothing(spatialdims) || spatialdims == (1, 2) ||
        throw(ArgumentError("a matrix uses spatialdims=(1, 2)"))
    steps = isnothing(spacing) ? (1.0, 1.0) : spacing
    grid = RectilinearGrid{(1, 2)}(axes(A), axes(A), steps)
    return A, grid
end

# The Raster adapter is the only code that understands DD dimensions and
# lookups.  The oracle is extension-friendly; DGG would add a Cells method.

isspatialdim(::Type) = false
isspatialdim(::Type{<:Rasters.XDim}) = true
isspatialdim(::Type{<:Rasters.YDim}) = true

function spatialparts(r::Raster; spatialdims=nothing, spacing=nothing)
    ndims(r) == 2 || throw(ArgumentError("the PoC expects one 2D surface"))
    selected = if isnothing(spatialdims)
        filter(d -> isspatialdim(typeof(d)), Rasters.dims(r))
    else
        Rasters.dims(r, spatialdims)
    end
    length(selected) == 2 || throw(ArgumentError("expected one X and one Y dimension"))

    xdim = only(filter(d -> d isa Rasters.XDim, selected))
    ydim = only(filter(d -> d isa Rasters.YDim, selected))
    P = Rasters.dimnum(r, (Rasters.XDim, Rasters.YDim))
    lookups = (Rasters.lookup(xdim), Rasters.lookup(ydim))
    steps = isnothing(spacing) ? (Float64(step(xdim)), Float64(step(ydim))) : spacing
    grid = RectilinearGrid{P}((xdim, ydim), lookups, steps)
    return r, grid
end

# ## The neighborhood vocabulary
#
# Logical offsets are always `(dx, dy)`: positive X is east and positive Y is
# north.  A rectilinear grid maps them to storage offsets using both the type
# parameter `P` and the sign of each coordinate step.

const ADJACENT_8 = (
    (-1, -1), (0, -1), (1, -1),
    (-1,  0),          (1,  0),
    (-1,  1), (0,  1), (1,  1),
)

function storageoffset(grid::RectilinearGrid{P}, (dx, dy)) where {P}
    xaxis, yaxis = P
    sx, sy = grid.spacing
    return CartesianIndex(ntuple(2) do axis
        axis == xaxis ? dx * Int(sign(sx)) : dy * Int(sign(sy))
    end)
end

function mapneighbors(f, ras, grid::RectilinearGrid, ::Type{T}) where {T}
    out = similar(ras, T)
    sx, sy = abs.(grid.spacing)
    for I in CartesianIndices(ras)
        neighbors = (
            let J = I + storageoffset(grid, offset),
                dx = offset[1], dy = offset[2]
                (; index=J, value=ras[J], distance=hypot(dx * sx, dy * sy),
                   bearing=mod(atand(dx * sx, dy * sy), 360.0))
            end
            for offset in ADJACENT_8
            if checkbounds(Bool, ras, I + storageoffset(grid, offset))
        )
        out[I] = f(I, ras[I], neighbors)
    end
    return out
end

function mapneighbors(f, ras::AbstractVector, grid::CellGrid, ::Type{T}) where {T}
    out = similar(ras, T)
    for i in eachindex(ras)
        x1, y1 = grid.centres[i]
        neighbors = (
            let x2 = grid.centres[j][1], y2 = grid.centres[j][2]
                (; index=j, value=ras[j], distance=hypot(x2 - x1, y2 - y1),
                   bearing=mod(atand(x2 - x1, y2 - y1), 360.0))
            end
            for j in grid.adjacency[i]
        )
        out[i] = f(i, ras[i], neighbors)
    end
    return out
end

# ## Two algorithms, one vocabulary

function steepest(value, neighbors)
    best_gradient = 0.0
    best_bearing = NaN
    for neighbor in neighbors
        gradient = (value - neighbor.value) / neighbor.distance
        if gradient > best_gradient
            best_gradient = gradient
            best_bearing = neighbor.bearing
        end
    end
    return best_gradient, best_bearing
end

_steepest_slope(ras, grid) = mapneighbors(ras, grid, Float64) do _, value, neighbors
    gradient, _ = steepest(value, neighbors)
    atand(gradient)
end

_flow_direction(ras, grid) = mapneighbors(ras, grid, Float64) do _, value, neighbors
    _, bearing = steepest(value, neighbors)
    bearing # clockwise from north; NaN marks a pit or edge outlet
end

function steepest_slope(input; kwargs...)
    ras, grid = spatialparts(input; kwargs...)
    return _steepest_slope(ras, grid)
end

function flow_direction(input; kwargs...)
    ras, grid = spatialparts(input; kwargs...)
    return _flow_direction(ras, grid)
end

steepest_slope(ras, grid::AbstractGrid) = _steepest_slope(ras, grid)
flow_direction(ras, grid::AbstractGrid) = _flow_direction(ras, grid)

# ## Run the experiment
#
# The same physical plane is stored once as X,Y and once as Y,X.  Y is also
# reverse ordered and the two cell widths differ, so an accidental storage-axis
# assumption is visible.

xs = 0.0:2.0:8.0
ys = 20.0:-5.0:5.0
surface(x, y) = 100.0 - 3.0x

data_xy = [surface(x, y) for x in xs, y in ys]
raster_xy = Raster(data_xy, (X(xs), Y(ys)))
raster_yx = Raster(permutedims(data_xy), (Y(ys), X(xs)))

_, grid_xy = spatialparts(raster_xy)
_, grid_yx = spatialparts(raster_yx)
@assert grid_xy isa RectilinearGrid{(1, 2)}
@assert grid_yx isa RectilinearGrid{(2, 1)}

slope_xy = steepest_slope(raster_xy)
slope_yx = steepest_slope(raster_yx; spatialdims=(X, Y))
direction_xy = flow_direction(raster_xy)
direction_yx = flow_direction(raster_yx)

@assert slope_xy isa Raster
@assert parent(slope_xy) ≈ permutedims(parent(slope_yx))
@assert all(isequal.(parent(direction_xy), permutedims(parent(direction_yx))))
@assert all(d -> isnan(d) || d ≈ 90.0, parent(direction_xy))

# A plain matrix takes the same algorithm path through a synthetic grid.
matrix_slope = steepest_slope(data_xy; spacing=(2.0, -5.0))
@assert matrix_slope ≈ parent(slope_xy)

# Finally, a tiny stand-in cell topology exercises the same algorithms.  The
# real DGG extension would construct this grid from its one Cells lookup.
cell_values = [10.0, 7.0, 8.0, 3.0]
cell_grid = CellGrid(
    (:Cells,),
    1:4,
    ((2, 3), (1, 4), (1, 4), (2, 3)),
    ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)),
)
cell_slope = steepest_slope(cell_values, cell_grid)
cell_direction = flow_direction(cell_values, cell_grid)
@assert cell_slope[1] ≈ atand(3.0)
@assert cell_direction[1] ≈ 90.0

println("PoC passed: X,Y; Y,X; Matrix; and one-axis CellGrid agree with the grid interface.")
