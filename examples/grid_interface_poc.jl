# # A tiny grid interface for Geomorphometry
#
# This proof of concept keeps DimensionalData vocabulary at the boundary.
# Numerical algorithms receive only `(ras, grid)` and ask the grid to map
# neighborhoods.  It is deliberately small: eight-neighbor steepest descent,
# a slope angle, and a compass flow direction.

using Rasters
import DiscreteGlobalGrids as DGG
import Stencils

# ## Grid objects
#
# `P == (xaxis, yaxis)`.  The dimensions and lookups are retained for adapters
# and rebuilding, but the algorithms below never inspect them.

abstract type AbstractGridSpec end # no relation to DiscreteGlobalGrids.AbstractGrid ;)

struct RectilinearGrid{P,D,L,S} <: AbstractGridSpec
    dims::D
    lookups::L
    spacing::S # signed semantic (x, y) spacing
end

function RectilinearGrid{P}(dims, lookups, spacing) where {P}
    P in ((1, 2), (2, 1)) || throw(ArgumentError("invalid axis map $P"))
    all(!iszero, spacing) || throw(ArgumentError("cell spacing must be nonzero"))
    return RectilinearGrid{P,typeof(dims),typeof(lookups),typeof(spacing)}(
        dims, lookups, spacing,
    )
end

# A cell grid retains the adapter details needed to ask DGG about topology and
# geometry.  It has no axis parameter: `Cells` must be the only array axis.

struct CellGrid{D,L,G} <: AbstractGridSpec
    dims::D
    lookup::L
    levelgrid::G
end

function CellGrid(dims)
    length(dims) == 1 || throw(ArgumentError("Cells must be the only axis"))
    dim = only(dims)
    lookup = Rasters.lookup(dim)
    lookup isa DGG.CellLookup ||
        throw(ArgumentError("Dimension must contain a DGG.CellLookup; found $(typeof(lookup))"))
    levelgrid = DGG.levelgrid(DGG.system(lookup), DGG.level(lookup))
    return CellGrid{typeof(dims),typeof(lookup),typeof(levelgrid)}(dims, lookup, levelgrid)
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
# lookups.  The oracle is extension-friendly: each integration registers its
# own dimension types here.

isspatialdim(::Type) = false
isspatialdim(::Type{<:Rasters.XDim}) = true
isspatialdim(::Type{<:Rasters.YDim}) = true
isspatialdim(::Type{<:DGG.Cells}) = true

_astuple(x::Tuple) = x
_astuple(x) = (x,)

function spatialparts(r::Raster; spatialdims=nothing, spacing=nothing)
    found_spatialdims = if isnothing(spatialdims)
        filter(d -> isspatialdim(typeof(d)), Rasters.dims(r))
    else
        Rasters.dims(r, spatialdims)
    end
    selected = _astuple(found_spatialdims)

    if length(selected) == 1 && lookup(only(selected)) isa DGG.CellLookup
        ndims(r) == 1 || throw(ArgumentError("Cells must be the only array axis"))
        isnothing(spacing) || throw(ArgumentError("CellGrid derives geometry from DGG"))
        return r, CellGrid(selected)
    end

    ndims(r) == 2 || throw(ArgumentError("a rectilinear surface must be 2D"))
    length(selected) == 2 ||
        throw(ArgumentError("expected one X and one Y dimension, or one Cells dimension"))

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
# `NORTH_UP_NEIGHBORS` is expressed in logical `(dx, dy)` coordinates: positive
# X is east and positive Y is north.  The rectilinear implementation turns it
# into a storage-ordered Stencils.jl stencil.  Algorithms never see that work.

const NORTH_UP_NEIGHBORS = (
    northwest=(-1,  1), north=(0,  1), northeast=(1,  1),
    west=(-1,  0),                       east=(1,  0),
    southwest=(-1, -1), south=(0, -1), southeast=(1, -1),
)

function storageoffset(grid::RectilinearGrid{P}, (dx, dy)) where {P}
    xaxis, yaxis = P
    sx, sy = grid.spacing
    return ntuple(2) do axis
        axis == xaxis ? dx * Int(sign(sx)) : dy * Int(sign(sy))
    end
end

northupstencil(grid::RectilinearGrid) =
    Stencils.NamedStencil(map(offset -> storageoffset(grid, offset), NORTH_UP_NEIGHBORS))

function mapneighbors(f, ras, grid::RectilinearGrid)
    sx, sy = abs.(grid.spacing)
    stencil = northupstencil(grid)
    return Stencils.mapstencil(
        stencil, ras, CartesianIndices(ras); boundary=Stencils.Remove(missing)
    ) do hood, I
        neighbors = (
            let dx = logical_offset[1], dy = logical_offset[2]
                (; value, distance=hypot(dx * sx, dy * sy),
                   bearing=mod(atand(dx * sx, dy * sy), 360.0))
            end
            for (value, logical_offset) in zip(
                Stencils.neighbors(hood),
                values(NORTH_UP_NEIGHBORS),
            )
            if !ismissing(value)
        )
        f(I, Stencils.center(hood), neighbors)
    end
end

# DGG's lookup maps a storage position directly to neighboring storage
# positions.  It preserves the topology's ring order and clips neighbors to
# lookup membership.  Geomorphometry enriches those positions with the same
# geometry record used by the rectilinear implementation.

const AUTHALIC_RADIUS_M = 6_371_007.180918475

function cellmetrics(grid::CellGrid, from_cell, to_cell)
    from = DGG.cell_centroid(grid.levelgrid, from_cell)
    to = DGG.cell_centroid(grid.levelgrid, to_cell)
    central_angle = acos(clamp(sum(from[k] * to[k] for k in 1:3), -1.0, 1.0))

    lon1, lon2 = atan(from[2], from[1]), atan(to[2], to[1])
    lat1, lat2 = asin(clamp(from[3], -1.0, 1.0)), asin(clamp(to[3], -1.0, 1.0))
    delta_lon = lon2 - lon1
    east = sin(delta_lon) * cos(lat2)
    north = cos(lat1) * sin(lat2) -
            sin(lat1) * cos(lat2) * cos(delta_lon)
    bearing = mod(rad2deg(atan(east, north)), 360.0)
    return central_angle * AUTHALIC_RADIUS_M, bearing
end

# DGG owns cell traversal and output inference.  Geomorphometry only translates
# its positioned handles into the same records the rectilinear path produces.
function mapneighbors(f, ras::AbstractVector, grid::CellGrid)
    return DGG.mapneighbors(ras; pass=DGG.Neighbors()) do cell, adjacent
        i = DGG.cellposition(cell)
        neighbors = (
            let (distance, bearing) =
                    cellmetrics(grid, DGG.cellid(cell), DGG.cellid(neighbor))
                (; value=ras[neighbor], distance, bearing)
            end
            for neighbor in adjacent
        )
        f(i, ras[cell], neighbors)
    end
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

function _steepest_slope(_, value, neighbors)
    gradient, _ = steepest(value, neighbors)
    return atand(gradient)
end

function _flow_direction(_, value, neighbors)
    _, bearing = steepest(value, neighbors)
    return bearing # clockwise from north; NaN marks a pit or edge outlet
end

steepest_slope(input; kwargs...) = steepest_slope(spatialparts(input; kwargs...)...)
flow_direction(input; kwargs...) = flow_direction(spatialparts(input; kwargs...)...)

steepest_slope(ras, grid::AbstractGridSpec) = mapneighbors(_steepest_slope, ras, grid)
flow_direction(ras, grid::AbstractGridSpec) = mapneighbors(_flow_direction, ras, grid)

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

# Named fields have the same geographic meaning in either storage layout.
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

# Finally, a real DGG raster takes the same algorithm path.  Elevation rises
# with the centroid's Z coordinate, giving most cells a downhill neighbor.
level = 2
dgg = DGG.levelgrid(DGG.IGeo7System(), level)
cell_lookup = DGG.CellLookup(dgg)
cell_values = [10_000.0 * DGG.cell_centroid(dgg, cell)[3] for cell in cell_lookup]
cell_raster = Raster(cell_values, (DGG.Cells(cell_lookup),))

_, cell_grid = spatialparts(cell_raster; spatialdims=DGG.Cells)
@assert cell_grid isa CellGrid
@assert length(cell_grid.dims) == 1
@assert cell_grid.lookup === cell_lookup

cell_slope = steepest_slope(cell_raster)
cell_direction = flow_direction(cell_raster)
@assert cell_slope isa Raster
@assert cell_direction isa Raster
@assert all(>=(0.0), parent(cell_slope))
@assert any(isfinite, parent(cell_direction))
@assert all(d -> isnan(d) || 0.0 <= d < 360.0, parent(cell_direction))

println("PoC passed: X,Y; Y,X; Matrix; and a real one-axis DGG CellGrid use one interface.")
