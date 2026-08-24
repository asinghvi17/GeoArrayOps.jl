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

#
# The `CellGrid(celldim)` constructor is defined by the DiscreteGlobalGrids
# extension: the struct is pure storage, but building one means asking DGG for
# the cell collection's system and level.

struct CellGrid{C,G} <: AbstractGridSpec
    manifold::Spherical{Float64}
    cells::C
    levelgrid::G
end

_cellmanifold(::Nothing) = Spherical(; radius=AUTHALIC_RADIUS_M)
_cellmanifold(m::Spherical) = Spherical(; radius=Float64(m.radius))
_cellmanifold(m::Manifold) = throw(ArgumentError(
    "a DGG cell grid is spherical; got $(typeof(m))"))

manifold(grid::CellGrid) = grid.manifold

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

# The *named* window, dilated by `radius`. This is the stencil the windowed
# derivative family reads by geographic name: `w.northeast` is northeast at any
# radius and in either storage order. It is deliberately not `northupstencil`
# with `NeighborRings(radius)` — that is the cumulative ring set, which loses the
# names; a dilated eight-neighborhood keeps them.
function northupwindow(grid::RectilinearGrid, radius::Int=1)
    radius >= 1 || throw(ArgumentError("window radius must be positive; got $radius"))
    return Stencils.NamedStencil(
        map(o -> storageoffset(grid, o .* radius), NORTH_UP_NEIGHBORS))
end

# These methods expose stable cell identifiers and convert them to storage
# indices without requiring algorithms to know the grid representation.

cellindices(ras, ::RectilinearGrid) = CartesianIndices(axes(ras))

storagekey(::RectilinearGrid, I) = I

# The `CellGrid` methods of both are in the DiscreteGlobalGrids extension.

# ## Positions and keys
#
# `cellindices` streams keys in storage-position order, which serves a
# whole-grid sweep. A *queue-driven* traversal that also needs geometry pops a
# storage position and must turn it back into the key `eachneighbor` and the
# array accept — and must turn a record's `Index` back into a position to reach
# its own flat work arrays. Neither direction is expressible through
# `cellindices` (a `CellGrid` yields a generator, not an indexable collection),
# so the bijection is named here.
#
# `neighbortable` is the other traversal primitive and needs neither, because it
# is position-keyed end to end — but it carries no geometry. An algorithm that
# needs *both* positions and distances (cost distance, DInf, FD8) is exactly the
# case these two methods exist for.

cellkey(ras, grid::RectilinearGrid, p::Int) = @inbounds CartesianIndices(gridsize(grid))[p]
storageposition(grid::RectilinearGrid, I::CartesianIndex) =
    @inbounds LinearIndices(gridsize(grid))[I]

# `cellkey(ras, ::CellGrid, p)` and `storageposition(::CellGrid, cell)` are in
# the DiscreteGlobalGrids extension.
