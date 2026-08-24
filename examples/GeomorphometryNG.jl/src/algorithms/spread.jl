# ## Friction-distance spread
#
# The fourth algorithm family (decisions §6): a global traversal in a
# data-dependent order. Unlike the sweep family it needs *geometry* as well as
# topology — the cost of an edge is its length times the mean friction across it
# — so it is the case that pins down what the two traversal primitives can and
# cannot do.
#
# `neighbortable` is position-keyed but carries no distances;
# `neighborgeometry`/`eachneighbor` carries distances but is keyed by the
# array's own index type. This algorithm uses the second, with `cellkey` and
# `storageposition` bridging back to positions for its flat work arrays, and the
# whole loop behind a function barrier taking the hoisted geometry — the shape
# decisions §6 sketches for this family, and the shape §8 requires of it.

"""Friction-distance [`spread`](@ref) by priority-queue search (Tomlin, 1983)."""
struct Tomlin end

const SPREAD_NEEDS = (Index(), Value(), Distance())

"""
    spread(points, initial, friction; method=Tomlin(), limit=Inf)

Total friction distance from the source cells named by `points`, which is either
a mask shaped like `friction` (nonzero marks a source) or a vector of cell
indices. `initial` is the value at the sources: a scalar, or an array shaped
like `friction`. Unreached cells hold `limit`.

The cost of stepping between two cells is the distance between their centers
times the mean of their two friction values, so the result is in
`friction × length` units on both backends.
"""
spread(points, initial, friction, grid::AbstractGridSpec; method=Tomlin(), limit=Inf) =
    _spread(method, points, initial, friction, grid; limit)

_spread(method, points, initial, friction, grid::AbstractGridSpec; limit=Inf) =
    throw(ArgumentError(
        "spread has no rule for $(nameof(typeof(method))). The supported method is " *
        "Tomlin(). Geomorphometry's Eastman() is a pushbroom sweep in *storage* " *
        "order, which is a spatial sweep only on a rectilinear grid; its " *
        "FastSweeping() is an exported singleton with no method at all."))

# The source set, as storage positions.
_spreadseeds(points::AbstractVector{<:CartesianIndex}, grid, n) =
    [storageposition(grid, I) for I in points]
_spreadseeds(points::AbstractArray{Bool}, grid, n) =
    (_checkspreadsize(length(points), n); findall(vec(_data(points))))
_spreadseeds(points::AbstractArray{<:Real}, grid, n) =
    (_checkspreadsize(length(points), n); findall(>(0), vec(_data(points))))

_checkspreadsize(got, n) = got == n ||
    throw(ArgumentError("`points` has $got entries for $n cells"))

_seedvalues(initial::Real, seeds, n) = fill(Float64(initial), length(seeds))
function _seedvalues(initial, seeds, n)
    values = vec(_data(initial))
    _checkspreadsize(length(values), n)
    return [Float64(values[s]) for s in seeds]
end

function _spread(::Tomlin, points, initial, friction, grid::AbstractGridSpec; limit=Inf)
    frictions = vec(_data(friction))
    n = length(frictions)
    result = fill(Float64(limit), n)
    n == 0 && return reshape(result, size(_data(friction)))
    seeds = _spreadseeds(points, grid, n)
    seedvalues = _seedvalues(initial, seeds, n)
    settled = falses(n)
    queue = GM.FastPriorityQueue{Float64}(n)
    @inbounds for (s, v) in zip(seeds, seedvalues)
        result[s] = v
        queue[s] = v
    end
    geom = neighborgeometry(grid, NeighborRings(1), SPREAD_NEEDS)
    _spreadloop!(result, settled, queue, frictions, friction, grid, geom)
    return reshape(result, size(_data(friction)))
end

# THE FUNCTION BARRIER: `geom`'s concrete type is runtime-determined, exactly as
# `table`'s is in `_floodsweep!`.
function _spreadloop!(result, settled, queue, frictions, ras, grid, geom)
    while !isempty(queue)
        p = GM.dequeue!(queue)
        @inbounds settled[p] && continue
        @inbounds settled[p] = true
        base = @inbounds result[p]
        here = @inbounds frictions[p]
        for neighbor in eachneighbor(geom, ras, cellkey(ras, grid, p))
            q = storageposition(grid, neighbor.index)
            @inbounds settled[q] && continue
            candidate = base + (here + neighbor.value) * neighbor.distance / 2
            @inbounds candidate < result[q] || continue
            @inbounds result[q] = candidate
            queue[q] = candidate
        end
    end
    return result
end
