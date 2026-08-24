# ## The sweep family
#
# A priority flood is random access in a data-dependent order, so this family is
# written against [`neighbortable`](@ref) — the slot-stable, position-keyed
# traversal primitive in `interface/neighbortable.jl` — and not against the
# record API, which serves whole-grid sweeps in storage order.
#
# Everything below is keyed by storage position (`1:ncells`) on both backends,
# never by cell id. The rectilinear and cell paths differ only in
# `neighbortable`, `_fillcellareas!` and the direction codec (`flowdir.jl`); the
# flood, the settling pass and the accumulation pass are one implementation each.

# The sweep family's D8 rule is "lowest neighbor", not "steepest gradient", so it
# reads no geometry whatsoever — `neighbortable` serves it instead of the record
# API. This is what a multi-direction method would have to request, and why
# `DInf()` and `FD8()` are refused there rather than silently downgraded.
const MULTIDIRECTION_NEEDS = (Index(), Value(), Distance(), Bearing())

function _boundaryseeds(tbl, closedv)
    seeds = Int[]
    @inbounds for p in 1:length(tbl)
        (closedv[p] || !isboundary(tbl, p)) && continue
        push!(seeds, p)
    end
    return seeds
end

# ### The flood
#
# `order` is the pop sequence and `down[p]` the position that first reached `p`.
# On the surface the flood implies — the depression-filled one — `down[p]` is a
# downstream neighbor, which is Barnes' "flow directions from Priority-Flood":
# routing out of depressions is a property of the pop order, not of a separate
# filling pass.
struct FloodSweep{P<:Integer}
    order::Vector{P} # Pop sequence, truncated to the number of cells visited
    down::Vector{P}  # Downstream position; 0 at an outlet or an unvisited cell
    nseeds::Int
end

_closedmask(::Nothing, n) = falses(n)
function _closedmask(closed, n)
    mask = collect(Bool, vec(_data(closed)))
    length(mask) == n ||
        throw(ArgumentError("`closed` has $(length(mask)) entries for $n cells"))
    return mask
end

"""
    floodsweep(ras, grid; closed=nothing, seeds=nothing, table=neighbortable(grid))

Priority-flood traversal from the domain boundary inward, taking the lowest
queued value each step. Returns a `FloodSweep`.

`closed` cells are never enqueued, never visited and never receive or emit flow.
`seeds` overrides the boundary seed set with real hydrological outlets. `table`
is threaded through so a composite operation builds the adjacency once.

Serial by nature: the queue order is the algorithm. Allocates `order`, `down`, a
`closed` scratch and the queue, and nothing else that is O(ncells).
"""
function floodsweep(ras, grid::AbstractGridSpec; closed=nothing, seeds=nothing,
        table=neighbortable(grid, NeighborRings(1)))
    z = vec(_data(ras))
    n = length(z)
    closedv = _closedmask(closed, n)
    # Positions index the grid, so `Int32` carries any grid below 2.1e9 cells at
    # half the residency. The branch is answered once, above two barriers.
    return n <= typemax(Int32) ? _floodsweep(Int32, z, closedv, table, seeds) :
        _floodsweep(Int64, z, closedv, table, seeds)
end

function _floodsweep(::Type{P}, z, closedv, table, seeds) where {P}
    seedv = seeds === nothing ? _boundaryseeds(table, closedv) : collect(Int, seeds)
    isempty(seedv) && (seedv = _noboundaryseeds(z, closedv))
    order = Vector{P}(undef, length(z) - count(closedv))
    down = zeros(P, length(z))
    nvisited, nseeded = _floodsweep!(order, down, closedv, z, table, seedv)
    resize!(order, nvisited) # Unreachable cells leave a tail nobody may read
    return FloodSweep(order, down, nseeded)
end

# THE FUNCTION BARRIER. `table`'s concrete type is runtime-determined, so this
# loop cannot live at the hoist site: inlining a traversal there cost 35x time
# and ~170x allocations in the v2 measurements.
function _floodsweep!(order::Vector{P}, down::Vector{P}, closedv, z, tbl, seeds) where {P}
    # Array-backed, keyed by position: no hash and no `Dict` per enqueue.
    open = GM.FastPriorityQueue{eltype(z)}(length(z))
    nseeded = 0
    @inbounds for s in seeds
        p = Int(s)
        closedv[p] && continue # A seed the caller also closed is not a seed
        closedv[p] = true
        nseeded += 1
        GM.enqueue!(open, p, z[p])
    end
    nord = length(order)
    i = 0
    @inbounds while !isempty(open)
        p = GM.dequeue!(open)
        i += 1
        i <= nord || _sweepoverrun(i, nord)
        order[i] = P(p)
        for (_, q) in slots(tbl, p)
            (q == 0 || closedv[q]) && continue
            closedv[q] = true
            down[q] = P(p)
            GM.enqueue!(open, q, z[q])
        end
    end
    return i, nseeded
end

@noinline _sweepoverrun(i, n) = error(
    "the flood visited $i cells but `order` holds $n: the seed set and `closed` " *
    "disagree about how many cells are reachable")

# A complete level grid has no boundary, so there is no outlet to seed and the
# honest answers are "seed the lowest cell" or "refuse". Seeding the global
# minimum keeps the call total and conserves area; what it must not do is return
# the silent `acc == cellarea` answer an empty queue produces.
function _noboundaryseeds(z, closedv)
    isempty(z) && return Int[]
    minima = _globalminima(z, closedv)
    if isempty(minima)
        @warn "no seed cells: every cell is closed or non-finite. Nothing is visited."
    else
        @warn "no domain boundary: seeding the global minimum. A complete level grid " *
              "has no outlet — pass `seeds=` for a real one." nseeds = length(minima)
    end
    return minima
end

function _globalminima(z::AbstractVector{T}, closedv) where {T}
    found = false
    best = zero(T)
    @inbounds for p in eachindex(z)
        (closedv[p] || !isfinite(z[p])) && continue
        if !found || z[p] < best
            best = z[p]
            found = true
        end
    end
    found || return Int[]
    return [p for p in eachindex(z) if !closedv[p] && z[p] == best]
end

# ### `settle`: the depression-filled surface

"""
    settle(ras, grid; closed=nothing, sweep=floodsweep(ras, grid; closed))

Depression-filled ("settled") elevations: `settled[p] = max(z[p], settled[down[p]])`,
with outlets keeping their own elevation. This is the minimax (spill-point) fill,
so it is value-identical to Geomorphometry's `filldepressions` — flats are filled
to exactly the spill elevation, with no epsilon.

Non-finite input stays non-finite and does not propagate upstream; a cell the
flood never visited is `NaN`. Integer input is promoted to `float(eltype)`.
Serial: it walks the pop order forward.
"""
function settle(ras, grid::AbstractGridSpec; closed=nothing,
        sweep=floodsweep(ras, grid; closed))
    return reshape(_settledvalues(ras, sweep), size(_data(ras)))
end

# The flat settled surface, which is what the traversals downstream of `settle`
# — `height_above_nearest_drainage`, `depression_depth` — actually consume.
function _settledvalues(ras, sweep)
    z = _floatvalues(vec(_data(ras)))
    out = similar(z)
    _settlepass!(out, z, sweep.order, sweep.down)
    return out
end

_floatvalues(z::AbstractVector{<:AbstractFloat}) = z
_floatvalues(z::AbstractVector) = float.(z) # So that "unvisited" is expressible

function _settlepass!(out, z, order, down)
    fill!(out, NaN)
    @inbounds for p in order # Forward: a parent always pops before its children
        zp = z[p]
        d = down[p]
        out[p] = (d != 0 && isfinite(zp) && isfinite(out[d])) ? max(zp, out[d]) : zp
    end
    return out
end

# ### Accumulation

"""
    flowaccumulation(ras, grid; method=D8(), closed=nothing) -> (acc, directions)

`acc[c]` is the total upstream area draining through `c`, in the grid's area
units, including `c`'s own cell area. `directions` is the matching
[`flowdirection`](@ref) output.

`acc` is `Float64`: at level 13 the largest accumulation is ~2.7e9 m² built from
16.2M additions of ~526 m² each, and a `Float32` mantissa resolves only ~161 m²
there — a third of a cell lost per addition near the outlet.

The flood and the downhill pass are serial; the cell-area fill and the direction
encoding are threaded over contiguous chunks.
"""
flowaccumulation(ras, grid::AbstractGridSpec; method=D8(), closed=nothing) =
    _flowaccumulation(method, ras, grid, closed)

function _flowaccumulation(::D8, ras, grid::AbstractGridSpec, closed)
    table = neighbortable(grid, NeighborRings(1)) # Built once, threaded through
    sweep = floodsweep(ras, grid; closed, table)
    acc = Vector{Float64}(undef, length(_data(ras)))
    _fillcellareas!(acc, grid)
    _accumulatedown!(acc, sweep.order, sweep.down)
    dirs = flowdirection(ras, grid; sweep, table)
    return reshape(acc, size(_data(ras))), dirs
end

_flowaccumulation(method::FlowDirectionMethod, ras, grid, closed) = throw(ArgumentError(
    "flowaccumulation has no rule for $(nameof(typeof(method))). The supported " *
    "methods are D8(), DInf() and FD8()."))

# Reverse pop order is a topological order: a child always pops after its parent,
# so reversing guarantees a cell's own inflow is complete before it drains.
function _accumulatedown!(acc, order, down)
    @inbounds for j in reverse(eachindex(order))
        p = order[j]
        d = down[p]
        d == 0 && continue
        acc[d] += acc[p]
    end
    return acc
end

# `cellarea` supplies the only unit-carrying input to the whole family. Nothing
# else in it needs geometry: no distances, no bearings, no manifold dispatch.
_fillcellareas!(acc, grid::RectilinearGrid) = _copychunks!(acc, cellarea(grid))
_fillcellareas!(acc, grid::CellGrid) =
    _cellareachunks!(acc, grid.cells, grid.levelgrid, manifold(grid).radius)

function _copychunks!(acc, areas)
    Threads.@threads for r in _chunkranges(length(acc))
        @inbounds for p in r
            acc[p] = areas[p]
        end
    end
    return acc
end

# Each chunk walks `p` ascending, which is what a `CellVector` window cursor wants.
function _cellareachunks!(acc, cells, levelgrid, radius)
    Threads.@threads for r in _chunkranges(length(acc))
        @inbounds for p in r
            acc[p] = _cellareaof(levelgrid, cells[p], radius)
        end
    end
    return acc
end

# ### Depressions
#
# Two readings of the same difference: `settle` minus the surface, per cell and
# integrated. Both are grid-generic because `settle` and `cellarea` are.

"""
    depression_depth(dem; closed=nothing, filled=settle(dem))

How deep each cell sits below the depression-filled surface. Zero outside
depressions.
"""
function depression_depth(ras, grid::AbstractGridSpec; closed=nothing,
        sweep=floodsweep(ras, grid; closed), filled=nothing)
    surface = isnothing(filled) ? _settledvalues(ras, sweep) :
        _floatvalues(vec(_data(filled)))
    z = _floatvalues(vec(_data(ras)))
    return reshape(surface .- z, size(_data(ras)))
end

"""
    depression_volume(dem; closed=nothing) -> Real

The total volume of every depression in the grid: the depth integrated against
[`cellarea`](@ref), so it is in the grid's length³ units on both backends.
"""
function depression_volume(ras, grid::AbstractGridSpec; closed=nothing,
        sweep=floodsweep(ras, grid; closed), filled=nothing)
    depth = vec(depression_depth(ras, grid; closed, sweep, filled))
    areas = Vector{Float64}(undef, length(depth))
    _fillcellareas!(areas, grid)
    return _volumesum(depth, areas)
end

function _volumesum(depth, areas)
    total = 0.0
    @inbounds for p in eachindex(depth)
        isfinite(depth[p]) || continue
        total += depth[p] * areas[p]
    end
    return total
end

# ### Accumulation-derived indices
#
# Three one-liners over `slope` and `flowaccumulation`. They are grid-generic
# because both inputs are; `slopemethod` is the seam that picks a windowed or a
# record estimator per backend, exactly as `slope` does on its own.

for (name, expr, doc) in (
        (:topographic_wetness_index, :(log(a / tand(s))),
            "Topographic Wetness Index: `log(accumulation / tan(slope))`."),
        (:stream_power_index, :(log(a * tand(s))),
            "Stream Power Index: `log(accumulation * tan(slope))`."),
        (:drainage_potential, :(sind(s) / log1p(a)),
            "How well a cell drains: `sin(slope) / log1p(accumulation)`."))
    @eval begin
        """
            $($(QuoteNode(name)))(dem; method=D8(), slopemethod=..., closed=nothing)

        $($doc)
        """
        function $name(ras, grid::AbstractGridSpec; method=D8(), closed=nothing,
                slopemethod=defaultmethod(slope, grid))
            s = _data(slope(ras, grid; method=slopemethod))
            acc, _ = flowaccumulation(ras, grid; method, closed)
            return map((a, s) -> $expr, acc, s)
        end
    end
end
const TWI = topographic_wetness_index
const SPI = stream_power_index

# ### Height above nearest drainage
#
# One flood, one adjacency table, `down` positions throughout. The stream mask
# is burned *down* the tree in reverse pop order, so a stream stays a stream all
# the way to its outlet; the height is then integrated *up* the tree in forward
# pop order, over the settled surface so that a depression does not report a
# negative height.

"""
    height_above_nearest_drainage(dem; threshold=100, streams=nothing, closed=nothing)

Height above the nearest downstream drainage cell (Nobre et al., 2011).

Streams are cells whose upstream area reaches `threshold` (in the grid's area
units), or the cells of an explicit boolean `streams` mask. With a mask the raw
surface is used, matching Geomorphometry's two-argument form; with a threshold
the depression-filled surface is used, matching its keyword form.

`D8` only: a multi-direction method gives a cell several downstream neighbors
and therefore several candidate heights, which is a different algorithm.
"""
function height_above_nearest_drainage(ras, grid::AbstractGridSpec; threshold=100,
        streams=nothing, closed=nothing, table=neighbortable(grid, NeighborRings(1)),
        sweep=floodsweep(ras, grid; closed, table))
    n = length(_data(ras))
    surface, mask = _handsurface(streams, threshold, ras, grid, sweep, n)
    out = Vector{Float64}(undef, n)
    _handpass!(out, sweep.order, sweep.down, surface, mask)
    return reshape(out, size(_data(ras)))
end

function _handsurface(::Nothing, threshold, ras, grid, sweep, n)
    acc = Vector{Float64}(undef, n)
    _fillcellareas!(acc, grid)
    _accumulatedown!(acc, sweep.order, sweep.down)
    mask = acc .>= threshold
    _burnstreams!(mask, sweep.order, sweep.down)
    return _settledvalues(ras, sweep), mask
end

function _handsurface(streams, threshold, ras, grid, sweep, n)
    mask = collect(Bool, vec(_data(streams)))
    length(mask) == n ||
        throw(ArgumentError("`streams` has $(length(mask)) entries for $n cells"))
    return _floatvalues(vec(_data(ras))), mask
end

# Reverse pop order is a topological order, so one pass carries every stream
# cell's mark all the way down its own chain.
function _burnstreams!(mask, order, down)
    @inbounds for j in reverse(eachindex(order))
        p = order[j]
        mask[p] || continue
        d = down[p]
        d == 0 && continue
        mask[d] = true
    end
    return mask
end

function _handpass!(out, order, down, surface, mask)
    fill!(out, NaN)
    @inbounds for p in order
        d = down[p]
        if mask[p] || d == 0
            out[p] = 0.0
        elseif isfinite(surface[p]) && isfinite(surface[d]) && isfinite(out[d])
            out[p] = out[d] + surface[p] - surface[d]
        else
            out[p] = NaN
        end
    end
    @inbounds for p in eachindex(out)
        out[p] = max(0.0, out[p])
    end
    return out
end

# ### Multi-direction flow
#
# `down` carries one downstream position per cell, which is the whole of D8 and
# none of DInf or FD8: those split a cell's flow between several neighbors. The
# sweep result therefore grows a second, *ragged* structure alongside `down`
# rather than replacing it — `FloodSweep` stayed opaque so this is an addition,
# not a break.
#
# The exclusion rule is the one Geomorphometry expresses with a `visited` array
# swept in reverse pop order: a cell may only send flow to a neighbor that
# popped *before* it. Stated as a rank comparison it is position-keyed, order
# independent and cheap, which is what lets the partition be built in one pass
# over positions instead of interleaved with the accumulation.

"""
Ragged downstream partition in CSR form: cell `p`'s outgoing edges are
`starts[p]:starts[p+1]-1`, each an entry of `targets` and a matching `weights`.
A pit or an outlet has no edges; every other cell's weights sum to 1.
"""
struct FlowPartition{P<:Integer}
    starts::Vector{Int32}
    targets::Vector{P}
    weights::Vector{Float64}
end

Base.length(part::FlowPartition) = length(part.starts) - 1
@inline partitionrange(part::FlowPartition, p::Int) =
    (@inbounds Int(part.starts[p])):(@inbounds Int(part.starts[p + 1]) - 1)

"""
    flowpartition(method, ras, grid; closed, table, sweep) -> FlowPartition

The downstream flow partition a multi-direction `method` implies on the surface
`sweep` settled. `D8()` gives the degenerate one-edge partition.
"""
function flowpartition(method::FlowDirectionMethod, ras, grid::AbstractGridSpec;
        closed=nothing, table=neighbortable(grid, NeighborRings(1)),
        sweep=floodsweep(ras, grid; closed, table))
    n = length(_data(ras))
    rank = zeros(Int32, n)
    @inbounds for (j, p) in enumerate(sweep.order)
        rank[Int(p)] = Int32(j)
    end
    starts = Vector{Int32}(undef, n + 1)
    targets = similar(sweep.down, 0)
    weights = Float64[]
    _partitionpass!(method, starts, targets, weights, rank, sweep.down,
        vec(_data(ras)), _methodaspects(method, ras, grid), table, slotgeometry(grid))
    return FlowPartition(starts, targets, weights)
end

# DInf routes by the surface aspect, so it takes whichever aspect estimator the
# grid's `defaultmethod` names — Horn on a rectilinear grid, the tangent-plane
# fit on a cell grid. That is exactly Geomorphometry's own split, but here it is
# the estimator seam (decisions §7) rather than two hard-coded branches.
_methodaspects(::DInf, ras, grid) = vec(_data(aspect(ras, grid)))
_methodaspects(::FlowDirectionMethod, ras, grid) = nothing

# THE FUNCTION BARRIER, for the same reason `_floodsweep!` has one: `table` and
# `sg` are runtime-determined types.
function _partitionpass!(method, starts, targets::Vector{P}, weights, rank, down, z,
        aspects, tbl, sg) where {P}
    n = length(rank)
    @inbounds for p in 1:n
        starts[p] = Int32(length(targets) + 1)
        d = Int(down[p])
        d == 0 && continue
        _partitioncell!(method, targets, weights, rank, down, p, d, z, aspects, tbl, sg)
    end
    starts[n + 1] = Int32(length(targets) + 1)
    return nothing
end

# A neighbor may receive flow when it popped earlier — `rank` is the pop index —
# or when it is an outlet, which emits nothing and therefore cannot close a
# cycle whatever its rank. Rank 0 means "never visited", which covers closed and
# unreachable cells at once, and those are never valid targets.
#
# Geomorphometry spells the same rule as a `visited` array swept in reverse pop
# order, but it never marks an outlet visited (it `continue`s before the mark),
# so an outlet is permanently "unvisited" there. That accident is the second
# clause here, made deliberate: without it a border outlet that happens to pop
# late is refused, and the two implementations disagree on ~1% of cells.
@inline function _isdownstream(rank, down, q::Int, p::Int)
    q == 0 && return false
    @inbounds rq = rank[q]
    rq == 0 && return false
    return rq < (@inbounds rank[p]) || (@inbounds down[q]) == 0
end

@inline function _singleedge!(targets::Vector{P}, weights, d) where {P}
    push!(targets, P(d))
    push!(weights, 1.0)
    return nothing
end

function _partitioncell!(::D8, targets, weights, rank, down, p, d, z, aspects, tbl, sg)
    return _singleedge!(targets, weights, d)
end

# DInf: bracket the aspect with the two angularly nearest neighbors and split
# the flow between them in proportion to the angular gaps (Tarboton, 1997).
function _partitioncell!(::DInf, targets::Vector{P}, weights, rank, down, p, d, z,
        aspects, tbl, sg) where {P}
    asp = @inbounds aspects[p]
    isfinite(asp) || return _singleedge!(targets, weights, d)
    origin = cellorigin(sg, p)
    lowgap = upgap = 360.0
    lowq = upq = 0
    for (k, q) in slots(tbl, p)
        q == 0 && continue
        b = edgeat(sg, origin, p, k, q).bearing
        g1 = mod(asp - b, 360.0)
        g2 = mod(b - asp, 360.0)
        g1 < lowgap && ((lowgap, lowq) = (g1, q))
        g2 < upgap && ((upgap, upq) = (g2, q))
    end
    # A cell whose bracketing pair does not contain its own downstream neighbor
    # is inside a depression on the raw surface: the aspect points somewhere the
    # flood does not go, so the flood wins.
    (lowq == 0 || (lowq != d && upq != d)) && return _singleedge!(targets, weights, d)
    wa, wb = if lowq == upq
        (1.0, 0.0)
    else
        gap = lowgap + upgap
        (upgap / gap, lowgap / gap)
    end
    # At most one of the two can be upstream, because one of them is `d`, and
    # `d` popped before `p` by construction.
    _isdownstream(rank, down, lowq, p) || ((wa, wb) = (0.0, 1.0))
    _isdownstream(rank, down, upq, p) || ((wa, wb) = (1.0, 0.0))
    wa > 0 && (push!(targets, P(lowq)); push!(weights, wa))
    wb > 0 && (push!(targets, P(upq)); push!(weights, wb))
    return nothing
end

# FD8: every downslope neighbor receives flow, weighted by
# `(Δz/distance · contour_length)^p` (Quinn et al., 1991). The contour length is
# the angular one — half the distance times the tangents of the half-gaps to the
# angular neighbors — which is defined on any ring, not a table of eight
# constants.
function _partitioncell!(fd8::FD8, targets::Vector{P}, weights, rank, down, p, d, z,
        aspects, tbl, sg) where {P}
    origin = cellorigin(sg, p)
    zp = @inbounds z[p]
    total = 0.0
    dweight = 0.0
    for (k, q) in slots(tbl, p)
        w = _fd8weight(fd8, zp, rank, down, p, k, q, z, tbl, sg, origin)
        total += w
        q == d && (dweight = w)
    end
    (total <= 0 || dweight <= 0) && return _singleedge!(targets, weights, d)
    for (k, q) in slots(tbl, p)
        w = _fd8weight(fd8, zp, rank, down, p, k, q, z, tbl, sg, origin)
        w > 0 || continue
        push!(targets, P(q))
        push!(weights, w / total)
    end
    return nothing
end

function _fd8weight(fd8, zp, rank, down, p::Int, k::Int, q::Int, z, tbl, sg, origin)
    _isdownstream(rank, down, q, p) || return 0.0
    difference = zp - (@inbounds z[q])
    (difference < 0 || isnan(difference)) && return 0.0
    edge = edgeat(sg, origin, p, k, q)
    edge.distance > 0 || return 0.0
    previous = next = Inf
    for (k2, q2) in slots(tbl, p)
        (q2 == 0 || q2 == q) && continue
        other = edgeat(sg, origin, p, k2, q2).bearing
        back = mod(edge.bearing - other, 360.0)
        forward = mod(other - edge.bearing, 360.0)
        iszero(back) || (previous = min(previous, back))
        iszero(forward) || (next = min(next, forward))
    end
    (isfinite(previous) && isfinite(next)) || return 0.0
    contour = edge.distance / 2 * (tand(previous / 2) + tand(next / 2))
    contour > 0 || return 0.0
    return (difference / edge.distance * contour)^fd8.p
end

# The same reverse-order pass `_accumulatedown!` makes, with the single
# downstream position replaced by the cell's edge list.
function _accumulateweighted!(acc, order, part::FlowPartition)
    @inbounds for j in reverse(eachindex(order))
        p = Int(order[j])
        a = acc[p]
        for e in partitionrange(part, p)
            acc[Int(part.targets[e])] += part.weights[e] * a
        end
    end
    return acc
end

function _flowaccumulation(method::Union{DInf,FD8}, ras, grid::AbstractGridSpec, closed)
    table = neighbortable(grid, NeighborRings(1))
    sweep = floodsweep(ras, grid; closed, table)
    part = flowpartition(method, ras, grid; table, sweep)
    acc = Vector{Float64}(undef, length(_data(ras)))
    _fillcellareas!(acc, grid)
    _accumulateweighted!(acc, sweep.order, part)
    dirs = _multidirections(ras, grid, table, part)
    return reshape(acc, size(_data(ras))), dirs
end
