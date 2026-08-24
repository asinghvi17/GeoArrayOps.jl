# ### (f) The sweep family
#
# Geomorphometry is the reference. Where the answer is uniquely determined — the
# settled surface always, the flood tree whenever no two cells tie — the
# assertion is `==`; where it is not, the assertion is the invariant, never a
# loosened tolerance.

# The PoC collected the deliberate no-outlet warnings with a hand-rolled logger
# so they were asserted rather than leaked into the output; here `@test_logs`
# does both jobs, and `silently` (in `runtests.jl`) covers the one place a
# warning only needs suppressing.

# One flood per fixture, so `down`, `acc` and the directions all describe the
# same tree — and so the tests exercise the same `table`-threading a composite
# operation uses.
function sweepparts(ras, grid; closed=nothing)
    table = neighbortable(grid)
    sweep = floodsweep(ras, grid; closed, table)
    acc = Vector{Float64}(undef, length(_data(ras)))
    _fillcellareas!(acc, grid)
    _accumulatedown!(acc, sweep.order, sweep.down)
    dirs = flowdirection(ras, grid; sweep, table)
    return sweep, reshape(acc, size(_data(ras))), dirs
end

# Walk `down` to its root, refusing to loop forever.
function downroot(down, p)
    steps = 0
    while down[p] != 0
        p = Int(down[p])
        (steps += 1) <= length(down) || error("`down` contains a cycle")
    end
    return p
end

# The linear `getindex` the traversal work arrays use must agree with the
# Cartesian one the record API uses.
@test all(I -> ca[LinearIndices(size(ca))[I]] == ca[I[1], I[2]],
    CartesianIndices(size(ca)))
@test all(I -> ca_geo[LinearIndices(size(ca_geo))[I]] == ca_geo[I[1], I[2]],
    CartesianIndices(size(ca_geo)))
@test Base.IndexStyle(typeof(ca)) isa Base.IndexLinear
@test Base.IndexStyle(typeof(ca_geo)) isa Base.IndexLinear
# A cell grid's areas are lazy: no residency, whatever the grid size.
@test cellarea(cell_grid) isa CellAreas
@test sizeof(cellarea(cell_grid)) < 8 * length(cell_raster)
