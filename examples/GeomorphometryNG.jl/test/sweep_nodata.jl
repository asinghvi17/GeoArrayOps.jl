# #### (4, 6) Nodata
#
# By default NaN participates, exactly as Geomorphometry does: NaN sorts last in
# the queue, so those cells are visited at the end, get a parent, and inject their
# own area into the accumulation. `closed = isnan.(z)` is the one-liner that
# excludes them instead.
nan_z = Float64[20.0 - i - 0.3j for i in 1:8, j in 1:8]
nan_z[3:5, 3:5] .= NaN # A 3x3 hole, so its center's whole neighborhood is nodata
_, nan_grid = spatialparts(nan_z; spacing=(2.0, -5.0))
nan_table = neighbortable(nan_grid)
nan_holes = findall(isnan, vec(nan_z))
@test all(((k, q),) -> q == 0 || isnan(nan_z[q]), slots(nan_table, LinearIndices(nan_z)[4, 4]))

nan_sweep, nan_acc, nan_dirs = sweepparts(nan_z, nan_grid)
nan_settled = settle(nan_z, nan_grid; sweep=nan_sweep)
@test length(nan_sweep.order) == length(nan_z)
@test all(p -> nan_sweep.down[p] != 0, nan_holes) # Visited, and given a parent
@test all(p -> isnan(nan_settled[p]), nan_holes)  # But with no settled elevation
@test all(p -> isfinite(nan_settled[p]), setdiff(eachindex(vec(nan_z)), nan_holes))
@test sum(vec(nan_acc)[nan_sweep.down .== 0]) ≈ sum(cellarea(nan_grid))

closed_sweep, closed_acc, closed_dirs = sweepparts(nan_z, nan_grid; closed=isnan.(nan_z))
@test length(closed_sweep.order) == length(nan_z) - length(nan_holes)
@test all(p -> closed_sweep.down[p] == 0, nan_holes)
@test all(p -> closed_acc[p] == 10.0, nan_holes) # Closed cells emit and receive nothing
@test all(p -> ispit(closed_dirs[p]), nan_holes)
let reached = (vec(closed_sweep.down) .== 0) .& .!isnan.(vec(nan_z))
    @test sum(vec(closed_acc)[reached]) ≈ sum(cellarea(nan_grid)) - length(nan_holes) * 10.0
end

# (6) A raster with nothing to seed says so instead of looping or silently
# returning `acc == cellarea`.
allnan_z = fill(NaN, 4, 4)
_, allnan_grid = spatialparts(allnan_z)
allnan_sweep = @test_logs (:warn, r"every cell is closed") match_mode=:any begin
    floodsweep(allnan_z, allnan_grid; closed=isnan.(allnan_z))
end
@test isempty(allnan_sweep.order) && all(==(0), allnan_sweep.down)

# #### (11) `closed` disconnecting a component
#
# A ring of closed cells encloses the bowl's four inner cells, so fewer cells are
# reachable than `order` was sized for. The truncation is load-bearing: without
# it the untruncated tail of ones would inflate one arbitrary cell.
_, bowl_grid = spatialparts(bowl_z; spacing=(2.0, -5.0))
bowl_ring = falses(6, 6)
bowl_ring[2:5, 2] .= true
bowl_ring[2:5, 5] .= true
bowl_ring[2, 2:5] .= true
bowl_ring[5, 2:5] .= true
cut_sweep, cut_acc, cut_dirs = sweepparts(bowl_z, bowl_grid; closed=bowl_ring)
cut_inner = vec(LinearIndices((6, 6))[3:4, 3:4])
@test length(cut_sweep.order) == length(bowl_z) - count(bowl_ring) - length(cut_inner)
@test length(cut_sweep.order) < length(bowl_z) - count(bowl_ring)
@test all(p -> cut_sweep.down[p] == 0, cut_inner)
@test all(p -> cut_acc[p] == 10.0, cut_inner) # Never touched by the `ones` tail
let reached = (vec(cut_sweep.down) .== 0) .& .!vec(bowl_ring)
    reached[cut_inner] .= false
    @test sum(vec(cut_acc)[reached]) ≈ length(cut_sweep.order) * 10.0
end
