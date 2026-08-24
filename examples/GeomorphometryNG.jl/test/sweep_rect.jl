# #### Fixture 1 — a bowl with a flat and a pit
#
# Two rectilinear DEMs, each run in `(X, Y)` and `(Y, X)` storage with unequal,
# sign-flipped spacing. `pit_z` has a uniform rim, so the whole interior lies
# below it and the fill raises everything; `bowl_z` has a graded rim, so only the
# enclosed cells rise — the two-cell flat and the single-cell pit spill over the
# ring of 9s.
pit_z = Float64[
    10 10 10 10 10 10
    10 8 7 7 6 10
    10 7 3 3 5 10
    10 6 3 1 4 10
    10 5 4 4 4 10
    10 10 10 10 10 10]
bowl_z = Float64[
    1 2 3 4 5 6
    2 9 9 9 9 7
    3 9 4 4 9 8
    4 9 4 1 9 9
    5 9 9 9 9 10
    6 7 8 9 10 11]
sweep_xs, sweep_ys = 0.0:2.0:10.0, 25.0:-5.0:0.0
asxy(z) = Raster(z, (X(sweep_xs), Y(sweep_ys)))
asyx(z) = Raster(permutedims(z), (Y(sweep_ys), X(sweep_xs)))

for z in (pit_z, bowl_z)
    # (2) Nested depressions: `settle` is the minimax fill, so it equals
    # `filldepressions` exactly. On Float64 there is no tolerance to hide in.
    @test parent(settle(asxy(z))) == GM.filldepressions(z)
    # (8) Both storage orders. Settled elevation is a geographic quantity, so the
    # two layouts are transposes of each other.
    @test parent(settle(asxy(z))) == permutedims(parent(settle(asyx(z))))
    @test all(parent(settle(asxy(z))) .>= z)
    # (7) Integer elevations take the same path, promoted to float on the way in.
    @test parent(settle(asxy(Int.(z)))) == GM.filldepressions(Int.(z))
    @test parent(settle(asxy(UInt8.(z)))) == GM.filldepressions(UInt8.(z))
end
@test parent(settle(asxy(pit_z)))[4, 4] == 10.0 # A uniform rim fills the whole bowl
@test parent(settle(asxy(bowl_z)))[4, 4] == 9.0 # The pit spills over the 9-ring

function checkrectsweep(z)
    _, g = spatialparts(asxy(z))
    sweep, area, direction = sweepparts(asxy(z), g)
    @test length(sweep.order) == length(z) # Every cell reachable, every slot used
    # Conservation: all area reaches a cell with no downstream neighbor.
    @test sum(vec(area)[sweep.down .== 0]) ≈ sum(cellarea(g))
    @test all(area .>= 10.0)
    # (3) No cycles: every `down` chain terminates at an outlet.
    @test all(p -> sweep.down[downroot(sweep.down, p)] == 0, eachindex(vec(z)))
    # Pits are exactly the seeds; no interior pit survives the fill.
    @test count(ispit, direction) == sweep.nseeds ==
            count(p -> isboundary(neighbortable(g), p), eachindex(vec(z)))
    return nothing
end
checkrectsweep(pit_z)
checkrectsweep(bowl_z)

# #### Fixture 2 — a tie-free surface, exact parity with Geomorphometry
#
# Ties make the flood tree non-unique, and the two priority queues break them
# differently (Geomorphometry's is `Dict`-backed; this one is array-backed). With
# every value distinct the tree is unique, so the direction rasters must agree
# cell for cell. `splitmix01` stands in for a seeded RNG, which is not a
# dependency here.
function splitmix01(k::Integer)
    h = UInt64(k) * 0x9E3779B97F4A7C15
    h ⊻= h >> 29
    h *= 0xBF58476D1CE4E5B9
    h ⊻= h >> 32
    h *= 0x94D049BB133111EB
    h ⊻= h >> 31
    return Float64(h >> 11) * (1 / 9007199254740992.0)
end
tf_n = 200
tf_z = [100.0 * splitmix01(i + tf_n * (j - 1)) + 1e-9 * (i + tf_n * (j - 1))
        for i in 1:tf_n, j in 1:tf_n]
@test length(unique(tf_z)) == length(tf_z) # The premise of every `==` below
tf_ras = Raster(tf_z, (X(range(0.0; step=2.0, length=tf_n)),
    Y(range(1000.0; step=-5.0, length=tf_n))))
_, tf_grid = spatialparts(tf_ras)
tf_sweep, tf_acc, tf_dirs = sweepparts(tf_ras, tf_grid)
tf_settled = settle(tf_ras, tf_grid; sweep=tf_sweep)
gm_tf_acc, gm_tf_dirs = GM.flowaccumulation(tf_z; method=D8(), cellsize=(2.0, -5.0))

@test tf_settled == GM.filldepressions(tf_z)
# The flood tree itself, compared through the LDD codec. No `_orient`, no
# `cellsize` sign logic — `storageoffset` already absorbed storage order.
@test Int.(tf_dirs) == Int.(gm_tf_dirs)
# Geomorphometry accumulates in Float32 over the identical addition sequence, so
# the difference is its rounding, not a different answer.
@test maximum(abs.(tf_acc .- Float64.(gm_tf_acc))) / maximum(tf_acc) < 1e-5
@test all(isapprox.(Float32.(tf_acc), gm_tf_acc; rtol=1e-5))
@test sum(vec(tf_acc)[tf_sweep.down .== 0]) ≈ sum(cellarea(tf_grid))
# Only the algorithm's own output and work arrays: `order` and `down` at 4 bytes
# each, `closed` at a bit, `acc` at 8, the directions at 1, and the queue. The
# rectilinear neighbor table is arithmetic and holds no memory at all.
@test sizeof(neighbortable(tf_grid)) < 1024
flowaccumulation(tf_ras, tf_grid) # Warm up, then measure the steady state
@test (@allocated flowaccumulation(tf_ras, tf_grid)) < 128 * length(tf_z)

# (8) Storage order, on the fixture where accumulation is determined: identical
# LDD codes, not permuted ones, because the codes are geographic.
tf_yx = Raster(permutedims(tf_z), (Y(range(1000.0; step=-5.0, length=tf_n)),
    X(range(0.0; step=2.0, length=tf_n))))
_, tf_grid_yx = spatialparts(tf_yx)
tf_sweep_yx, tf_acc_yx, tf_dirs_yx = sweepparts(tf_yx, tf_grid_yx)
@test tf_acc == permutedims(tf_acc_yx)
@test Int.(tf_dirs) == permutedims(Int.(tf_dirs_yx))
@test tf_settled == permutedims(settle(tf_yx, tf_grid_yx; sweep=tf_sweep_yx))

# (1) A single-cell pit: filled to its spill elevation, and routed rather than
# left a sink. Its direction points uphill on the raw DEM, which is what correct
# depression routing looks like.
tf_table = neighbortable(tf_grid)
tf_pits = findall(eachindex(tf_z)) do p
    !isboundary(tf_table, p) && all(((k, q),) -> tf_z[q] > tf_z[p], slots(tf_table, p))
end
@test !isempty(tf_pits)
let p = first(tf_pits)
    @test tf_settled[p] > tf_z[p]
    @test !ispit(tf_dirs[p])
    @test tf_settled[p] == max(tf_z[p], tf_settled[tf_sweep.down[p]])
end

# #### (3) A flat draining to two outlets
#
# A plateau at 5 inside a rim at 9, with two low notches in the rim. Which flat
# cell goes to which notch is heap-order dependent, so the assertions are the
# partition and the conservation, not the assignment.
flat_z = fill(5.0, 5, 7)
flat_z[1, :] .= 9.0
flat_z[5, :] .= 9.0
flat_z[:, 1] .= 9.0
flat_z[:, 7] .= 9.0
flat_z[1, 3] = 1.0
flat_z[5, 5] = 2.0
_, flat_grid = spatialparts(flat_z; spacing=(2.0, -5.0))
flat_sweep, flat_acc, flat_dirs = sweepparts(flat_z, flat_grid)
flat_lin = LinearIndices(flat_z)
flat_interior = vec(flat_lin[2:4, 2:6])
@test all(p -> flat_sweep.down[downroot(flat_sweep.down, p)] == 0, eachindex(vec(flat_z)))
flat_roots = unique(downroot(flat_sweep.down, p) for p in flat_interior)
@test Set(flat_roots) == Set((flat_lin[1, 3], flat_lin[5, 5])) # The two notches, both used
@test flat_acc[1, 3] + flat_acc[5, 5] ≈ (length(flat_interior) + 2) * 10.0
@test sum(vec(flat_acc)[flat_sweep.down .== 0]) ≈ sum(cellarea(flat_grid))
@test all(parent(settle(flat_z, flat_grid; sweep=flat_sweep))[flat_interior] .== 5.0)
