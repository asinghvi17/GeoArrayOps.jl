# ### The rest of the hydrology
#
# `tf_ras` (200×200, every value distinct) makes the flood tree unique, so
# Geomorphometry and this package settle the same surface and the comparisons
# below are about the algorithms rather than about heap order.

# #### Depressions
@test parent(depression_depth(tf_ras)) == GM.filldepressions(tf_z) .- tf_z
@test all(parent(depression_depth(tf_ras)) .>= 0.0)
# `depression_volume` integrates against `cellarea`, so it carries units — and
# it is a number, which the façade must hand back unwrapped. Geomorphometry's
# matrix path assumes unit cells, so that is the comparison; `tf_ras`'s own
# 2 × 5 cells make it ten times larger.
@test depression_volume(tf_z) ≈ GM.depression_volume(tf_z)
@test depression_volume(tf_z) isa Real
@test depression_volume(tf_ras) ≈ 10 * depression_volume(tf_z)
# A bowl with a known rim has a known volume.
bowl_depth = parent(depression_depth(asxy(pit_z)))
@test sum(bowl_depth) ≈ sum(GM.filldepressions(pit_z) .- pit_z)
@test depression_volume(asxy(pit_z)) ≈ sum(bowl_depth) * 10.0

# #### Height above nearest drainage
#
# One flood, one adjacency table, `down` positions throughout. Geomorphometry
# reuses its `Float32` accumulation buffer as the filled surface, so the heights
# it reports carry `Float32` precision; that is the whole of the tolerance here.
hand_ng = height_above_nearest_drainage(tf_ras; threshold=1000)
hand_gm = GM.height_above_nearest_drainage(tf_z; method=D8(), threshold=1000,
    cellsize=(2.0, -5.0))
@test all(isapprox.(parent(hand_ng), hand_gm; atol=1e-3))
@test all(>=(0.0), parent(hand_ng))
@test Rasters.name(hand_ng) == :height_above_nearest_drainage
# With an explicit stream mask the surface is the raw one, which is
# Geomorphometry's two-argument form exactly — Float64 on both sides, so exact.
hand_mask = falses(size(tf_z))
hand_mask[100, :] .= true
@test parent(height_above_nearest_drainage(tf_ras; streams=hand_mask)) ==
      GM.height_above_nearest_drainage(tf_z, hand_mask)
# Every stream cell is at height zero, and so is every outlet.
hand_streamed = parent(height_above_nearest_drainage(tf_ras; streams=hand_mask))
@test all(hand_streamed[hand_mask] .== 0.0)
@test all(hand_streamed[vec(tf_sweep.down) .== 0] .== 0.0)
# A lower threshold marks more streams, so no cell can be further from one.
hand_low = parent(height_above_nearest_drainage(tf_ras; threshold=200))
@test all(hand_low .<= parent(hand_ng) .+ 1e-9)
# A mask the wrong size is refused rather than silently clipped.
@test (try height_above_nearest_drainage(tf_ras; streams=falses(3, 3))
    nothing catch e; e end) isa ArgumentError

# #### Accumulation-derived indices
hyd_int(A) = view(Array(A), 2:(tf_n - 1), 2:(tf_n - 1))
hyd_twi = parent(topographic_wetness_index(tf_ras))
hyd_spi = parent(stream_power_index(tf_ras))
hyd_dp = parent(drainage_potential(tf_ras))
@test all(isapprox.(hyd_int(hyd_twi),
    hyd_int(GM.topographic_wetness_index(tf_z; method=D8(), cellsize=(2.0, -5.0)));
    atol=1e-4))
@test all(isapprox.(hyd_int(hyd_spi),
    hyd_int(GM.stream_power_index(tf_z; method=D8(), cellsize=(2.0, -5.0))); atol=1e-4))
@test all(isapprox.(hyd_int(hyd_dp),
    hyd_int(GM.drainage_potential(tf_z; method=D8(), cellsize=(2.0, -5.0))); atol=1e-4))
# TWI and SPI differ only in the sign of the slope term, so where both are
# finite they sum to twice the log of the accumulation.
hyd_pair = hyd_twi .+ hyd_spi .- 2 .* log.(tf_acc)
@test all(<(1e-9), abs.(filter(isfinite, hyd_pair)))
@test count(isfinite, hyd_pair) > length(tf_z) - 4 * tf_n
@test Rasters.name(topographic_wetness_index(tf_ras)) == :topographic_wetness_index
@test TWI === topographic_wetness_index && SPI === stream_power_index

# #### Cell grids
#
# A complete sphere has no outlet, so the no-outlet warning is the contract and
# is asserted rather than silenced.
cell_hand = @test_logs (:warn, r"no domain boundary") match_mode = :any begin
    height_above_nearest_drainage(cell_raster; threshold=1e14)
end
@test all(>=(0.0), parent(cell_hand))
@test all(isfinite, parent(cell_hand))
@test count(==(0.0), parent(cell_hand)) >= 1 # At least the seeded minimum
# The subtree has a rim, so it floods to real outlets and needs no warning.
sub_hand = height_above_nearest_drainage(sub_raster; threshold=1e11)
@test all(>=(0.0), parent(sub_hand))
@test all(parent(sub_hand)[collect(DGG.border(sub_cells))] .== 0.0)
# HAND is measured on the settled surface, so it is zero exactly where the
# stream mask is — the invariant that survives on a grid with no baseline.
sub_stream = falses(length(sub_cells))
sub_stream[1:10] .= true
sub_hand_masked = height_above_nearest_drainage(sub_raster; streams=sub_stream)
@test all(parent(sub_hand_masked)[1:10] .== 0.0)

silently() do
    @test all(parent(depression_depth(cell_raster)) .>= -1e-9)
    @test depression_volume(cell_raster) >= 0.0
    @test all(isfinite, parent(topographic_wetness_index(cell_raster)))
end
# A region of a cell grid settles into real depressions, and the volume is the
# depth integrated against the true spherical cell areas.
sub_depth = parent(depression_depth(sub_raster))
@test all(sub_depth .>= -1e-9)
@test depression_volume(sub_raster) ≈ sum(sub_depth .* cellarea(sub_grid))
