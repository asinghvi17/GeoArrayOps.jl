# ### The rest of the neighborhood statistics
#
# Reuses `td_ras` (north-up, unequal spacing) from the terrain section.

# #### Counting statistics — exact, both backends
@test tdin(parent(prominence(td_ras))) == tdin(GM.prominence(td_z))
@test tdin(parent(percentile_elevation(td_ras))) == tdin(GM.percentile_elevation(td_z))
@test tdin(parent(pitremoval(td_ras))) == tdin(GM.pitremoval(td_z))
@test eltype(parent(prominence(td_ras))) == Int8
# A pit is raised to its lowest neighbor; everything else keeps its value.
pit_patch = fill(5.0, 5, 5)
pit_patch[3, 3] = 1.0
pit_ras = Raster(pit_patch, (X(1.0:5.0), Y(1.0:5.0)))
@test parent(pitremoval(pit_ras))[3, 3] == 5.0
@test parent(pitremoval(pit_ras; limit=10.0))[3, 3] == 1.0 # Too shallow to count
# Geomorphometry returns `typemax(eltype)` for a cell whose neighbors are all
# exactly equal to it, because its accumulator starts there and a tie never
# clears the pit flag. Here such a cell keeps its own elevation.
flat_patch = Raster(zeros(5, 5), (X(1.0:5.0), Y(1.0:5.0)))
@test parent(pitremoval(flat_patch))[3, 3] == 0.0
@test isinf(GM.pitremoval(zeros(5, 5))[3, 3])

# #### Residual roughness
#
# Two record passes, so the second pass reads the first pass's frame — the
# comparison is two cells in rather than one.
@test all(isapprox.(tdin(parent(roughness_index_elevation(td_ras)), 2),
    tdin(GM.roughness_index_elevation(td_z), 2); atol=1e-9))
@test all(>=(0.0), parent(roughness_index_elevation(td_ras)))
@test all(parent(roughness_index_elevation(flat_patch)) .== 0.0)

# #### Rugosity
#
# The fan is built from bearings, so it generalizes off the square grid — and on
# a square grid it is Geomorphometry's number. Geomorphometry accumulates the
# eight elevation differences in a `Float32` buffer, hence the tolerance.
td_square = Raster(td_z, (X(range(0.0; step=2.0, length=td_n)),
    Y(range(0.0; step=2.0, length=td_n))))
@test all(isapprox.(tdin(parent(rugosity(td_square))),
    tdin(GM.rugosity(td_z; cellsize=(2.0, 2.0))); rtol=1e-6))
# Flat ground is exactly 1 whatever the cell shape, and rugosity never dips
# below 1: a surface cannot be smaller than its own projection.
@test all(parent(rugosity(Raster(fill(4.0, 6, 6), (X(1.0:6.0), Y(1.0:6.0))))) .== 1.0)
@test all(parent(rugosity(td_ras)) .>= 1.0 - 1e-12)
# With unequal cell sizes the two disagree, because Geomorphometry gives its
# edge-neighbor vectors the *other* axis's length: `(δx, 0, dz)` for a step
# along dim2 and `(δy, 0, dz)` for a step along dim1.
@test !all(isapprox.(tdin(parent(rugosity(td_ras))),
    tdin(GM.rugosity(td_z; cellsize=(td_dx, td_dy))); rtol=1e-3))
# Storage order does not enter a ratio of areas.
@test isequal(reverse(parent(rugosity(td_flip)); dims=2), parent(rugosity(td_ras)))

# #### Entropy
#
# Geomorphometry cannot run this at all: its histogram buffers hold nine entries
# while its window is 5×5, so any neighborhood with more than nine distinct
# binned values indexes past the end.
gm_entropy_error = try GM.entropy(td_z); nothing catch e; e end
@test gm_entropy_error isa BoundsError
td_entropy = entropy(td_ras)
@test all(isfinite, parent(td_entropy))
@test all(>=(0.0), parent(td_entropy))
# A neighborhood of one repeated value has zero entropy; 25 distinct values
# reach the 25-sample maximum.
@test parent(entropy(Raster(fill(2.0, 9, 9), (X(1.0:9.0), Y(1.0:9.0)))))[5, 5] == 0.0
distinct_ras = Raster(Float64.(reshape(1:81, 9, 9)), (X(1.0:9.0), Y(1.0:9.0)))
@test parent(entropy(distinct_ras; step=nothing))[5, 5] ≈ log(25)
# Binning is what makes it a histogram rather than a cell count.
@test parent(entropy(distinct_ras; step=100.0))[5, 5] < 0.7
@test parent(entropy(distinct_ras; step=nothing))[5, 5] >
      parent(entropy(distinct_ras; step=4.0))[5, 5]
# One ring is a 3×3 window, two is Geomorphometry's 5×5.
@test parent(entropy(distinct_ras; step=nothing, rings=NeighborRings(1)))[5, 5] ≈ log(9)
# NaN is skipped, not binned, and a NaN center is NaN out.
entropy_nan = copy(td_z)
entropy_nan[5, 5] = NaN
entropy_nan_ras = Raster(entropy_nan, (X(1.0:td_n), Y(1.0:td_n)))
@test isnan(parent(entropy(entropy_nan_ras))[5, 5])
@test isfinite(parent(entropy(entropy_nan_ras))[6, 6])

# #### The one neighborhood the vocabulary cannot name
bpi_error = try bathymetric_position_index(td_ras); nothing catch e; e end
@test bpi_error isa ArgumentError
@test occursin("annulus", sprint(showerror, bpi_error))
@test occursin("cumulative", sprint(showerror, bpi_error))
@test (try bathymetric_position_index(cell_raster); nothing catch e; e end) isa ArgumentError

# #### Cell grids
#
# Everything above except BPI reads records only, so all of it runs on the
# sphere. The invariants are the ones that hold on any tessellation.
@test extrema(parent(prominence(cell_raster))) == (1, 6) # No pit, no peak on a tilt
@test all(x -> 0.0 <= x <= 1.0, parent(percentile_elevation(cell_raster)))
@test all(isfinite, parent(roughness_index_elevation(cell_raster)))
@test all(>=(0.0), parent(roughness_index_elevation(cell_raster)))
cell_rugosity = rugosity(cell_raster)
@test all(parent(cell_rugosity) .>= 1.0 - 1e-9)
# A hexagonal ring is a fan too: a smooth field over a sphere is very nearly
# flat at level 2, so its rugosity sits just above 1.
@test all(parent(cell_rugosity) .< 1.001)
@test all(parent(rugosity(Raster(fill(3.0, length(cell_raster)),
    (DGG.Cells(cell_lookup),)))) .≈ 1.0)
cell_entropy = entropy(cell_raster)
@test all(isfinite, parent(cell_entropy))
@test all(>=(0.0), parent(cell_entropy))
# Two rings of a hexagon is 18 neighbors, so 19 samples cap the entropy.
@test all(parent(cell_entropy) .<= log(19) + 1e-9)
@test all(parent(entropy(Raster(fill(3.0, length(cell_raster)),
    (DGG.Cells(cell_lookup),)))) .== 0.0)
# `pitremoval` on a tilted sphere touches only the handful of true pits.
@test count(!=(0.0), parent(pitremoval(cell_raster)) .- cell_values) <= 8
