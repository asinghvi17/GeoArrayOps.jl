# ### (e) Local statistics against hand-computed values
#
# `bump` is flat except for one raised cell, so every neighborhood statistic has
# an obvious closed form.
bump = zeros(5, 5)
bump[3, 3] = 10.0
bump_ras = Raster(bump, (X(1.0:5.0), Y(1.0:5.0)))
tpi_bump = topographic_position_index(bump_ras)
tri_bump = terrain_ruggedness_index(bump_ras)
rough_bump = roughness(bump_ras)
@test Rasters.name(tpi_bump) == :topographic_position_index
# The peak sits 10 above the mean of its eight zero neighbors.
@test parent(tpi_bump)[3, 3] ≈ 10.0
# A cell orthogonally adjacent to the peak has 8 neighbors, one of them the peak.
@test parent(tpi_bump)[3, 2] ≈ -10.0 / 8
# A corner cell has only 3 neighbors, none of them the peak.
@test parent(tpi_bump)[1, 1] ≈ 0.0
@test parent(tri_bump)[3, 3] ≈ sqrt(8 * 10.0^2)
@test parent(tri_bump)[3, 2] ≈ 10.0
@test parent(tri_bump)[1, 1] ≈ 0.0
@test parent(terrain_ruggedness_index(bump_ras; normalize=true))[3, 3] ≈ sqrt(8 * 100 / 8)
@test parent(terrain_ruggedness_index(bump_ras; normalize=true, squared=false))[3, 3] ≈ 10.0
@test parent(rough_bump)[3, 3] ≈ 10.0
@test parent(rough_bump)[3, 2] ≈ 10.0
@test parent(rough_bump)[1, 1] ≈ 0.0
# Storage order does not change a local statistic.
bump_yx = Raster(permutedims(bump), (Y(1.0:5.0), X(1.0:5.0)))
@test parent(topographic_position_index(bump_yx)) ≈ permutedims(parent(tpi_bump))
@test parent(roughness(bump_yx)) ≈ permutedims(parent(rough_bump))
