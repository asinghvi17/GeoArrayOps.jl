# ### (d) Nodata semantics
#
# A NaN center yields NaN; a NaN neighbor is never selected; every other cell is
# untouched.
nan_data = copy(data_xy)
nan_data[3, 2] = NaN
nan_ras = Raster(nan_data, (X(xs), Y(ys)))
nan_slope = steepest_slope(nan_ras)
nan_dir = flow_direction(nan_ras)
nan_pf = slope(nan_ras; method=PlaneFit())
nan_horn = slope(nan_ras)
@test isnan(parent(nan_slope)[3, 2])
@test isnan(parent(nan_dir)[3, 2])
@test isnan(parent(nan_pf)[3, 2])
@test isnan(parent(nan_horn)[3, 2])
# A nodata neighbor is skipped rather than selected, and never poisons the
# result. The cell west of the hole would have flowed straight east into it; it
# now reports the steepest *available* descent, which is the southeast diagonal.
@test isfinite(parent(nan_slope)[2, 2])
@test parent(nan_slope)[2, 2] ≈ atand(6.0 / hypot(2.0, 5.0))
@test parent(nan_dir)[2, 2] ≈ mod(atand(2.0, -5.0), 360.0)
# The NaN cell's eastern neighbor would have flowed east anyway; the point is
# that it never picks the NaN cell, whose gradient comparison is skipped.
@test parent(nan_slope)[4, 2] ≈ parent(slope_xy)[4, 2]
@test parent(nan_dir)[4, 2] ≈ 90.0
# Cells that do not touch the hole are bit-identical to the NaN-free run.
for column in (1, 5)
    @test all(isequal.(parent(nan_slope)[column, :], parent(slope_xy)[column, :]))
    @test all(isequal.(parent(nan_dir)[column, :], parent(direction_xy)[column, :]))
    @test all(isequal.(parent(nan_pf)[column, :], parent(pf)[column, :]))
end
# Integer rasters have no NaN and keep working through the same kernels.
@test steepest_slope(Int.(data_xy); spacing=(2.0, -5.0)) ≈ matrix_slope
@test eltype(roughness(Int.(data_xy))) == Int
