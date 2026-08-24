# ### The windowed derivative family and its generic counterparts
#
# The fixture is north-up with X and Y both *ascending* and unequal spacing, so
# Geomorphometry's raw axis convention (dim1 east, dim2 north) coincides with
# this package's geographic one and a value-for-value comparison is meaningful.
# Storage-order independence is then checked separately, on the layout
# Geomorphometry cannot express.

td_n = 12
td_z = [100.0 + 5sin(x / 2.3) * cos(y / 1.7) + 0.7x - 0.4y + 0.01x * y
        for x in 1:td_n, y in 1:td_n]
td_dx, td_dy = 2.0, 3.0
td_ras = Raster(td_z, (X(range(0.0; step=td_dx, length=td_n)),
    Y(range(0.0; step=td_dy, length=td_n))))
td_cell = (td_dx, td_dy)
# Geomorphometry pads its window at the frame; this package answers NaN there.
# Every comparison below is therefore interior-only.
tdin(A, k=1) = view(Array(A), (1 + k):(td_n - k), (1 + k):(td_n - k))
# Geomorphometry computes the first-derivative family in Float32.
tdapprox(a, b; rtol=1e-5, k=1) = all(isapprox.(tdin(a, k), tdin(b, k); rtol, atol=1e-5))

# #### Slope estimators
@test tdapprox(parent(slope(td_ras)), GM.slope(td_z; cellsize=td_cell))
@test tdapprox(parent(slope(td_ras; method=ZevenbergenThorne())),
    GM.slope(td_z; cellsize=td_cell, method=GM.ZevenbergenThorne()))
# The frame is the one thing a windowed estimator cannot answer.
@test all(isnan, parent(slope(td_ras))[1, :])
@test all(isnan, parent(slope(td_ras))[:, end])
# MDG counts an uphill neighbor too, so it never reports less than the downhill
# half of the same quantity — and it reads records, so it runs on any grid.
@test all(parent(slope(td_ras; method=MDG())) .>= parent(steepest_slope(td_ras)) .- 1e-9)
# `exaggeration` scales the gradient, not the angle.
@test all(isapprox.(tand.(parent(slope(td_ras; exaggeration=2.5))),
    2.5 .* tand.(parent(slope(td_ras))); atol=1e-9, nans=true))
@test isequal(parent(pssm(td_ras)), parent(slope(td_ras; exaggeration=2.3)))
@test tdapprox(parent(pssm(td_ras)), GM.pssm(td_z; cellsize=td_cell))

# #### Aspect
#
# One convention for every estimator: the downslope bearing, clockwise from
# north. Geomorphometry's `aspect` agrees; its `hillshade` reads the same angle
# without the `compass` conversion that produces it, which is why the two are
# ported through one shared definition here.
@test tdapprox(parent(aspect(td_ras)), GM.aspect(td_z; cellsize=td_cell))
@test tdapprox(parent(aspect(td_ras; method=ZevenbergenThorne())),
    GM.aspect(td_z; cellsize=td_cell, method=GM.ZevenbergenThorne()))
@test Rasters.name(aspect(td_ras)) == :aspect
# A plane has one exact aspect, so every estimator must land on it — including
# the record-based one, which is the cell grids' default.
tilt_z = [500.0 - 3.0x - 1.0y for x in 0.0:2.0:20.0, y in 0.0:2.0:20.0]
tilt_ras = Raster(tilt_z, (X(0.0:2.0:20.0), Y(0.0:2.0:20.0)))
tilt_bearing = mod(atand(3.0, 1.0), 360.0) # Descends east and north
tiltin(A) = view(Array(A), 2:(size(tilt_z, 1) - 1), 2:(size(tilt_z, 2) - 1))
for m in (Horn(), ZevenbergenThorne(), PlaneFit())
    @test all(isapprox.(tiltin(parent(aspect(tilt_ras; method=m))), tilt_bearing; atol=1e-9))
end
# A flat surface has no aspect at all. Geomorphometry answers `0.0` — north.
@test all(isnan, parent(aspect(Raster(fill(7.0, 6, 6), (X(1.0:6.0), Y(1.0:6.0))))))

# #### Curvature
#
# Exact: the second-derivative family is computed in Float64 on both sides, and
# the formulas are the same formulas.
@test tdin(parent(laplacian(td_ras))) == tdin(GM.laplacian(td_z; cellsize=td_cell))
@test tdin(parent(profile_curvature(td_ras))) ==
      tdin(GM.profile_curvature(td_z; cellsize=td_cell))
@test tdin(parent(tangential_curvature(td_ras))) ==
      tdin(GM.tangential_curvature(td_z; cellsize=td_cell))
@test tdin(parent(plan_curvature(td_ras))) == tdin(GM.plan_curvature(td_z; cellsize=td_cell))
@test tdin(parent(laplacian(td_ras; gis=true))) == 100 .* tdin(parent(laplacian(td_ras)))

# A dilated window is still a derivative. On a quadratic surface the second
# derivative is constant, so widening the window must not change the answer —
# which it does in Geomorphometry, whose `scaled8nb(radius)` dilates the offsets
# but keeps dividing by the undilated cell size.
quad_z = [0.5x^2 for x in 0.0:1.0:14.0, _ in 1:15]
quad_ras = Raster(quad_z, (X(0.0:1.0:14.0), Y(0.0:1.0:14.0)))
@test parent(laplacian(quad_ras; radius=1))[8, 8] ≈ parent(laplacian(quad_ras; radius=2))[8, 8]
@test parent(laplacian(quad_ras; radius=1))[8, 8] ≈
      parent(laplacian(quad_ras; radius=3))[8, 8]
@test GM.laplacian(quad_z; radius=1)[8, 8] != GM.laplacian(quad_z; radius=2)[8, 8]
# The frame a radius-2 window cannot answer is two cells wide, not one.
@test all(isnan, parent(laplacian(quad_ras; radius=2))[2, :])
@test !any(isnan, parent(laplacian(quad_ras; radius=2))[3:end-2, 3:end-2])

# #### Illumination
#
# Geomorphometry rounds to `UInt8`, so agreement to half a level is exact
# agreement before the rounding.
@test all(abs.(tdin(parent(hillshade(td_ras))) .-
                Float64.(tdin(GM.hillshade(td_z; cellsize=td_cell)))) .<= 0.5)
@test all(abs.(tdin(parent(multihillshade(td_ras))) .-
                Float64.(tdin(GM.multihillshade(td_z; cellsize=td_cell)))) .<= 0.5)
@test all(x -> isnan(x) || 0.0 <= x <= 255.0, parent(hillshade(td_ras)))
# The azimuth is a real argument, not a decoration.
@test parent(hillshade(td_ras; azimuth=135.0)) != parent(hillshade(td_ras))
@test Rasters.name(hillshade(td_ras)) == :hillshade

# #### Storage order
#
# The two layouts Geomorphometry cannot distinguish. Every quantity above is
# geographic, so flipping the Y lookup must flip the output and nothing else.
td_yx = Raster(permutedims(td_z), (Y(range(0.0; step=td_dy, length=td_n)),
    X(range(0.0; step=td_dx, length=td_n))))
td_flip = Raster(reverse(td_z; dims=2), (X(range(0.0; step=td_dx, length=td_n)),
    Y(range((td_n - 1) * td_dy; step=-td_dy, length=td_n))))
for f in (slope, aspect, laplacian, plan_curvature, profile_curvature,
        tangential_curvature, hillshade, multihillshade, pssm)
    reference = parent(f(td_ras))
    @test isequal(permutedims(parent(f(td_yx))), reference)
    @test isequal(reverse(parent(f(td_flip)); dims=2), reference)
end

# #### Nodata
td_nan = copy(td_z)
td_nan[6, 6] = NaN
td_nan_ras = Raster(td_nan, (X(range(0.0; step=td_dx, length=td_n)),
    Y(range(0.0; step=td_dy, length=td_n))))
for f in (slope, aspect, laplacian, hillshade, plan_curvature)
    @test isnan(parent(f(td_nan_ras))[6, 6])
end
# The record estimators skip a nodata neighbor rather than propagating it.
@test isfinite(parent(slope(td_nan_ras; method=PlaneFit()))[5, 6])
@test isfinite(parent(aspect(td_nan_ras; method=PlaneFit()))[5, 6])

# #### Cell grids
#
# The windowed family refuses, naming the generic estimator; the record family
# runs.
for f in (laplacian, plan_curvature, profile_curvature, tangential_curvature)
    err = try f(cell_raster); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("rectilinear", sprint(showerror, err))
end
for m in (Horn(), ZevenbergenThorne())
    err = try aspect(cell_raster; method=m); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("PlaneFit", sprint(showerror, err))
end
@test defaultmethod(aspect, cell_grid) isa PlaneFit
@test defaultmethod(aspect, grid_xy) isa Horn
cell_aspect = aspect(cell_raster)
@test cell_aspect isa Raster && length(cell_aspect) == length(cell_raster)
@test all(a -> isnan(a) || 0.0 <= a < 360.0, parent(cell_aspect))
@test all(isfinite, parent(cell_aspect)) # A tilted sphere has an aspect everywhere
cell_shade = hillshade(cell_raster)
@test all(x -> isnan(x) || 0.0 <= x <= 255.0, parent(cell_shade))
@test count(isfinite, parent(cell_shade)) == length(cell_raster)
@test all(isfinite, parent(multihillshade(cell_raster)))
@test all(isfinite, parent(pssm(cell_raster)))
@test all(parent(slope(cell_raster; method=MDG())) .>=
        parent(steepest_slope(cell_raster)) .- 1e-9)
