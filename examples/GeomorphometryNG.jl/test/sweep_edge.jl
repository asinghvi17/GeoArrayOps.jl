# #### (10) Degenerate grids
one_acc, one_dirs = flowaccumulation(fill(5.0, 1, 1))
@test one_acc == fill(1.0, 1, 1) # Its own outlet, carrying its own area
@test ispit(only(one_dirs))
@test only(floodsweep(fill(5.0, 1, 1), last(spatialparts(fill(5.0, 1, 1)))).order) == 1
empty_acc, empty_dirs = flowaccumulation(zeros(0, 0))
@test isempty(empty_acc) && isempty(empty_dirs)

# #### (12) Multi-direction methods
#
# These were refused while the sweep carried one downstream position per cell.
# They are now served by the ragged partition beside it (see
# `multidirection.jl`); what is checked here is that the façade routes them and
# that they still conserve area on the fixture the D8 assertions above use.
dinf_acc, dinf_dirs = flowaccumulation(asxy(bowl_z); method=DInf())
@test eltype(dinf_dirs) == FlowDirection{D8D,UInt8}
@test sum(parent(dinf_acc)[Int.(parent(dinf_dirs)) .== 0]) ≈ sum(cellarea(last(spatialparts(asxy(bowl_z)))))
fd8_acc, fd8_dirs = flowaccumulation(asxy(bowl_z); method=FD8())
@test eltype(fd8_dirs) == FlowDirection{D8D,UInt8}
@test sum(parent(fd8_acc)[Int.(parent(fd8_dirs)) .== 0]) ≈ sum(cellarea(last(spatialparts(asxy(bowl_z)))))

# #### (13) Round-trip through the façade
matrix_acc, matrix_dirs = flowaccumulation(bowl_z; spacing=(2.0, -5.0))
@test matrix_acc isa Matrix{Float64}
@test matrix_dirs isa Matrix{FlowDirection{LDD,UInt8}}
raster_acc, raster_dirs = flowaccumulation(asxy(bowl_z))
@test raster_acc isa Raster && raster_dirs isa Raster
@test (Rasters.name(raster_acc), Rasters.name(raster_dirs)) ==
        (:flowaccumulation, :flowdirection)
@test isnan(Rasters.missingval(raster_acc))
@test isnothing(Rasters.missingval(raster_dirs))
@test Rasters.name(settle(asxy(bowl_z))) == :settle
@test Rasters.name(flowdirection(asxy(bowl_z))) == :flowdirection
# Algorithm keywords reach the inner method through the façade unchanged.
@test isnan(parent(settle(Raster(nan_z, (X(1.0:8.0), Y(8.0:-1.0:1.0)));
    closed=isnan.(nan_z)))[4, 4])
# Rule A: below the façade every output is a plain array. Rasters' missingval
# machinery is reached exactly once, in `rebuildoutput`.
inner_acc, inner_dirs = flowaccumulation(asxy(bowl_z), last(spatialparts(asxy(bowl_z))))
@test inner_acc isa Matrix{Float64} && !(inner_acc isa Raster)
@test inner_dirs isa Matrix{FlowDirection{LDD,UInt8}} && !(inner_dirs isa Raster)
# A geographic raster keeps its CRS, and its per-row cell areas still conserve.
geo_acc, geo_dirs = flowaccumulation(geo_raster)
@test Rasters.crs(geo_acc) == Rasters.crs(geo_raster)
@test sum(vec(parent(geo_acc))[floodsweep(geo_raster, grid_geo).down .== 0]) ≈
        sum(cellarea(grid_geo))
