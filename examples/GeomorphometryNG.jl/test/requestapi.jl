# ### The request API
#
# The canonical field order is independent of the order the request used, and
# a record carries exactly the requested fields — no more.
@test requestfields((Value(),)) === (Value(),)
@test requestfields((Bearing(), Value(), Index())) === (Index(), Value(), Bearing())
@test requestfields((Value(), Value())) === (Value(),)
@test requestfields(()) === ()
bad_needs = try requestfields((:value,)); nothing catch e; e end
@test bad_needs isa ArgumentError
bad_needs2 = try requestfields(Value()); nothing catch e; e end
@test bad_needs2 isa ArgumentError

record = neighborrecord(requestfields((Distance(), Value())), CartesianIndex(1, 1),
    3.0, (distance=2.0, bearing=90.0))
@test record === (value=3.0, distance=2.0)
@test keys(record) == (:value, :distance)
# (a) Reading a field that was not requested is a loud failure: Julia's own
# NamedTuple field error, naming the field and listing what is available.
missing_field = try record.bearing; nothing catch e; e end
@test missing_field isa Exception
@test occursin("bearing", sprint(showerror, missing_field))
missing_index = try record.index; nothing catch e; e end
@test missing_index isa Exception
@test occursin("index", sprint(showerror, missing_index))
# The same failure reaches a kernel: PlaneFit's request has no `index`.
kernel_field_error = try
    mapneighbors(raster_xy, grid_xy; needs=PLANEFIT_NEEDS) do _, value, neighbors
        sum(n -> n.index[1], neighbors)
    end
    nothing
catch e
    e
end
@test kernel_field_error isa Exception

# Named stencil fields refer to the same geographic directions in both storage
# layouts.
tagged(x, y) = 100.0x + y
tagged_xy = Raster([tagged(x, y) for x in xs, y in ys], (X(xs), Y(ys)))
tagged_yx = Raster(permutedims(parent(tagged_xy)), (Y(ys), X(xs)))
stencil_xy, stencil_yx = northupstencil(grid_xy), northupstencil(grid_yx)
hood_xy = Stencils.stencil(Stencils.StencilArray(tagged_xy, stencil_xy), (3, 2))
hood_yx = Stencils.stencil(Stencils.StencilArray(tagged_yx, stencil_yx), (2, 3))
@test all(name -> getproperty(hood_xy, name) == getproperty(hood_yx, name),
    keys(NORTH_UP_NEIGHBORS))
@test hood_xy.north == tagged(xs[3], ys[1])
@test hood_xy.east == tagged(xs[4], ys[2])

# Increasing the ring count expands the neighborhood without changing the
# callback. Index filtering removes out-of-bounds neighbors.
neighborcount(_, _, neighbors) = count(_ -> true, neighbors)
counts1_xy = mapneighbors(neighborcount, raster_xy, grid_xy)
@test eltype(counts1_xy) == Int # Boundary handling does not add `missing`
@test parent(counts1_xy)[1, 1] == 3
counts2_xy = mapneighbors(neighborcount, raster_xy, grid_xy, NeighborRings(2))
counts2_yx = mapneighbors(neighborcount, raster_yx, grid_yx, NeighborRings(2))
@test all(parent(counts2_xy) .>= parent(counts1_xy))
@test maximum(parent(counts2_xy)) > maximum(parent(counts1_xy))
@test parent(counts2_xy) == permutedims(parent(counts2_yx))

# A rectilinear geometry table is only built when distance or bearing is asked
# for; the default request leaves the payload a zero-size singleton.
@test neighborgeometry(grid_xy).payload.geometry isa NoGeometry
let withgeom = neighborgeometry(grid_xy, NeighborRings(1), (Value(), Distance()))
    @test withgeom.payload.geometry isa UniformGeometry
end

slope_xy = steepest_slope(raster_xy)
slope_yx = steepest_slope(raster_yx; spatialdims=(X, Y))
direction_xy = flow_direction(raster_xy)
direction_yx = flow_direction(raster_yx)

@test slope_xy isa Raster
# Boundary handling preserves the inferred output type.
@test eltype(slope_xy) == Float64
# The public method restores output metadata.
@test Rasters.name(slope_xy) == :steepest_slope
@test parent(slope_xy) ≈ permutedims(parent(slope_yx))
@test all(isequal.(parent(direction_xy), permutedims(parent(direction_yx))))
@test all(d -> isnan(d) || d ≈ 90.0, parent(direction_xy))

# Plain matrices use the same algorithms through an inferred rectilinear grid.
# An isbits padding value also supports integer input arrays.
matrix_slope = steepest_slope(data_xy; spacing=(2.0, -5.0))
@test matrix_slope isa Matrix{Float64} # Plain arrays produce plain-array outputs
@test matrix_slope ≈ parent(slope_xy)
int_slope = steepest_slope(Int.(data_xy); spacing=(2.0, -5.0))
@test int_slope ≈ matrix_slope

# Grid-construction keywords are consumed by `spatialparts`. Unknown algorithm
# keywords raise a `MethodError`, and unsupported manifolds raise an
# `ArgumentError` during grid construction.
bad_kwarg = try steepest_slope(raster_xy; bogus=1); nothing catch e; e end
@test bad_kwarg isa MethodError
geodesic_err = try spatialparts(raster_xy; manifold=Geodesic()); nothing catch e; e end
@test geodesic_err isa ArgumentError
# `Stencils.mapstencil` always threads, so the rectilinear driver has no
# `threaded` knob and says so by refusing the keyword.
no_threading = try steepest_slope(raster_xy; threaded=false); nothing catch e; e end
@test no_threading isa MethodError

# Horn and PlaneFit recover the same slope for an exact plane through their
# respective window and neighbor-record interfaces.
horn = slope(raster_xy) # Rectilinear grids use Horn by default
@test all(v -> v ≈ atand(3.0), parent(horn)[2:end-1, 2:end-1])
# A complete Horn window is unavailable at the border.
@test all(isnan, parent(horn)[1, :])
pf = slope(raster_xy; method=PlaneFit())
@test all(v -> v ≈ atand(3.0), parent(pf)) # Plane fitting uses available edge neighbors
pf_yx = slope(raster_yx; spatialdims=(X, Y), method=PlaneFit())
@test parent(pf) ≈ permutedims(parent(pf_yx))

# Every cell in this planar grid has the same area.
ca = cellarea(grid_xy)
@test ca[3, 2] == 10.0
@test sum(ca) == 10.0 * length(raster_xy)
@test cellarea(grid_xy, CartesianIndex(1, 1)) == 10.0

# Flow accumulation is the priority-flood sweep, and it returns two named
# products. All area reaches a cell with no downstream neighbor, which is an
# outlet; `missingval = nothing` on the direction raster is what keeps Rasters
# from cooking one up out of `typemax(FlowDirection)`.
acc, dirs = flowaccumulation(raster_xy)
@test Rasters.name(acc) == :flowaccumulation
@test Rasters.name(dirs) == :flowdirection
@test eltype(dirs) == FlowDirection{LDD,UInt8}
@test isnothing(Rasters.missingval(dirs))
@test isnan(Rasters.missingval(acc))
sweep_xy = floodsweep(raster_xy, grid_xy)
@test sum(parent(acc)[sweep_xy.down .== 0]) ≈ sum(cellarea(grid_xy))
@test all(parent(acc) .>= 10.0) # Every cell carries at least its own area
@test count(ispit, parent(dirs)) == sweep_xy.nseeds
