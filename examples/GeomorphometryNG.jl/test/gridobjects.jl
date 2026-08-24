# ## Behavior checks
#
# Store the same plane in X,Y and Y,X axis order. The reversed Y coordinates and
# unequal X and Y spacing expose calculations that incorrectly depend on
# storage-axis order.

xs = 0.0:2.0:8.0
ys = 20.0:-5.0:5.0
surface(x, y) = 100.0 - 3.0x

data_xy = [surface(x, y) for x in xs, y in ys]
raster_xy = Raster(data_xy, (X(xs), Y(ys)))
raster_yx = Raster(permutedims(data_xy), (Y(ys), X(xs)))

_, grid_xy = spatialparts(raster_xy)
_, grid_yx = spatialparts(raster_yx)
@test axismap(grid_xy) == (1, 2)
@test axismap(grid_yx) == (2, 1)
# A Raster without a CRS defaults to planar geometry.
@test grid_xy.manifold isa Planar
@test collect(cellindices(raster_xy, grid_xy)) == collect(CartesianIndices(raster_xy))
