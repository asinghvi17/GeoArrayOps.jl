# ### (c) `precompute` is a performance choice, not a semantic one
cell_geom_steep = neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
cell_geom_fit = neighborgeometry(cell_grid, NeighborRings(1), PLANEFIT_NEEDS)
@test precompute(cell_geom_steep).payload isa StoredCentroids
@test precompute(cell_geom_steep).fields === cell_geom_steep.fields
@test all(isequal.(
    parent(mapneighbors(_steepest_slope, cell_raster, cell_geom_steep)),
    parent(mapneighbors(_steepest_slope, cell_raster, precompute(cell_geom_steep)))))
@test all(isequal.(
    parent(mapneighbors(_planefit_slope, cell_raster, cell_geom_fit)),
    parent(mapneighbors(_planefit_slope, cell_raster, precompute(cell_geom_fit)))))
# `precompute` on a rectilinear grid is a no-op: its tables are already O(rows).
@test precompute(geom_geo).payload.geometry === geom_geo.payload.geometry

# A NaN cell propagates on the cell backend too.
nan_cell_values = copy(cell_values)
nan_cell_values[5] = NaN
nan_cell_raster = Raster(nan_cell_values, (DGG.Cells(cell_lookup),))
@test isnan(parent(steepest_slope(nan_cell_raster))[5])
@test isnan(parent(flow_direction(nan_cell_raster))[5])
@test isnan(parent(slope(nan_cell_raster))[5])
# A cell that is not a neighbor of the NaN cell is untouched.
untouched = findfirst(p -> p != 5 && !(5 in DGG.neighbors(cell_lookup, p, 1)),
    1:length(cell_values))
@test parent(steepest_slope(nan_cell_raster))[untouched] ==
        parent(cell_slope)[untouched]

# `eachneighbor` honors the request identically on both backends.
rect_records = collect(eachneighbor(
    neighborgeometry(grid_xy, NeighborRings(1), (Index(), Value())),
    raster_xy, CartesianIndex(1, 1)))
@test length(rect_records) == 3
@test all(r -> keys(r) == (:index, :value), rect_records)
cell_records = collect(eachneighbor(value_geom, cell_raster, cell_index))
@test all(r -> keys(r) == (:value,), cell_records)
@test length(cell_records) == length(DGG.neighbors(cell_lookup, 1, 1))
cell_records_full = collect(eachneighbor(
    neighborgeometry(cell_grid, NeighborRings(1), PLANEFIT_NEEDS), cell_raster, cell_index))
@test all(r -> keys(r) == (:value, :distance, :bearing), cell_records_full)
@test all(r -> r.distance > 0 && 0 <= r.bearing < 360, cell_records_full)
