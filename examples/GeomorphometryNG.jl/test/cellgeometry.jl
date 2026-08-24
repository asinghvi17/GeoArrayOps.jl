# ### (b) A value-only cell request builds no centroid table
#
# The payload for the default request is a zero-size singleton: there is no
# centroid storage to be found, at any grid size. When geometry *is* requested,
# the on-demand provider holds only the level grid and the radius.
ncells = length(cell_raster)
value_geom = neighborgeometry(cell_grid, NeighborRings(1), (Value(),))
ondemand_geom = neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
@test value_geom.payload isa NoCentroids
@test sizeof(value_geom.payload) == 0
@test ondemand_geom.payload isa OnDemandCentroids
@test _streamable(value_geom.fields) isa Requested   # -> DGG pass=Values()
@test _streamable(ondemand_geom.fields) isa Absent   # -> DGG pass=Neighbors()
neighborgeometry(cell_grid, NeighborRings(1), (Value(),))            # warm up
neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
value_geom_bytes = @allocated neighborgeometry(cell_grid, NeighborRings(1), (Value(),))
ondemand_bytes = @allocated neighborgeometry(cell_grid, NeighborRings(1), STEEPEST_NEEDS)
# A materialized centroid table is 24 bytes per cell (three Float64 direction
# cosines). Both on-demand paths stay orders of magnitude below that.
@test value_geom_bytes < 24 * ncells / 8
@test ondemand_bytes < 24 * ncells / 8
precompute(ondemand_geom) # warm up
precomputed_bytes = @allocated precompute(ondemand_geom)
@test precomputed_bytes >= 24 * ncells # The opt-in really does materialize
