# A DGG Raster uses the same public algorithms through `CellGrid`. Elevation is
# proportional to each centroid's Z coordinate, so most cells have a downhill
# neighbor.
level = 2
dgg = DGG.levelgrid(DGG.IGeo7System(), level)
cell_lookup = DGG.CellLookup(dgg)
cell_values = [10_000.0 * DGG.cell_centroid(dgg, cell)[3] for cell in cell_lookup]
cell_raster = Raster(cell_values, (DGG.Cells(cell_lookup),))

_, cell_grid = spatialparts(cell_raster; spatialdims=DGG.Cells)
@test cell_grid isa CellGrid
@test cell_grid.cells === cell_lookup.cells
@test manifold(cell_grid) isa Spherical
@test manifold(cell_grid).radius == AUTHALIC_RADIUS_M

# A user-provided radius rescales every derived quantity; on a unit sphere the
# cell areas sum to the full sphere's solid angle. Non-spherical manifolds are
# rejected.
_, unit_grid = spatialparts(cell_raster; manifold=Spherical(; radius=1.0))
@test manifold(unit_grid).radius == 1.0
@test sum(cellarea(unit_grid)) ≈ 4π
planar_cells_err = try spatialparts(cell_raster; manifold=Planar()); nothing catch e; e end
@test planar_cells_err isa ArgumentError

cell_index = first(cellindices(cell_raster, cell_grid))
@test cell_index isa DGG.SubsetIndexedCell
@test cell_raster[cell_index] == first(cell_values)

cell_counts2 = mapneighbors(neighborcount, cell_raster, cell_grid, NeighborRings(2))
@test parent(cell_counts2)[1] == length(DGG.neighbors(cell_lookup, 1, 2))

cell_areas = cellarea(cell_grid)
@test sum(cell_areas) ≈ 4π * AUTHALIC_RADIUS_M^2 # IGeo7 cells cover the sphere
@test cellarea(cell_grid, cell_index) == cell_areas[1]

cell_slope = steepest_slope(cell_raster)
cell_direction = flow_direction(cell_raster)
@test cell_slope isa Raster
@test all(>=(0.0), parent(cell_slope))
@test any(isfinite, parent(cell_direction))
@test all(d -> isnan(d) || 0.0 <= d < 360.0, parent(cell_direction))

pf_cells = slope(cell_raster) # Cell grids use PlaneFit by default
@test all(isfinite, parent(pf_cells))
horn_err = try slope(cell_raster; method=Horn()); nothing catch e; e end
@test horn_err isa ArgumentError

# The same local statistics run on the cell backend, through DGG's streaming
# value pass, and reproduce the same hand-checkable relationships.
cell_tpi = topographic_position_index(cell_raster)
cell_tri = terrain_ruggedness_index(cell_raster)
cell_rough = roughness(cell_raster)
@test cell_tpi isa Raster && length(cell_tpi) == length(cell_raster)
@test Rasters.name(cell_rough) == :roughness
@test all(>=(0.0), parent(cell_rough))
@test all(isfinite, parent(cell_tpi))
# TRI with no normalization is the root of the summed squares of the same
# differences roughness takes the maximum of, so it is never the smaller one.
@test all(parent(cell_tri) .>= parent(cell_rough) .- 1e-9)
# `threaded` and `order` reach DGG.mapneighbors; the answer must not depend on
# either.
@test parent(topographic_position_index(cell_raster, cell_grid; threaded=false)) ==
        parent(cell_tpi)
@test parent(mapneighbors(_tpi_kernel, cell_raster, cell_grid, NeighborRings();
    needs=LOCAL_NEEDS, order=reverse(1:length(cell_raster)))) == parent(cell_tpi)
