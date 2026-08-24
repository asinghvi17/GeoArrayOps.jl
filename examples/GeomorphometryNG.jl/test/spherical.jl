# An EPSG:4326 Raster uses spherical geometry. East-west distances and cell
# areas vary by latitude row.
lons = 0.0:1.0:5.0
lats = 60.0:-1.0:55.0
geo_data = [1000.0 - 10.0 * lon for lon in lons, lat in lats]
geo_raster = Raster(geo_data, (X(lons), Y(lats)); crs=EPSG(4326))
_, grid_geo = spatialparts(geo_raster)
@test grid_geo.manifold isa Spherical

geom_geo = neighborgeometry(grid_geo, NeighborRings(1), (Value(), Distance(), Bearing()))
@test geom_geo.payload.geometry isa RowGeometry
east = findfirst(==(:east), collect(keys(NORTH_UP_NEIGHBORS)))
t60 = geometryat(geom_geo, CartesianIndex(3, 1))
t55 = geometryat(geom_geo, CartesianIndex(3, 6))
@test t60[east].bearing == 90.0
@test t60[east].distance ≈ deg2rad(1.0) * cosd(60.0) * AUTHALIC_RADIUS_M
@test t55[east].distance > t60[east].distance # Longitude spacing widens toward the equator

slope_geo = steepest_slope(geo_raster)
@test parent(slope_geo)[3, 1] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(60.0) * AUTHALIC_RADIUS_M))
@test parent(slope_geo)[3, 6] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(55.0) * AUTHALIC_RADIUS_M))
horn_geo = slope(geo_raster)
@test parent(horn_geo)[3, 2] ≈
    atand(10.0 / (deg2rad(1.0) * cosd(59.0) * AUTHALIC_RADIUS_M))

ca_geo = cellarea(grid_geo)
@test ca_geo[1, 1] ≈ AUTHALIC_RADIUS_M^2 * deg2rad(1.0) * (sind(60.5) - sind(59.5))
@test ca_geo[1, 6] > ca_geo[1, 1]
