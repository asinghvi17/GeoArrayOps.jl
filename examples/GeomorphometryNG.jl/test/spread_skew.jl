# ### Cost distance and the distribution-only algorithm
#
# The two ends of the family map: `spread` is the traversal that needs geometry
# *and* positions at once, and `skewness_balancing` needs no grid at all.

# #### spread
#
# Uniform friction on a unit grid makes the cost the chebyshev-ish path length,
# and Geomorphometry's Tomlin is the same Dijkstra — exactly, in Float64.
sp_n = 15
sp_friction = fill(1.0, sp_n, sp_n)
sp_ras = Raster(sp_friction, (X(1.0:sp_n), Y(1.0:sp_n)))
sp_points = falses(sp_n, sp_n)
sp_points[8, 8] = true
sp_ng = spread(sp_points, 0.0, sp_ras)
@test parent(sp_ng) == GM.spread(sp_points, 0.0, sp_friction; cellsize=(1.0, 1.0))
@test parent(sp_ng)[8, 8] == 0.0
@test parent(sp_ng)[8, 9] ≈ 1.0
@test parent(sp_ng)[9, 9] ≈ sqrt(2)
@test Rasters.name(sp_ng) == :spread

# Variable friction, two sources, and a nonzero initial value — the shape of a
# real call.
sp_var = [1.0 + 0.5 * abs(i - 8) for i in 1:sp_n, _ in 1:sp_n]
sp_var_ras = Raster(sp_var, (X(1.0:sp_n), Y(1.0:sp_n)))
sp_two = falses(sp_n, sp_n)
sp_two[2, 2] = true
sp_two[14, 14] = true
sp_initial = fill(5.0, sp_n, sp_n)
# Six of the 225 cells differ in the last bit: the two queues break equal-cost
# ties in a different order, so the additions happen in a different order.
@test parent(spread(sp_two, sp_initial, sp_var_ras)) ≈
      GM.spread(findall(sp_two), sp_initial, sp_var; cellsize=(1.0, 1.0))
@test parent(spread(sp_two, 5.0, sp_var_ras)) == parent(spread(sp_two, sp_initial, sp_var_ras))
@test parent(spread(findall(sp_two), sp_initial, sp_var_ras)) ==
      parent(spread(sp_two, sp_initial, sp_var_ras))
# Unreached cells hold `limit`, and a source list that reaches nothing leaves
# every cell there.
sp_none = spread(falses(sp_n, sp_n), 0.0, sp_ras; limit=-1.0)
@test all(parent(sp_none) .== -1.0)
# Cell size is real: doubling the spacing doubles every distance.
@test parent(spread(sp_points, 0.0, sp_friction; spacing=(2.0, 2.0))) ≈
      2 .* parent(sp_ng)
# Storage order does not move a source.
sp_flip = Raster(reverse(sp_friction; dims=2),
    (X(1.0:sp_n), Y(Float64(sp_n):-1.0:1.0)))
@test reverse(parent(spread(reverse(sp_points; dims=2), 0.0, sp_flip)); dims=2) ≈
      parent(sp_ng)
# A mask of the wrong size is refused.
@test (try spread(falses(3, 3), 0.0, sp_ras); nothing catch e; e end) isa ArgumentError
# Only Tomlin is implemented, and the refusal says why the other two are not.
sp_method_error = try spread(sp_points, 0.0, sp_ras; method=:eastman)
    nothing catch e; e end
@test sp_method_error isa ArgumentError
@test occursin("Tomlin", sprint(showerror, sp_method_error))

# On a cell grid the same traversal runs off `eachneighbor`, with `cellkey` and
# `storageposition` bridging the queue's positions back to cell handles.
sp_cell_friction = Raster(fill(1.0, length(cell_raster)), (DGG.Cells(cell_lookup),))
sp_cell_points = falses(length(cell_raster))
sp_cell_points[1] = true
sp_cell = spread(sp_cell_points, 0.0, sp_cell_friction)
@test all(isfinite, parent(sp_cell))
@test parent(sp_cell)[1] == 0.0
@test all(>=(0.0), parent(sp_cell))
# The cost of the first ring is the great-circle distance to those centroids.
sp_cell_geom = neighborgeometry(cell_grid, NeighborRings(1), (Index(), Distance()))
sp_first_ring = collect(eachneighbor(sp_cell_geom, sp_cell_friction,
    first(cellindices(sp_cell_friction, cell_grid))))
@test all(sp_first_ring) do n
    isapprox(parent(sp_cell)[DGG.localindex(n.index)], n.distance; rtol=1e-9)
end
# The farthest cell is about half a great circle away, at friction 1.
# The path is a chain of centroid-to-centroid hops, so it is longer than the
# great circle it approximates but of the same order.
@test maximum(parent(sp_cell)) < 1.5 * π * AUTHALIC_RADIUS_M
@test maximum(parent(sp_cell)) > 0.4 * π * AUTHALIC_RADIUS_M

# #### skewness_balancing
#
# The only algorithm here that reads no neighborhood at all, so it is identical
# on both backends by construction and exact against Geomorphometry.
@test parent(skewness_balancing(tf_ras)) == GM.skewness_balancing(tf_z)
@test eltype(parent(skewness_balancing(tf_ras))) == Bool
@test Rasters.name(skewness_balancing(tf_ras)) == :skewness_balancing
# A mask of a monotone ramp keeps the low end and drops the high end.
sk_ramp = Raster(Float64.(reshape(1:400, 20, 20)), (X(1.0:20.0), Y(1.0:20.0)))
sk_mask = parent(skewness_balancing(sk_ramp))
@test sk_mask == GM.skewness_balancing(Float64.(reshape(1:400, 20, 20)))
@test count(sk_mask) > 0
# Non-finite cells sort to the top of the search window, so they end up on the
# object side of any threshold that falls below the top — including this
# positively skewed fixture, where the threshold really does fall inside.
sk_skewed = reshape(vcat(fill(1.0, 360), collect(50.0:1.0:89.0)), 20, 20)
sk_skewed[5] = NaN
sk_nan_mask = parent(skewness_balancing(Raster(sk_skewed, (X(1.0:20.0), Y(1.0:20.0)))))
@test !sk_nan_mask[5]
@test sk_nan_mask == GM.skewness_balancing(sk_skewed)
@test 0 < count(sk_nan_mask) < length(sk_skewed)
# On a cell grid it runs unchanged and agrees with the same call on the bare
# value vector — the point being that the grid never entered.
@test parent(skewness_balancing(cell_raster)) == GM.skewness_balancing(cell_values)
