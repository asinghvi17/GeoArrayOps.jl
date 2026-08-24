# #### Fixture 3 — cell grids
#
# `cell_raster` is a *complete* level-2 grid, and a sphere has no edges: the
# boundary seed set is empty, so the no-outlet fallback must fire rather than let
# an empty queue return `acc == cellarea` everywhere.
sphere_table = neighbortable(cell_grid)
@test !any(p -> isboundary(sphere_table, p), 1:length(sphere_table))
@test isempty(_boundaryseeds(sphere_table, falses(length(cell_raster))))
sphere_sweep = @test_logs (:warn, r"no domain boundary") match_mode=:any begin
    floodsweep(cell_raster, cell_grid; table=sphere_table)
end
@test length(sphere_sweep.order) == length(cell_raster)
sphere_acc = let a = Vector{Float64}(undef, length(cell_raster))
    _fillcellareas!(a, cell_grid)
    _accumulatedown!(a, sphere_sweep.order, sphere_sweep.down)
end
@test sum(sphere_acc[sphere_sweep.down .== 0]) ≈ sum(cellarea(cell_grid)) ≈
        4π * AUTHALIC_RADIUS_M^2
# A one-ring table is the contract; a wider request is refused, not clipped.
wide_table_error = try neighbortable(cell_grid, NeighborRings(2)); nothing catch e; e end
@test wide_table_error isa ArgumentError
# The adjacency table is built once per public call and threaded through the
# flood and the direction encoding — a second build would double this.
silently() do
    neighbortable(cell_grid)
    flowaccumulation(cell_raster, cell_grid) # Warm up before measuring
    @test (@allocated flowaccumulation(cell_raster, cell_grid)) <
            2 * (@allocated neighbortable(cell_grid))
end

# A *region* of a cell grid has a rim, so Geomorphometry can run on it too. The
# values carry a per-position ramp to make the flood tree unique.
sub_region = DGG.subtree(DGG.IGeo7System(),
    DGG.cellindex(DGG.levelgrid(DGG.IGeo7System(), 0), 5), 3)
sub_lookup = DGG.CellLookup(sub_region)
sub_cells = sub_lookup.cells
sub_vals = [10_000.0 * DGG.cell_centroid(sub_region, c)[3] + 1e-6 * p
            for (p, c) in enumerate(sub_cells)]
@test length(unique(sub_vals)) == length(sub_vals)
sub_raster = Raster(sub_vals, (DGG.Cells(sub_lookup),))
_, sub_grid = spatialparts(sub_raster; spatialdims=DGG.Cells)
sub_table = neighbortable(sub_grid)

# The seed set read off the zero slots is exactly DGG's own rim, and exactly the
# rim Geomorphometry rebuilds per cell from `neighborcount(complete, cv[p])`.
@test _boundaryseeds(sub_table, falses(length(sub_cells))) ==
        sort(collect(DGG.border(sub_cells)))
# The invariant the IGeo7 codec rests on: a complete-width row preserves slot
# identity, and on IGeo7 slot `k` is direction code `k`. That is why no ring is
# rebuilt and no relative cell is constructed.
@test all(1:length(sub_cells)) do p
    all(((k, q),) -> q == 0 || DGG.directioncode(sub_cells[q] - sub_cells[p]) == k,
        slots(sub_table, p))
end

# `nslots` is the complete degree, which is what makes a slot code portable
# between regions: a clipped row would renumber it.
@test nslots(tf_table, 1) == 8 == length(collect(slots(tf_table, 1)))
@test all(p -> nslots(sub_table, p) == length(collect(slots(sub_table, p))),
    1:length(sub_cells))
@test Set(nslots(sub_table, p) for p in 1:length(sub_cells)) ⊆ Set((5, 6))

sub_sweep, sub_acc, sub_dirs = sweepparts(sub_raster, sub_grid)
gm_sub_acc, gm_sub_dirs = GM.flowaccumulation(sub_raster; method=D8())
@test eltype(sub_dirs) == FlowDirection{LDD,UInt8}
@test Int.(sub_dirs) == Int.(parent(gm_sub_dirs))
@test maximum(abs.(sub_acc .- Float64.(parent(gm_sub_acc)))) / maximum(sub_acc) < 1e-5
@test sum(sub_acc[sub_sweep.down .== 0]) ≈ sum(cellarea(sub_grid))
@test count(ispit, sub_dirs) == sub_sweep.nseeds == length(collect(DGG.border(sub_cells)))

# Geomorphometry has no `filldepressions` for cell grids, so `settle` is checked
# against its defining invariants instead.
sub_settled = settle(sub_raster, sub_grid; sweep=sub_sweep)
@test all(sub_settled .>= sub_vals)
@test all(eachindex(sub_vals)) do p
    d = sub_sweep.down[p]
    d == 0 ? sub_settled[p] == sub_vals[p] :
    sub_settled[p] == max(sub_vals[p], sub_settled[d])
end

# The generic seam. A cell system with no relative-cell arithmetic gets the ring
# slot itself, which decodes back to a position in O(1) off the same row — no
# ring rebuild and no `findfirst`. It is exercised here through the IGeo7 table
# because slot `k` is what the IGeo7 codec is built on in the first place.
slot_dirs = Vector{FlowDirection{RingSlot,UInt8}}(undef, length(sub_cells))
_cellencode!(slot_dirs, sub_sweep.down, sub_table, FlowDirection{RingSlot,UInt8})
@test all(p -> downstreamposition(sub_table, p, slot_dirs[p]) == sub_sweep.down[p],
    eachindex(sub_cells))
@test count(ispit, slot_dirs) == sub_sweep.nseeds
@test Int.(sub_dirs) == [IGEO7_TO_LDD[Int(d) + 1] for d in slot_dirs]
@test sprint(show, first(filter(ispit, slot_dirs))) == "·"
@test sprint(show, first(filter(!ispit, slot_dirs))) != "·"

# (5) The tutorial's nodata overhang, in miniature: nodata cells sitting on the
# region rim, excluded with `closed`. This is the configuration in which
# Geomorphometry seeds rim cells it was told were closed and then writes past the
# end of `order`.
sub_nan = copy(sub_vals)
sub_nan[first(collect(DGG.border(sub_cells)), 10)] .= NaN
sub_nan_raster = Raster(sub_nan, (DGG.Cells(sub_lookup),))
sub_nan_sweep, sub_nan_acc, _ = sweepparts(sub_nan_raster, sub_grid; closed=isnan.(sub_nan))
@test length(sub_nan_sweep.order) == length(sub_nan) - count(isnan, sub_nan)
@test all(p -> sub_nan_sweep.down[p] == 0, findall(isnan, sub_nan))
