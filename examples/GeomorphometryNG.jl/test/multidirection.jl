# ### Multi-direction flow: DInf and FD8
#
# `down` carries one downstream position per cell, which is the whole of D8.
# DInf and FD8 split a cell's flow, so the sweep result grows a *ragged*
# partition beside `down` — an addition, not a replacement, which is what
# `FloodSweep`'s opacity was for.
#
# A separate, smaller tie-free fixture keeps the two Geomorphometry runs cheap.

md_n = 40
md_z = [100.0 * splitmix01(i + md_n * (j - 1)) + 1e-9 * (i + md_n * (j - 1))
        for i in 1:md_n, j in 1:md_n]
@test length(unique(md_z)) == length(md_z)
# North-up and ascending in both axes, so Geomorphometry's `_cellbearing` and
# its Horn aspect agree with this package's geographic ones.
md_ras = Raster(md_z, (X(range(0.0; step=1.0, length=md_n)),
    Y(range(0.0; step=1.0, length=md_n))))
_, md_grid = spatialparts(md_ras)

md_acc8, md_dir8 = flowaccumulation(md_ras)
md_accinf, md_dirinf = flowaccumulation(md_ras; method=DInf())
md_accfd, md_dirfd = flowaccumulation(md_ras; method=FD8())
gm_accinf, gm_dirinf = GM.flowaccumulation(md_z; method=DInf(), cellsize=(1.0, 1.0))
gm_accfd, gm_dirfd = GM.flowaccumulation(md_z; method=FD8(), cellsize=(1.0, 1.0))

# #### DInf: exact parity, directions and mass
#
# The direction raster is `D8D` on a rectilinear grid, so it is bit-comparable
# with Geomorphometry.
@test eltype(md_dirinf) == FlowDirection{D8D,UInt8}
@test Int.(parent(md_dirinf)) == Int.(gm_dirinf)
# Geomorphometry accumulates in Float32 over the same addition sequence.
@test maximum(abs.(parent(md_accinf) .- Float64.(gm_accinf))) / maximum(parent(md_accinf)) <
      1e-5
@test count(d -> GM.ndirections(d) > 1, parent(md_dirinf)) > 100 # It really does split

# #### FD8: the same direction set, a different contour length
#
# Geomorphometry's rectilinear FD8 uses a hard-coded table of contour lengths
# that is transposed with respect to δx/δy and twice Quinn's cardinal value; the
# port uses the angular rule its non-matrix path defines, which is the one a
# hexagonal ring can also answer. The *set* of receiving neighbors is identical,
# so that is what is compared.
@test eltype(md_dirfd) == FlowDirection{D8D,UInt8}
@test Int.(parent(md_dirfd)) == Int.(gm_dirfd)
@test maximum(abs.(parent(md_accfd) .- Float64.(gm_accfd))) / maximum(parent(md_accfd)) >
      1e-3
@test count(d -> GM.ndirections(d) > 1, parent(md_dirfd)) > 100

# #### The partition itself
md_part = flowpartition(DInf(), md_ras)
md_partfd = flowpartition(FD8(), md_ras)
md_part8 = flowpartition(D8(), md_ras)
_, md_g = spatialparts(md_ras)
md_sweep = floodsweep(md_ras, md_g)
for part in (md_part, md_partfd, md_part8)
    # Every cell either drains (weights summing to one) or is an outlet.
    @test all(eachindex(md_z)) do p
        r = partitionrange(part, p)
        isempty(r) ? md_sweep.down[p] == 0 : sum(part.weights[r]) ≈ 1.0
    end
    @test all(w -> 0.0 < w <= 1.0, part.weights)
    # Mass conservation: everything reaches a cell with no outgoing edge.
    local md_acc = Vector{Float64}(undef, length(md_z))
    _fillcellareas!(md_acc, md_g)
    _accumulateweighted!(md_acc, md_sweep.order, part)
    local md_sinks = [p for p in eachindex(md_z) if isempty(partitionrange(part, p))]
    @test sum(md_acc[md_sinks]) ≈ sum(cellarea(md_g))
end
# The degenerate partition is D8: one edge of weight one, pointing where `down`
# points, and accumulating to the same numbers as the dedicated pass.
@test all(p -> md_sweep.down[p] == 0 ||
        (length(partitionrange(md_part8, p)) == 1 &&
         md_part8.targets[first(partitionrange(md_part8, p))] == md_sweep.down[p]),
    eachindex(md_z))
md_acc8b = Vector{Float64}(undef, length(md_z))
_fillcellareas!(md_acc8b, md_g)
_accumulateweighted!(md_acc8b, md_sweep.order, md_part8)
@test reshape(md_acc8b, size(md_z)) ≈ parent(md_acc8)
# Where DInf does not split, it is D8.
md_unsplit = [p for p in eachindex(md_z) if length(partitionrange(md_part, p)) == 1]
@test !isempty(md_unsplit)
@test all(p -> md_part.targets[first(partitionrange(md_part, p))] == md_sweep.down[p],
    md_unsplit)

# #### Storage order, and the façade
md_flip = Raster(reverse(md_z; dims=2), (X(range(0.0; step=1.0, length=md_n)),
    Y(range((md_n - 1) * 1.0; step=-1.0, length=md_n))))
@test reverse(parent(flowaccumulation(md_flip; method=DInf())[1]); dims=2) ≈
      parent(md_accinf)
@test reverse(Int.(parent(flowaccumulation(md_flip; method=DInf())[2])); dims=2) ==
      Int.(parent(md_dirinf))
@test (Rasters.name(md_accinf), Rasters.name(md_dirinf)) ==
      (:flowaccumulation, :flowdirection)
@test isnothing(Rasters.missingval(md_dirinf))
@test Int.(parent(flowdirection(md_ras; method=DInf()))) == Int.(parent(md_dirinf))
# Anything that is not one of the three named methods still refuses clearly.
struct MDTestMethod <: GM.FlowDirectionMethod end
md_error = try flowaccumulation(md_ras; method=MDTestMethod()); nothing catch e; e end
@test md_error isa ArgumentError
@test occursin("D8(), DInf() and FD8()", sprint(showerror, md_error))

# #### Cell grids
#
# Geomorphometry's own DInf crashes on a cell raster (it builds a
# `similar(::Raster, FlowDirection)` and Rasters cooks up a `typemax` for it),
# so there is no baseline. The invariants are the baseline.
sub_accinf, sub_dirinf = flowaccumulation(sub_raster; method=DInf())
sub_accfd, sub_dirfd = flowaccumulation(sub_raster; method=FD8())
sub_partinf = flowpartition(DInf(), sub_raster)
# A hexagonal ring has no numpad, so the set-valued direction is a bitmask over
# ring slots — the split-flow analogue of `RingSlot`.
@test eltype(sub_dirinf) == FlowDirection{RingMask,UInt8}
@test count(ispit, sub_dirinf) == sub_sweep.nseeds
@test all(eachindex(sub_vals)) do p
    r = partitionrange(sub_partinf, p)
    isempty(r) ? sub_sweep.down[p] == 0 : sum(sub_partinf.weights[r]) ≈ 1.0
end
@test count(p -> length(partitionrange(sub_partinf, p)) > 1, eachindex(sub_vals)) > 10
sub_sinks = [p for p in eachindex(sub_vals) if isempty(partitionrange(sub_partinf, p))]
@test sum(parent(sub_accinf)[sub_sinks]) ≈ sum(cellarea(sub_grid))
sub_partfd = flowpartition(FD8(), sub_raster)
@test sum(parent(sub_accfd)[[p for p in eachindex(sub_vals)
    if isempty(partitionrange(sub_partfd, p))]]) ≈ sum(cellarea(sub_grid))
# Every bit set names a real ring slot of that cell.
@test all(eachindex(sub_vals)) do p
    code = Int(parent(sub_dirinf)[p])
    all(k -> k <= nslots(sub_table, p),
        [k for k in 1:8 if (code >> (k - 1)) & 1 == 1])
end
# A one-slot mask shows as that slot's circled digit; a split shows as `✳`.
@test sprint(show, first(filter(ispit, parent(sub_dirinf)))) == "·"
@test sprint(show, FlowDirection{RingMask,UInt8}(0x04)) == "③"
@test sprint(show, FlowDirection{RingMask,UInt8}(0x05)) == "✳"
@test GM.ismulti(RingMask)
# The nodata convention: a closed cell neither receives nor emits.
sub_inf_closed = flowpartition(DInf(), sub_raster; closed=isnan.(sub_nan))
@test all(p -> isempty(partitionrange(sub_inf_closed, p)), findall(isnan, sub_nan))
@test all(eachindex(sub_vals)) do p
    all(e -> !isnan(sub_nan[Int(sub_inf_closed.targets[e])]),
        partitionrange(sub_inf_closed, p))
end
