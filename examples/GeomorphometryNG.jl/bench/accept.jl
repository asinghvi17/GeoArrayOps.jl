# Acceptance benchmark for GeomorphometryNG on a DGG cell grid (measurement only).
#   julia --project=examples/GeomorphometryNG.jl/bench accept.jl <level> [mode]
# mode = "full" (default, the whole matrix) or "sub" (rows 1, 2, 4 only).
#
# Everything runs in ONE session per level so JIT / thread state is shared.
include(joinpath(@__DIR__, "setup.jl"))
using Printf, Statistics
import Geomorphometry as GM
using GeomorphometryNG
const NG = GeomorphometryNG

const LEVEL = parse(Int, get(ARGS, 1, "13"))
const MODE  = get(ARGS, 2, "full")

# min of n reps + separate @allocated rep; keep every rep so spread is visible.
function benchn(f; n=3)
    f()                              # warm up / compile
    GC.gc(); GC.gc()
    ts = Float64[]
    for _ in 1:n
        GC.gc()
        push!(ts, @elapsed f())
    end
    GC.gc()
    a = @allocated f()
    GC.gc()
    tmin = minimum(ts)
    return (time=tmin, times=ts, alloc=a, spread=(maximum(ts)-tmin)/tmin)
end

println("julia $(VERSION)  nthreads=$(Threads.nthreads())  mode=$MODE")
tsetup = @elapsed elevation = tutorial_elevation(LEVEL)
_, g3 = NG.spatialparts(elevation)
const N = length(elevation)
nnan = count(isnan, parent(elevation))
@printf("level %d  ncells=%d  NaN=%d (%.3f%%)  setup=%.1fs\n", LEVEL, N, nnan, 100nnan/N, tsetup)
const BUDGET = 24 * N        # centroid-table size: the O(ncells) geometry budget
@printf("O(ncells) geometry budget (24 B x ncells) = %.1f MiB\n\n", BUDGET/2^20)

rows = NamedTuple[]
function run!(op, side, f; n=3, note="")
    b = benchn(f; n=n)
    push!(rows, (op=op, side=side, b=b, note=note))
    @printf("%-28s %-34s %10.1f ms / %9.2f MiB   (reps %s, spread %.1f%%)\n",
            op, side, b.time*1e3, b.alloc/2^20,
            join((@sprintf("%.1f", t*1e3) for t in b.times), " "), 100b.spread)
    return b
end

# ============================================================== row 1: TPI ====
# Raw baseline: DGG's own streaming pass, kernel identical to GM's TPI.
function gm_tpi(cell, value, values)
    total = 0.0; count = 0
    for v in values; total += v; count += 1; end
    return Float32(value - total / count)
end
tpi_base = run!("TPI", "raw DGG.mapneighbors(Values())",
    () -> DGG.mapneighbors(gm_tpi, elevation; pass=DGG.Values(), threaded=true); n=5)
tpi_gm = run!("TPI", "GM.topographic_position_index",
    () -> GM.topographic_position_index(elevation); n=5)
tpi_ng = run!("TPI", "NG topographic_position_index",
    () -> NG.topographic_position_index(elevation); n=5)
let a = parent(DGG.mapneighbors(gm_tpi, elevation; pass=DGG.Values(), threaded=true)),
    b = parent(NG.topographic_position_index(elevation)),
    c = parent(GM.topographic_position_index(elevation))
    @printf("  [agree] TPI  NG==baseline: %s   NG==GM: %s\n",
            all(isequal.(a, b)), all(isequal.(c, b)))
end
GC.gc(); GC.gc()

# ================================================== row 2/3: steepest_slope ====
steep_ng = run!("steepest_slope", "NG on-demand geometry",
    () -> NG.steepest_slope(elevation); n=3)
geom = NG.neighborgeometry(g3, NG.NeighborRings(1), NG.STEEPEST_NEEDS)
pre = NG.precompute(geom)
steep_pre = run!("steepest_slope", "NG precompute'd (kernel only)",
    () -> NG.mapneighbors(NG._steepest_slope, elevation, pre); n=3)
steep_build = run!("steepest_slope", "NG precompute build (one-off)",
    () -> NG.precompute(NG.neighborgeometry(g3, NG.NeighborRings(1), NG.STEEPEST_NEEDS)); n=3)
@printf("  [agree] on-demand == precomputed: %s\n",
        all(isequal.(parent(NG.steepest_slope(elevation)),
                     parent(NG.mapneighbors(NG._steepest_slope, elevation, pre)))))
@printf("  [check] NaN in %d -> NaN out %d\n", nnan, count(isnan, parent(NG.steepest_slope(elevation))))
pre = geom = nothing; GC.gc(); GC.gc()

# ================================================= row 4: flowaccumulation ====
fa_gm = run!("flowaccumulation D8", "GM.flowaccumulation(D8())",
    () -> GM.flowaccumulation(elevation; method=GM.D8()); n=3)
fa_ng = run!("flowaccumulation D8", "NG flowaccumulation",
    () -> NG.flowaccumulation(elevation); n=3)

# parity spot stats
let (gacc, gdir) = GM.flowaccumulation(elevation; method=GM.D8()),
    (vacc, vdir) = NG.flowaccumulation(elevation)
    ga, va = Float64.(parent(gacc)), parent(vacc)
    gd, vd = Int.(parent(gdir)), Int.(parent(vdir))
    @printf("  [agree] identical direction codes: %.4f%%   pits GM=%d NG=%d\n",
            100count(gd .== vd)/N, count(==(5), gd), count(==(5), vd))
    x = log1p.(ga); y = log1p.(va)
    @printf("  [agree] cor(log1p acc) = %.6f   max rel diff = %.3e\n",
            cor(x, y), maximum(abs.(ga .- va))/maximum(va))
end
GC.gc(); GC.gc()

# ============================================ rows 5: context (no hard gate) ==
if MODE == "full"
    run!("flow_direction", "NG flow_direction", () -> NG.flow_direction(elevation); n=3)
    run!("slope (PlaneFit)", "NG slope", () -> NG.slope(elevation); n=3)
    GC.gc()
end

# ==================================================================== gates ====
println("\n", "="^96)
@printf("ACCEPT 1  TPI value-only: NG %.1f ms vs raw baseline %.1f ms = %.2fx  (<= 1.50x)   %s\n",
        tpi_ng.time*1e3, tpi_base.time*1e3, tpi_ng.time/tpi_base.time,
        tpi_ng.time/tpi_base.time <= 1.5 ? "PASS" : "FAIL")
@printf("          absolute %.3f s (<= 0.700 s)                                       %s\n",
        tpi_ng.time, tpi_ng.time <= 0.7 ? "PASS" : "FAIL")
@printf("          (context) NG vs GM.TPI %.1f ms = %.2fx\n", tpi_gm.time*1e3, tpi_ng.time/tpi_gm.time)
@printf("ACCEPT 2  steepest_slope on-demand %.3f s (<= 5.000 s)                        %s\n",
        steep_ng.time, steep_ng.time <= 5.0 ? "PASS" : "FAIL")
@printf("          stretch (<= 2.000 s)                                                %s\n",
        steep_ng.time <= 2.0 ? "PASS" : "FAIL")
@printf("          no O(ncells) geometry: %.1f MiB = %.1f B/cell (output floor 8.0 B/cell,\n",
        steep_ng.alloc/2^20, steep_ng.alloc/N)
@printf("          centroid table would be 24.0 B/cell = %.1f MiB); overhead above output\n", BUDGET/2^20)
@printf("          = %.2f B/cell = %.3f x the centroid table                            %s\n",
        (steep_ng.alloc - 8N)/N, (steep_ng.alloc - 8N)/BUDGET,
        (steep_ng.alloc - 8N) < BUDGET/8 ? "PASS" : "FAIL")
@printf("ACCEPT 3  steepest_slope precompute'd %.3f s / %.1f MiB (opt-in; build %.3f s / %.1f MiB)\n",
        steep_pre.time, steep_pre.alloc/2^20, steep_build.time, steep_build.alloc/2^20)
@printf("ACCEPT 4  flowaccumulation NG %.3f s vs GM %.3f s = %.2fx  (NG <= GM)          %s\n",
        fa_ng.time, fa_gm.time, fa_ng.time/fa_gm.time, fa_ng.time <= fa_gm.time ? "PASS" : "FAIL")
@printf("          absolute %.3f s (<= 10.543 s)                                       %s\n",
        fa_ng.time, fa_ng.time <= 10.543 ? "PASS" : "FAIL")
@printf("          memory NG %.1f MiB vs GM %.1f MiB = %.2fx\n",
        fa_ng.alloc/2^20, fa_gm.alloc/2^20, fa_ng.alloc/fa_gm.alloc)
println("="^96)

println("\nMATRIX (level $LEVEL)")
@printf("%-24s %-36s %12s %12s %8s\n", "op", "implementation", "min time", "alloc MiB", "spread")
for r in rows
    @printf("%-24s %-36s %12s %12.2f %7.1f%%\n", r.op, r.side,
            r.b.time < 1 ? @sprintf("%.2f ms", r.b.time*1e3) : @sprintf("%.3f s", r.b.time),
            r.b.alloc/2^20, 100r.b.spread)
end
