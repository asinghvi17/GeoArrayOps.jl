# Rectilinear regression check for GeomorphometryNG (measurement only).
#   julia --project=examples/GeomorphometryNG.jl/bench accept_rect.jl
# Reproduces the two 1000x1000 Float64 DEMs the PoC's logs used:
#   * smooth DEM  -> steepest_slope / flow_direction / PlaneFit
#   * random DEM  -> flowaccumulation vs Geomorphometry
using Rasters, Printf
import Geomorphometry as GM
using GeomorphometryNG
const NG = GeomorphometryNG

function benchn(f; n=7)
    f(); GC.gc(); GC.gc()
    ts = Float64[]
    for _ in 1:n; GC.gc(); push!(ts, @elapsed f()); end
    GC.gc(); a = @allocated f(); GC.gc()
    tmin = minimum(ts)
    return (time=tmin, times=ts, alloc=a, spread=(maximum(ts)-tmin)/tmin)
end
rows = NamedTuple[]
function run!(name, f; n=7)
    b = benchn(f; n=n); push!(rows, (name=name, b=b))
    @printf("%-44s %9.2f ms / %9.2f MiB   (spread %.1f%%)\n",
            name, b.time*1e3, b.alloc/2^20, 100b.spread)
    return b
end

println("julia $(VERSION)  nthreads=$(Threads.nthreads())")

# ---------------------------------------------- smooth DEM (bench_rect.jl) ----
const NR = 1000
xs = range(0.0, 10.0; length=NR); ys = range(10.0, 0.0; length=NR)
A = [Float64(sinpi(3x)*cospi(2y) * 100 + 7x - 3y) for x in xs, y in ys]
ras = Raster(A, (X(xs), Y(ys)))
println("smooth $(NR)x$(NR) Float64")
s3 = run!("steepest_slope", () -> NG.steepest_slope(ras))
run!("flow_direction", () -> NG.flow_direction(ras))
run!("slope PlaneFit", () -> NG.slope(ras; method=NG.PlaneFit()))
run!("TPI (Value only)", () -> NG.topographic_position_index(ras))

# ------------------------------------ random DEM (bench_rect_sweep.jl) --------
function splitmix01(k::Integer)
    h = UInt64(k) * 0x9E3779B97F4A7C15
    h ⊻= h >> 29; h *= 0xBF58476D1CE4E5B9
    h ⊻= h >> 32; h *= 0x94D049BB133111EB
    h ⊻= h >> 31
    return Float64(h >> 11) * (1 / 9007199254740992.0)
end
z = [200.0 * splitmix01(i + NR*(j-1)) + 60.0*sinpi(2i/NR)*cospi(3j/NR) + 1e-9*(i + NR*(j-1))
     for i in 1:NR, j in 1:NR]
rras = Raster(z, (X(range(0.0; step=2.0, length=NR)), Y(range(5.0NR; step=-5.0, length=NR))))
_, g = NG.spatialparts(rras)
println("random $(NR)x$(NR) Float64 (sweep-family DEM)")
fgm = run!("GM.flowaccumulation(D8)",
           () -> GM.flowaccumulation(z; method=GM.D8(), cellsize=(2.0,-5.0)); n=5)
fng = run!("flowaccumulation(D8)", () -> NG.flowaccumulation(rras, g); n=5)
let (vacc, vdir) = NG.flowaccumulation(rras, g),
    (gacc, gdir) = GM.flowaccumulation(z; method=GM.D8(), cellsize=(2.0,-5.0))
    @printf("  [agree] dirs identical: %s   max rel acc diff: %.3e\n",
            Int.(vdir) == Int.(gdir), maximum(abs.(vacc .- Float64.(gacc)))/maximum(vacc))
end

println("\n", "="^90)
@printf("REGRESSION  steepest_slope %.2f ms (logged ~3.5 ms)\n", s3.time*1e3)
@printf("REGRESSION  flowaccumulation %.1f ms / %.1f MiB  vs GM %.1f ms / %.1f MiB = %.2fx\n",
        fng.time*1e3, fng.alloc/2^20, fgm.time*1e3, fgm.alloc/2^20, fgm.time/fng.time)
@printf("            (logged: NG 135.6 ms / 43.1 MiB, GM 511.3 ms / 120.7 MiB)\n")
println("="^90)
