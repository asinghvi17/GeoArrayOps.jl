# Shared setup: synthetic stand-in for the tutorial's Copernicus GLO-30 tile,
# regridded to IGeo7 exactly as the tutorial does.
import DiscreteGlobalGrids as DGG
using Rasters

"Synthetic 1x1 deg, 1-arcsec (3600x3600) Float32 DEM over the Alps tile N46E010."
function synth_dem(n=3600; lon0=10.0, lat0=46.0)
    xs = range(lon0, lon0 + 1.0; length=n)
    ys = range(lat0 + 1.0, lat0; length=n)   # north-up, descending Y
    A = Array{Float32}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        u = (i - 1) / (n - 1); v = (j - 1) / (n - 1)
        z = 1400.0
        z += 900.0 * sinpi(2u) * cospi(2v)
        z += 450.0 * sinpi(6u + 0.7) * cospi(5v - 0.3)
        z += 210.0 * sinpi(13u - 1.1) * cospi(11v + 0.9)
        z += 95.0  * sinpi(29u + 0.4) * cospi(31v + 1.7)
        z += 40.0  * sinpi(61u) * cospi(59v)
        A[i, j] = Float32(clamp(z, 180.0, 4100.0))
    end
    return Raster(A, (X(xs), Y(ys)); crs=Rasters.EPSG(4326), missingval=Float32(NaN))
end

using Serialization
# Regridded values are cached beside this script. The directory is scratch: the
# caches are reproducible, just slow (minutes at level 13), so they are ignored
# by git and rebuilt on demand.
const CACHE = joinpath(@__DIR__, "cache")
mkpath(CACHE)

# Reproduces the tutorial pipeline: MultiOrderCoverage query + regrid.  The
# regridded VALUES are cached to disk so repeated benchmark processes reuse
# them; the cell axis itself is rebuilt by re-running the query (seconds).
function tutorial_elevation(level::Int; n=3600)
    sys = DGG.IGeo7System()
    ext = Rasters.Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))
    region = DGG.query(sys, DGG.MultiOrderCoverage(ext); level)
    lk = DGG.CellLookup(region)
    f = joinpath(CACHE, "elev_L$(level)_n$(n).jls")
    vals = if isfile(f)
        deserialize(f)::Vector{Float32}
    else
        dem = synth_dem(n)
        igeo7_dem = DGG.regrid(dem; to = region)
        v = Vector{Float32}(collect(parent(igeo7_dem)))
        serialize(f, v)
        v
    end
    @assert length(vals) == length(lk)
    return Raster(vals, (DGG.Cells(lk),); missingval = Float32(NaN))
end

"min-of-n elapsed + post-warmup allocations"
function bench(f; n=5)
    f()                       # warm up / compile
    GC.gc(); GC.gc()
    t = Inf
    for _ in 1:n
        GC.gc()
        t = min(t, @elapsed f())
    end
    GC.gc()
    a = @allocated f()
    GC.gc()
    return (time=t, alloc=a)
end
fmt(b) = string(round(b.time*1e3; digits=2), " ms / ", round(b.alloc/2^20; digits=2), " MiB")
