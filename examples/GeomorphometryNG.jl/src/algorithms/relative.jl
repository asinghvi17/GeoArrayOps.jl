# ## Local statistics
#
# These match Geomorphometry's definitions (`src/relative.jl`) and read values
# only, so on a cell grid they run on DGG's streaming pass. Their nodata
# behavior is Geomorphometry's: NaN propagates through the arithmetic rather
# than being skipped, so a NaN neighbor makes the whole neighborhood NaN.

const LOCAL_NEEDS = (Value(),)

# TPI: the center minus the mean of its neighbors.
function _tpi_kernel(_, value, neighbors)
    total = 0.0
    count = 0
    for n in neighbors
        total += n.value
        count += 1
    end
    return oftype(float(value), value - total / count)
end

topographic_position_index(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_tpi_kernel, ras, grid, NeighborRings(); needs=LOCAL_NEEDS, kw...)
const TPI = topographic_position_index

# Roughness: the largest absolute difference between the center and a neighbor.
function _roughness_kernel(_, value, neighbors)
    o = zero(value)
    for n in neighbors
        o = max(o, abs(n.value - value))
    end
    return o
end

roughness(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_roughness_kernel, ras, grid, NeighborRings(); needs=LOCAL_NEEDS, kw...)

# TRI: the root of the summed squared differences to the neighbors. `normalize`
# divides by the neighbor count first; `squared=false` gives the mean absolute
# difference instead.
struct RuggednessKernel
    normalize::Bool
    squared::Bool
end

function (k::RuggednessKernel)(_, value, neighbors)
    total = 0.0
    count = 0
    for n in neighbors
        difference = abs(n.value - value)
        total += k.squared ? difference^2 : difference
        count += 1
    end
    k.normalize && (total /= count)
    return oftype(float(value), k.squared ? sqrt(total) : total)
end

function terrain_ruggedness_index(ras, grid::AbstractGridSpec;
        normalize=false, squared=true, kw...)
    if !normalize && !squared
        @warn "TRI: normalize=false and squared=false is not recommended."
    end
    return mapneighbors(RuggednessKernel(normalize, squared), ras, grid, NeighborRings();
        needs=LOCAL_NEEDS, kw...)
end
const TRI = terrain_ruggedness_index

# ## Neighborhood shape statistics
#
# The rest of Geomorphometry's `relative.jl`. Everything below reads records and
# nothing else, so it runs on both backends — with one exception,
# `bathymetric_position_index`, whose *neighborhood* is the thing the vocabulary
# cannot name.

"""
    prominence(dem)

The number of neighbors lower than or equal to the center: 8 is a peak on a
square grid, 0 a pit.
"""
function _prominence_kernel(_, value, neighbors)
    count = 0
    for n in neighbors
        count += n.value <= value
    end
    return Int8(count)
end

prominence(ras, grid::AbstractGridSpec; rings=NeighborRings(), kw...) =
    mapneighbors(_prominence_kernel, ras, grid, rings; needs=LOCAL_NEEDS, kw...)

"""
    percentile_elevation(dem; rings=NeighborRings(1))

The fraction of the neighbors that lie below the center cell. 0 is a local
minimum, 1 a local maximum.
"""
function _percentile_kernel(_, value, neighbors)
    lower = 0
    total = 0
    for n in neighbors
        total += 1
        lower += n.value < value
    end
    return iszero(total) ? NaN : lower / total
end

percentile_elevation(ras, grid::AbstractGridSpec; rings=NeighborRings(), kw...) =
    mapneighbors(_percentile_kernel, ras, grid, rings; needs=LOCAL_NEEDS, kw...)

"""
    pitremoval(dem; limit=0.0)

Raise a cell to its lowest neighbor when *every* neighbor stands more than
`limit` above it. Cells that are not pits keep their own elevation.

Geomorphometry returns `typemax(eltype)` — `Inf` for a float DEM — when every
neighbor is exactly equal to the center, because its accumulator starts at
`typemax` and a tie never clears the pit flag. Here a pit with no strictly
higher neighbor keeps its own value.
"""
struct PitRemovalKernel{T}
    limit::T
end

function (k::PitRemovalKernel)(_, value, neighbors)
    _isnodata(value) && return value
    ispit = true
    lowest = value
    found = false
    for n in neighbors
        n.value == value && continue
        if (n.value - value) > k.limit
            if !found || n.value < lowest
                lowest = n.value
                found = true
            end
        else
            ispit = false
        end
    end
    return (ispit && found) ? lowest : value
end

pitremoval(ras, grid::AbstractGridSpec; limit=zero(eltype(_data(ras))), kw...) =
    mapneighbors(PitRemovalKernel(limit), ras, grid, NeighborRings();
        needs=LOCAL_NEEDS, kw...)

"""
    roughness_index_elevation(dem)

The standard deviation of the residual topography (Cavalli et al., 2008): the
neighborhood spread of "elevation minus local mean". Two record passes, so it
runs on either backend.
"""
function _residual_kernel(_, value, neighbors)
    total = float(value)
    count = 1
    for n in neighbors
        total += n.value
        count += 1
    end
    return float(value) - total / count
end

# The sample standard deviation over the center and its neighbors, matching
# `Statistics.std`'s `n - 1` denominator.
function _neighborhoodstd(_, value, neighbors)
    total = float(value)
    count = 1
    for n in neighbors
        total += n.value
        count += 1
    end
    mean = total / count
    ss = (float(value) - mean)^2
    for n in neighbors
        ss += (n.value - mean)^2
    end
    return count > 1 ? sqrt(ss / (count - 1)) : 0.0
end

function roughness_index_elevation(ras, grid::AbstractGridSpec; kw...)
    residual = mapneighbors(_residual_kernel, ras, grid, NeighborRings();
        needs=LOCAL_NEEDS, kw...)
    return mapneighbors(_neighborhoodstd, residual, grid, NeighborRings();
        needs=LOCAL_NEEDS, kw...)
end

"""
    bathymetric_position_index(dem; inner=2, outer=3)

Not available. BPI is [`topographic_position_index`](@ref) over an *annulus*,
and the neighborhood vocabulary here has no annulus: [`NeighborRings`](@ref) is
cumulative — `NeighborRings(k)` is rings `1..k`, matching
`DGG.neighbors(grid, cell, k)` — so the ring range `inner+1 .. outer` cannot be
requested on either backend. Adding it means a ring-range neighborhood in the
vocabulary and a `k`-th-ring-only query in DiscreteGlobalGrids.
"""
bathymetric_position_index(ras, grid::AbstractGridSpec; inner=2, outer=3, kw...) =
    throw(ArgumentError(
        "bathymetric_position_index needs an annulus neighborhood (rings " *
        "$(inner + 1)..$outer). NeighborRings is cumulative — NeighborRings(k) is " *
        "rings 1..k on both backends — so a ring range has no expression in the " *
        "request vocabulary, and DGG.neighbors is cumulative too. This is a " *
        "vocabulary gap, not a grid limitation."))
const BPI = bathymetric_position_index

# ## Rugosity
#
# The ratio of surface area to planimetric area (Jenness, 2004): fan the ring
# into triangles between angularly adjacent neighbors, and compare the summed
# 3-D triangle areas with the summed projections.
#
# Geomorphometry's version hard-codes the eight corner/edge pairs of a square
# window and divides by `δx * δy`. Taking the *ratio* of the two sums instead
# generalizes: it is exactly Geomorphometry's number on a square grid (the
# projected sum is `δx * δy` up to the shared constant), and it is 1 on flat
# ground for any tessellation, which the fixed constant is not.
#
# The fan needs the neighbors in angular order, which the record stream does not
# promise — so the kernel finds each neighbor's angular successor itself, in
# `O(k²)` with no buffer, which keeps it pure under `Stencils.mapstencil`'s
# unconditional threading.

@inline _crossnorm(a1, a2, a3, b1, b2, b3) =
    hypot(a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)

function _rugosity_kernel(_, value, neighbors)
    _isnodata(value) && return NaN
    surface = 0.0
    planar = 0.0
    for n in neighbors
        _isnodata(n.value) && return NaN
        e1 = n.distance * sind(n.bearing)
        n1 = n.distance * cosd(n.bearing)
        z1 = n.value - value
        # The angular successor: the smallest strictly positive bearing gap.
        gap = 360.0
        e2 = n2 = z2 = 0.0
        found = false
        for o in neighbors
            candidate = mod(o.bearing - n.bearing, 360.0)
            (candidate <= 0 || candidate >= gap) && continue
            gap = candidate
            e2 = o.distance * sind(o.bearing)
            n2 = o.distance * cosd(o.bearing)
            z2 = o.value - value
            found = true
        end
        found || continue
        surface += _crossnorm(e1, n1, z1, e2, n2, z2)
        planar += _crossnorm(e1, n1, 0.0, e2, n2, 0.0)
    end
    return planar > 0 ? surface / planar : NaN
end

"""
    rugosity(dem)

Surface area divided by planimetric area (Jenness, 2004). 1 is flat ground.
"""
rugosity(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_rugosity_kernel, ras, grid, NeighborRings();
        needs=PLANEFIT_NEEDS, kw...)

# ## Entropy
#
# The Shannon entropy of the neighborhood's value histogram, binned to `step`.
#
# Geomorphometry keeps the histogram in two `@MVector`s allocated *once* and
# reused across cells, which decisions §8 already flagged as unportable under
# `Stencils.mapstencil`'s unconditional threading. Those vectors also hold nine
# entries while the window is 5×5, so `Geomorphometry.entropy` throws a
# `BoundsError` on any DEM with more than nine distinct binned values in a 5×5
# window — every real DEM.
#
# The port carries no buffer at all: it counts each distinct value in place, in
# `O(k²)`, which is pure and allocation-free at any ring count.

@inline _roundstep(x, ::Nothing) = float(x)
@inline _roundstep(x, step) = round(float(x) / step) * step

struct EntropyKernel{S}
    step::S
end

function (k::EntropyKernel)(_, value, neighbors)
    _isnodata(value) && return NaN
    v0 = _roundstep(value, k.step)
    samples = 1
    for n in neighbors
        _isnodata(n.value) || (samples += 1)
    end
    # The center's own bin, then every neighbor bin not already counted.
    hits = 1
    for n in neighbors
        _isnodata(n.value) && continue
        _roundstep(n.value, k.step) == v0 && (hits += 1)
    end
    p = hits / samples
    entropy = -p * log(p)
    i = 0
    for n in neighbors
        i += 1
        _isnodata(n.value) && continue
        vi = _roundstep(n.value, k.step)
        vi == v0 && continue
        seen = false
        j = 0
        for other in neighbors
            j += 1
            j >= i && break
            _isnodata(other.value) && continue
            if _roundstep(other.value, k.step) == vi
                seen = true
                break
            end
        end
        seen && continue
        hits = 0
        for other in neighbors
            _isnodata(other.value) && continue
            _roundstep(other.value, k.step) == vi && (hits += 1)
        end
        p = hits / samples
        entropy -= p * log(p)
    end
    return entropy
end

"""
    entropy(dem; step=0.5, rings=NeighborRings(2))

Shannon entropy of the neighborhood's elevation histogram, with values binned to
`step` (`nothing` to skip binning). `rings=NeighborRings(2)` is Geomorphometry's
5×5 window on a rectilinear grid. Nodata neighbors are skipped, not binned.
"""
entropy(ras, grid::AbstractGridSpec; step=0.5, rings=NeighborRings(2), kw...) =
    mapneighbors(EntropyKernel(step), ras, grid, rings; needs=LOCAL_NEEDS, kw...)
