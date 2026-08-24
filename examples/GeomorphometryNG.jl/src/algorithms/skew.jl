# ## Skewness balancing
#
# The one algorithm in Geomorphometry that is grid-generic for free: it reads no
# neighborhood at all, only the *distribution* of the values. It therefore needs
# nothing from the grid — no topology, no geometry, not even a shape beyond
# "the mask comes back shaped like the input". It is included precisely because
# that is a useful corner of the map: not every terrain algorithm is a
# neighborhood algorithm, and the interface costs this one nothing.

# Geomorphometry gets this from StatsBase; the definition is short enough that
# reproducing it here is cheaper than a dependency, and it keeps the two
# implementations comparable term by term.
function _skewness(v)
    n = length(v)
    n == 0 && return NaN
    mean = sum(v) / n
    m2 = 0.0
    m3 = 0.0
    @inbounds for x in v
        d = x - mean
        m2 += d * d
        m3 += d * d * d
    end
    m2 /= n
    m3 /= n
    return iszero(m2) ? NaN : m3 / m2^1.5
end

"""
    skewness_balancing(dem; mask=nothing)

Ground mask by skewness balancing (Bartels et al., 2006): the threshold below
which the elevation distribution stops being positively skewed. `true` marks a
ground (allowed) cell.

Grid-independent — this reads the value distribution and no neighborhood, so it
is identical on both backends.
"""
function skewness_balancing(ras, grid::AbstractGridSpec)
    values = vec(_data(ras))
    return reshape(_skewnessbalance(values), size(_data(ras)))
end

function _skewnessbalance(input)
    nonfinite = .!isfinite.(input)
    nbad = count(nonfinite)
    values = if nbad > 0
        # Non-finite cells sort to the top and are excluded from the search
        # window, which is what keeps them out of the ground mask.
        v = collect(float.(input))
        v[nonfinite] .= maxintfloat(eltype(v))
        v
    else
        input
    end
    order = sortperm(values)
    sorted = values[order]
    mask = trues(length(input))
    len = length(sorted) - nbad
    len < 2 && return mask

    # Binary search for the largest prefix of the sorted values whose skewness
    # is still non-positive.
    skew = 1.0
    splitby = 2
    step = i = len
    while step >= 1
        skew = _skewness(view(sorted, firstindex(sorted):i))
        step = len ÷ splitby
        splitby <<= 1
        i += skew > 0 ? -step : step
        1 <= i <= len || break
    end
    # `i` may have walked past either end; Geomorphometry lets that stand, and
    # an index past the top simply marks nothing as object, so it stands here
    # too. Only the lower end is clamped, where Geomorphometry would throw.
    (skew <= 0 || i < 1) && (i += 1)
    i = max(i, 1)
    @inbounds for j in i:length(order)
        mask[order[j]] = false
    end
    return mask
end
