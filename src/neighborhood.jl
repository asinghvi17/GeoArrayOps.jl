# Fallbacks visit every index and delegate lookup to `neighbors(dem, cell)`.

"""
    neighbors(dem)

Iterate over `(cell, neighbors(dem, cell))` for every cell in `dem`. Neighbor
iterators contain only indices within `dem`.
"""
neighbors(dem) = ((cell, neighbors(dem, cell)) for cell in eachindex(dem))

_neighborvalues(dem, cell) = (dem[neighbor] for neighbor in neighbors(dem, cell))

"""
    mapneighbors(f, dem; order = nothing, threaded = false)

Apply `f(cell, value, neighbor_values)` to every cell in `dem` and return the
results in a grid shaped like `dem`. `neighbor_values` contains samples from
[`neighbors`](@ref). Calls to `f` must be independent: backends may traverse in
`order` or run calls concurrently when `threaded` is `true`. The fallback
ignores both keywords.
"""
function mapneighbors(f::F, dem; order = nothing, threaded = false) where {F}
    cells = eachindex(dem)
    values = map(cell -> f(cell, dem[cell], _neighborvalues(dem, cell)), cells)
    dst = similar(dem, eltype(values))
    i = 0
    for cell in cells
        dst[cell] = values[i += 1]
    end
    return dst
end
