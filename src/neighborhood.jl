# Generic neighborhood forms over an indexed grid. The fallbacks here run the
# per-cell loop the kernels used to inline; backends with structure to exploit
# (compressed DGGS subsets) override them with amortized traversals.

"""
    neighbors(dem)

Iterator over `(cell, nbrs)` for every cell of `dem`, where `nbrs` is
`neighbors(dem, cell)` — in-domain cells only. Backends may specialize the
one-arg form to amortize per-cell lookup work; the fallback defers to the
two-arg form.
"""
neighbors(dem) = ((cell, neighbors(dem, cell)) for cell in eachindex(dem))

_neighborvalues(dem, cell) = (dem[neighbor] for neighbor in neighbors(dem, cell))

"""
    mapneighbors(f, dem; order = nothing, threaded = false)

Apply `f(cell, value, values)` to every cell of `dem` — `value` the cell's own
sample, `values` an iterator of its neighbors' samples (in-domain only, per the
[`neighbors`](@ref) contract) — and collect the scalar results into a grid
shaped like `dem`. `f` must be pure per cell: backends may thread the traversal
(`threaded`) or reorder it (`order`); the fallback runs a sequential loop and
ignores both.
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
