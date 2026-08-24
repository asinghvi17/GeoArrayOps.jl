# ## Algorithms: the shared helpers
#
# Each algorithm declares the fields it reads. Nothing else is computed for it.
#
# This file holds the three pieces every algorithm file leans on: the nodata
# convention, the unwrap from a `Raster` to its parent array, and the chunking
# the threaded passes share.
#
# Nodata follows Geomorphometry's convention: NaN in, NaN out at that cell. The
# test is on the value, not on the grid, so integer rasters — which have no NaN
# — keep the identical code path with the predicate folded away.

@inline _isnodata(x::AbstractFloat) = isnan(x)
@inline _isnodata(x) = false

_data(A::AbstractArray) = A
_data(r::Raster) = parent(r)

function _chunkranges(n::Int, nchunks::Int=Threads.nthreads())
    n <= 0 && return UnitRange{Int}[]
    width = cld(n, max(1, min(nchunks, n)))
    return [i:min(i + width - 1, n) for i in 1:width:n]
end
