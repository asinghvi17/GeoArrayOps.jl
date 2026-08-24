# ## Adapters
#
# Matrices use axis 1 for X and axis 2 for Y. Without explicit spacing or a
# manifold, they use unit spacing and planar geometry.

function spatialparts(A::AbstractMatrix; spatialdims=nothing, spacing=nothing, manifold=nothing)
    isnothing(spatialdims) || spatialdims == (1, 2) ||
        throw(ArgumentError("a matrix uses spatialdims=(1, 2)"))
    steps = isnothing(spacing) ? (1.0, 1.0) : spacing
    m = isnothing(manifold) ? Planar() : manifold
    grid = RectilinearGrid(m, (1, 2), axes(A), steps, nothing)
    return A, grid
end

isspatialdim(::Type) = false
isspatialdim(::Type{<:Rasters.XDim}) = true
isspatialdim(::Type{<:Rasters.YDim}) = true
# `isspatialdim(::Type{<:DGG.Cells})` is added by the DiscreteGlobalGrids
# extension, which is the only thing that knows the cell dimension exists.

# The cell backend answers `true` here for a `DGG.CellLookup`. Without
# DiscreteGlobalGrids loaded no lookup is a cell lookup, so a `Cells` dimension
# is not a spatial dimension either and the rectilinear path reports normally.
_iscelllookup(::Any) = false

_astuple(x::Tuple) = x
_astuple(x) = (x,)

# These provisional checks classify a Raster CRS once during grid construction.
# A complete implementation would use the CRS utilities in the Rasters
# extension.
_isgeographiccrs(::Nothing) = false
_isgeographiccrs(crs::EPSG) = GeoFormatTypes.val(crs) == 4326
_isgeographiccrs(crs) = occursin(r"GEOGCS|GEOGCRS"i, string(GeoFormatTypes.val(crs)))

function spatialparts(r::Raster; spatialdims=nothing, spacing=nothing, manifold=nothing)
    found_spatialdims = if isnothing(spatialdims)
        filter(d -> isspatialdim(typeof(d)), Rasters.dims(r))
    else
        Rasters.dims(r, spatialdims)
    end
    selected = _astuple(found_spatialdims)

    if length(selected) == 1 && _iscelllookup(lookup(only(selected)))
        ndims(r) == 1 || throw(ArgumentError("Cells must be the only array axis"))
        isnothing(spacing) || throw(ArgumentError("CellGrid derives geometry from DGG"))
        return r, CellGrid(only(selected); manifold)
    end

    ndims(r) == 2 || throw(ArgumentError("a rectilinear surface must be 2D"))
    length(selected) == 2 ||
        throw(ArgumentError("expected one X and one Y dimension, or one Cells dimension"))

    xdim = only(filter(d -> d isa Rasters.XDim, selected))
    ydim = only(filter(d -> d isa Rasters.YDim, selected))
    P = Rasters.dimnum(r, (Rasters.XDim, Rasters.YDim))
    lookups = (Rasters.lookup(xdim), Rasters.lookup(ydim))
    steps = isnothing(spacing) ? (Float64(step(xdim)), Float64(step(ydim))) : spacing
    rascrs = Rasters.crs(r)
    m = if !isnothing(manifold)
        manifold
    elseif _isgeographiccrs(rascrs)
        Spherical(; radius=AUTHALIC_RADIUS_M)
    else
        Planar()
    end
    grid = RectilinearGrid(m, P, lookups, steps, rascrs)
    return r, grid
end

# ## The façade boundary
#
# This signature defines all grid-construction keywords. Other keywords are
# forwarded to the algorithm method, where Julia reports unsupported keywords
# normally. `spatialdims`, `spacing`, and `manifold` are reserved by the public
# API.
#
# `needs` is deliberately *not* a façade keyword: every named algorithm knows
# its own request. Users writing their own kernels call `mapneighbors` or
# `eachneighbor` and pass `needs` explicitly.

splitspatial(; spatialdims=nothing, spacing=nothing, manifold=nothing, kwargs...) =
    (; spatialdims, spacing, manifold), kwargs

# Algorithms return data without rebuilding input metadata. These methods
# preserve plain-array outputs and restore Raster metadata, including the output
# name and floating-point missing-value convention.

rebuildoutput(input, grid, data; name) = data

# One `rebuild`, not two. An intermediate Raster whose `missingval` does not
# match its new eltype makes Rasters cook one up, and cooking one up is
# `typemin`/`typemax` on the eltype — which `FlowDirection` does not define.
# This is the single place the missingval machinery is reached, which is what
# structurally disarms that crash rather than working around it per call site.
function rebuildoutput(input::Raster, grid, data; name)
    mv = eltype(data) <: AbstractFloat ? eltype(data)(NaN) : nothing
    return data isa Raster ? Rasters.rebuild(data; name, missingval=mv) :
        Rasters.rebuild(input; data, name, missingval=mv)
end

# A reduction has no grid to rebuild onto.
rebuildoutput(input, grid, data::Number, names::Tuple) = data

# An algorithm may return several products; the façade names them positionally.
rebuildoutput(input, grid, data::Tuple, names::Tuple) =
    map((d, nm) -> rebuildoutput(input, grid, d; name=nm), data, names)
rebuildoutput(input, grid, data, names::Tuple) =
    rebuildoutput(input, grid, data; name=first(names))

const OUTPUTNAMES = Dict(
    :flowaccumulation => (:flowaccumulation, :flowdirection),
    :flowdirection => (:flowdirection,),
    :settle => (:settle,),
)
