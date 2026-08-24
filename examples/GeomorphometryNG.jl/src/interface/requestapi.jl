# ## The request: which per-neighbor quantities does this kernel read?
#
# Every per-neighbor quantity is named by a singleton type under a shared
# supertype, and a call requests a tuple of them. The tuple's *type* carries
# the request, so it specializes through `Stencils.mapstencil` and
# `DGG.mapneighbors` — both of which are function barriers — and the driver
# never branches on the request at run time.

abstract type NeighborField end

"""Position of the neighbor: a `CartesianIndex` or a positioned DGG cell."""
struct Index <: NeighborField end
"""The neighbor's value in the array being mapped."""
struct Value <: NeighborField end
"""Distance from the center cell to the neighbor, in the manifold's units."""
struct Distance <: NeighborField end
"""Tangent-plane bearing to the neighbor, degrees clockwise from north."""
struct Bearing <: NeighborField end

# Records always list their fields in this order, whatever order the request
# used, so a kernel written against `(Value(), Distance())` also accepts records
# built from `(Distance(), Value())`.
const CANONICAL_FIELDS = (Index(), Value(), Distance(), Bearing())

# `Index` is deliberately *not* in the default: it is the field that forces the
# `CellGrid` driver off DGG's streaming value pass.
const DEFAULT_NEEDS = (Value(),)

# Request membership is answered with singleton types rather than `Bool`s, so
# every branch below is a dispatch that resolves at compile time even when
# constant propagation does not fire.
struct Requested end
struct Absent end

@inline _asked(::Type{F}, ::Tuple{}) where {F<:NeighborField} = Absent()
@inline _asked(::Type{F}, needs::Tuple) where {F<:NeighborField} =
    _askedhead(F, first(needs), Base.tail(needs))
@inline _askedhead(::Type{F}, ::F, rest::Tuple) where {F<:NeighborField} = Requested()
@inline _askedhead(::Type{F}, ::NeighborField, rest::Tuple) where {F<:NeighborField} =
    _asked(F, rest)

@inline _either(::Absent, ::Absent) = Absent()
@inline _either(::Requested, ::Any) = Requested()
@inline _either(::Absent, ::Requested) = Requested()

# `Distance` and `Bearing` share one geometry provider, so they are gated
# together.
@inline wantsgeometry(fields::Tuple) =
    _either(_asked(Distance, fields), _asked(Bearing, fields))

@inline _keepfield(::Requested, f, rest::Tuple) = (f, rest...)
@inline _keepfield(::Absent, f, rest::Tuple) = rest

@inline _canonical(needs::Tuple, ::Tuple{}) = ()
@inline _canonical(needs::Tuple, fields::Tuple) =
    _keepfield(_asked(typeof(first(fields)), needs), first(fields),
        _canonical(needs, Base.tail(fields)))

function requestfields(needs)
    needs isa Tuple || throw(ArgumentError(
        "needs must be a tuple of NeighborField singletons, got $(typeof(needs))"))
    all(n -> n isa NeighborField, needs) || throw(ArgumentError(
        "needs must be a tuple of NeighborField singletons, " *
        "e.g. (Value(), Distance()); got $needs"))
    return _canonical(needs, CANONICAL_FIELDS)
end

# One record field per requested quantity. `geo` carries the pair the geometry
# provider produced, or `nothing` when no geometry was requested — in which case
# the `Distance`/`Bearing` methods are unreachable by construction.
@inline _recordfield(::Index, index, value, geo) = (; index)
@inline _recordfield(::Value, index, value, geo) = (; value)
@inline _recordfield(::Distance, index, value, geo) = (; distance=geo.distance)
@inline _recordfield(::Bearing, index, value, geo) = (; bearing=geo.bearing)

@inline neighborrecord(::Tuple{}, index, value, geo) = NamedTuple()
@inline neighborrecord(fields::Tuple, index, value, geo) =
    merge(_recordfield(first(fields), index, value, geo),
        neighborrecord(Base.tail(fields), index, value, geo))
