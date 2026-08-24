# ## Terrain derivatives
#
# The steepest-descent pair and the slope estimators. `steepest_slope` and
# `flow_direction` have no Geomorphometry counterpart; `slope` is
# Geomorphometry's, with `Horn` its rectilinear estimator.

const STEEPEST_NEEDS = (Value(), Distance())
# `flow_direction` reports the steepest downhill bearing, so it reads one more
# field than `steepest_slope` and no index.
const DIRECTION_NEEDS = (Value(), Distance(), Bearing())
const PLANEFIT_NEEDS = (Value(), Distance(), Bearing())

# A nodata center has no defined gradient; nodata neighbors are never selected.
function steepest(value, neighbors)
    _isnodata(value) && return (NaN, NaN)
    best_gradient = 0.0
    best_bearing = NaN
    for neighbor in neighbors
        _isnodata(neighbor.value) && continue
        gradient = (value - neighbor.value) / neighbor.distance
        if gradient > best_gradient
            best_gradient = gradient
            best_bearing = neighbor.bearing
        end
    end
    return best_gradient, best_bearing
end

# `steepest_slope` never reads a bearing, so it requests one field fewer than
# `flow_direction` and the driver skips the bearing arithmetic entirely.
function _steepest_slope(_, value, neighbors)
    _isnodata(value) && return NaN
    best_gradient = 0.0
    for neighbor in neighbors
        _isnodata(neighbor.value) && continue
        best_gradient = max(best_gradient, (value - neighbor.value) / neighbor.distance)
    end
    return atand(best_gradient)
end

function _flow_direction(_, value, neighbors)
    _, bearing = steepest(value, neighbors)
    return bearing # Clockwise from north; NaN indicates a pit, edge outlet or nodata
end

steepest_slope(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_steepest_slope, ras, grid, NeighborRings(); needs=STEEPEST_NEEDS, kw...)
flow_direction(ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_flow_direction, ras, grid, NeighborRings(); needs=DIRECTION_NEEDS, kw...)

# Slope dispatches to an estimator chosen for the grid type. Horn uses a named
# rectilinear window. PlaneFit converts each neighbor's distance and bearing to
# local east and north coordinates, so it works with either grid type.

"""Third-order finite difference over all eight neighbors (Horn, 1981). Windowed."""
struct Horn end
"""Second-order finite difference over the four orthogonal neighbors (Zevenbergen & Thorne, 1987). Windowed."""
struct ZevenbergenThorne end
"""Steepest absolute inter-cell gradient. Reads records, so it runs on any grid."""
struct MaximumDownwardGradient end
const MDG = MaximumDownwardGradient
"""Least-squares tangent plane through the neighbor ring. Reads records, so it runs on any grid."""
struct PlaneFit end

# The two windowed estimators are rectilinear by construction: they read named
# geographic directions out of a square window. The refusal on a cell grid names
# the generic alternative rather than reporting a bare `MethodError`.
_windowedonly(m) = throw(ArgumentError(
    "$(nameof(typeof(m)))() is a windowed estimator and requires a rectilinear grid; " *
    "use PlaneFit() (or MDG() for slope) on a cell grid"))

slope(ras, grid::AbstractGridSpec; method=defaultmethod(slope, grid), kw...) =
    _slope(method, ras, grid; kw...)

defaultmethod(::typeof(slope), ::RectilinearGrid) = Horn()
defaultmethod(::typeof(slope), ::CellGrid) = PlaneFit()

# ### The two windowed gradient estimators
#
# Both return `(∂z/∂east, ∂z/∂north)` in the manifold's units, read by
# geographic name. Nothing here compensates for storage order: `northupwindow`
# absorbed it once.

@inline function _horngradients(w, (sx, sy))
    dzdx = ((w.northeast + 2w.east + w.southeast) -
            (w.northwest + 2w.west + w.southwest)) / (8sx)
    dzdy = ((w.northeast + 2w.north + w.northwest) -
            (w.southeast + 2w.south + w.southwest)) / (8sy)
    return dzdx, dzdy
end

@inline _ztgradients(w, (sx, sy)) = ((w.east - w.west) / (2sx), (w.north - w.south) / (2sy))

@inline _gradients(::Horn, w, sp) = _horngradients(w, sp)
@inline _gradients(::ZevenbergenThorne, w, sp) = _ztgradients(w, sp)

function _horn_slope(w, (sx, sy))
    # Horn reads only the ring, so the center's nodata state must be tested
    # explicitly to keep the NaN-in/NaN-out convention.
    _isnodata(Stencils.center(w)) && return NaN
    dzdx, dzdy = _horngradients(w, (sx, sy))
    return atand(hypot(dzdx, dzdy))
end

# The tangent-plane fit is shared by `slope`, `aspect` and the generic
# `hillshade`: one least-squares pass over the ring, reported as
# `(∂z/∂east, ∂z/∂north)`. It is `(NaN, NaN)` when the fit is degenerate.
@inline function _planefitgradients(value, neighbors)
    see = sen = snn = sze = szn = 0.0
    for n in neighbors
        _isnodata(n.value) && continue
        east = n.distance * sind(n.bearing)
        north = n.distance * cosd(n.bearing)
        dz = n.value - value
        see += east * east; sen += east * north; snn += north * north
        sze += dz * east; szn += dz * north
    end
    det = see * snn - sen * sen
    iszero(det) && return (NaN, NaN)
    return ((sze * snn - szn * sen) / det, (szn * see - sze * sen) / det)
end

function _planefit_slope(_, value, neighbors)
    _isnodata(value) && return NaN
    ddeast, ddnorth = _planefitgradients(value, neighbors)
    return atand(hypot(ddeast, ddnorth))
end

# The steepest *absolute* inter-cell gradient. `steepest_slope` is the downhill
# half of this; MDG counts an uphill neighbor too, which is Geomorphometry's
# non-matrix `slope` default.
function _mdg_slope(_, value, neighbors)
    _isnodata(value) && return NaN
    best = 0.0
    for n in neighbors
        _isnodata(n.value) && continue
        best = max(best, abs(n.value - value) / n.distance)
    end
    return atand(best)
end

function _slope(method::Union{Horn,ZevenbergenThorne}, ras, grid::RectilinearGrid;
        exaggeration=1.0)
    return mapwindow(ras, grid) do w, sp
        _isnodata(Stencils.center(w)) && return NaN
        dzdx, dzdy = _gradients(method, w, sp)
        atand(hypot(dzdx, dzdy) * exaggeration)
    end
end
_slope(method::Union{Horn,ZevenbergenThorne}, ras, ::CellGrid; kw...) = _windowedonly(method)

function _slope(::PlaneFit, ras, grid::AbstractGridSpec; exaggeration=1.0, kw...)
    exaggeration == 1.0 && return mapneighbors(_planefit_slope, ras, grid, NeighborRings();
        needs=PLANEFIT_NEEDS, kw...)
    return mapneighbors(ras, grid, NeighborRings(); needs=PLANEFIT_NEEDS, kw...) do I, v, nbs
        s = _planefit_slope(I, v, nbs)
        isnan(s) ? s : atand(tand(s) * exaggeration)
    end
end

_slope(::MaximumDownwardGradient, ras, grid::AbstractGridSpec; exaggeration=1.0, kw...) =
    exaggeration == 1.0 ?
    mapneighbors(_mdg_slope, ras, grid, NeighborRings(); needs=STEEPEST_NEEDS, kw...) :
    mapneighbors(ras, grid, NeighborRings(); needs=STEEPEST_NEEDS, kw...) do I, v, nbs
        s = _mdg_slope(I, v, nbs)
        isnan(s) ? s : atand(tand(s) * exaggeration)
    end

# ## Aspect
#
# The downslope bearing, degrees clockwise from north. This is the same
# convention `flow_direction` reports, and the same one Geomorphometry's
# `aspect` and `_localaspect` report — but here it is one definition shared by
# the windowed and the record estimators rather than two that happen to agree.

aspect(ras, grid::AbstractGridSpec; method=defaultmethod(aspect, grid), kw...) =
    _aspect(method, ras, grid; kw...)

defaultmethod(::typeof(aspect), ::RectilinearGrid) = Horn()
defaultmethod(::typeof(aspect), ::CellGrid) = PlaneFit()

# A flat cell has no aspect. Geomorphometry returns `atand(0, 0) == 0` there,
# i.e. "north"; NaN is the honest answer and matches this package's `NaN` for a
# pit in `flow_direction`.
@inline _bearingof(dzdx, dzdy) =
    (iszero(dzdx) && iszero(dzdy)) ? NaN : mod(atand(-dzdx, -dzdy), 360.0)

function _aspect(method::Union{Horn,ZevenbergenThorne}, ras, grid::RectilinearGrid)
    return mapwindow(ras, grid) do w, sp
        _isnodata(Stencils.center(w)) && return NaN
        _bearingof(_gradients(method, w, sp)...)
    end
end
_aspect(method::Union{Horn,ZevenbergenThorne}, ras, ::CellGrid) = _windowedonly(method)

# The tangent-plane fit gives east and north gradients directly, so the aspect
# is the same two-argument `atand` the windowed estimators use.
function _planefit_aspect(_, value, neighbors)
    _isnodata(value) && return NaN
    ddeast, ddnorth = _planefitgradients(value, neighbors)
    isnan(ddeast) && return NaN
    return _bearingof(ddeast, ddnorth)
end

_aspect(::PlaneFit, ras, grid::AbstractGridSpec; kw...) =
    mapneighbors(_planefit_aspect, ras, grid, NeighborRings(); needs=PLANEFIT_NEEDS, kw...)


# ## Curvature
#
# The second-derivative family. Every member reads a *named* window, so the
# whole family is rectilinear by design (decisions §6, family 2): a quadratic
# surface fitted to a square window has no hexagonal counterpart, and the
# refusal on a cell grid says so rather than degrading silently.
#
# `radius` dilates the window. Unlike Geomorphometry's `scaled8nb(radius)`, the
# spacing handed to the kernel is dilated with it, so the derivatives stay in
# units of elevation per squared ground distance at any radius.

# The nine window samples, in Geomorphometry's `Z1..Z9` naming, so the
# coefficient formulas below can be read against the papers unchanged.
# Z1 southwest, Z2 south, Z3 southeast, Z4 west, Z5 center, Z6 east,
# Z7 northwest, Z8 north, Z9 northeast.
@inline _z(w) = (w.southwest, w.south, w.southeast, w.west, Stencils.center(w),
    w.east, w.northwest, w.north, w.northeast)

# Zevenbergen and Thorne's pure second derivatives, which is all `laplacian`
# reads of that surface.
@inline function _ztsecond(w, (sx, sy))
    Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9 = _z(w)
    return ((Z4 + Z6) / 2 - Z5) / sx^2, ((Z2 + Z8) / 2 - Z5) / sy^2
end

"""
    laplacian(dem; radius=1, gis=false)

The second derivative of elevation: positive on ridges, negative in valleys.
`gis=true` scales by 100, the convention some GIS tools report. Rectilinear
grids only.
"""
laplacian(ras, grid::AbstractGridSpec; radius=1, gis=false) =
    _laplacian(ras, grid; radius, gis)

function _laplacian(ras, grid::RectilinearGrid; radius=1, gis=false)
    scale = gis ? 100.0 : 1.0
    return mapwindow(ras, grid; radius) do w, sp
        _isnodata(Stencils.center(w)) && return NaN
        d, e = _ztsecond(w, sp)
        -2 * (d + e) * scale
    end
end
_laplacian(ras, ::CellGrid; kw...) = throw(ArgumentError(
    "laplacian fits a quadratic surface to a named square window and requires a " *
    "rectilinear grid; no cell-grid estimator of it is defined here"))

# LandSerf's quadratic `z = Ax² + By² + Cxy + Dx + Ey + F`, the surface all
# three curvatures below are read off. `d` and `e` are the first derivatives in
# east and north; `a`, `b`, `c` the second.
@inline function _landserf(w, (sx, sy))
    Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9 = _z(w)
    a = (Z6 + Z4 - 2Z5) / sx^2
    b = (Z8 + Z2 - 2Z5) / sy^2
    c = (Z9 - Z7 - Z3 + Z1) / (4 * sx * sy)
    d = (Z6 - Z4) / (2sx)
    e = (Z8 - Z2) / (2sy)
    return a, b, c, d, e
end

@inline _profilecurv(a, b, c, d, e) =
    (-2 * (a * d^2 + c * d * e + b * e^2)) / ((d^2 + e^2) * (1 + d^2 + e^2)^1.5)
@inline _tangentialcurv(a, b, c, d, e) =
    -2 * (a * e^2 - c * d * e + b * d^2) / ((d^2 + e^2) * sqrt(1 + d^2 + e^2))
@inline _plancurv(a, b, c, d, e) =
    (2 * (b * d^2 - c * d * e + a * e^2)) / ((1 + d^2 + e^2)^1.5)

for (name, kernel, what) in (
        (:profile_curvature, :_profilecurv, "normal slope line (profile) curvature"),
        (:tangential_curvature, :_tangentialcurv, "normal contour (tangential) curvature"),
        (:plan_curvature, :_plancurv, "projected contour (plan) curvature"))
    inner = Symbol(:_, name)
    @eval begin
        """
            $($(QuoteNode(name)))(dem; radius=1)

        $($what) as defined by Minár et al. (2020). Rectilinear grids only.
        """
        $name(ras, grid::AbstractGridSpec; radius=1) = $inner(ras, grid; radius)

        $inner(ras, grid::RectilinearGrid; radius=1) =
            mapwindow(ras, grid; radius) do w, sp
                _isnodata(Stencils.center(w)) && return NaN
                $kernel(_landserf(w, sp)...)
            end

        $inner(ras, ::CellGrid; kw...) = throw(ArgumentError(
            string($(QuoteNode(name)),
                " fits a quadratic surface to a named square window and requires a ",
                "rectilinear grid; no cell-grid estimator of it is defined here")))
    end
end

# ## Illumination
#
# `hillshade` needs a gradient and nothing else, so unlike the curvatures it has
# a generic estimator: the same `PlaneFit` that serves `slope` and `aspect`.
# Geomorphometry returns `Union{Missing,UInt8}`; this returns `Float64` in
# `0:255` so that the NaN-in/NaN-out convention has somewhere to land, and so
# that the border a windowed driver cannot answer is NaN rather than `missing`.

"""
    hillshade(dem; azimuth=315.0, zenith=45.0, method=defaultmethod(hillshade, grid))

Simulated illumination from a light source at `azimuth` and `zenith` degrees.
Returns `Float64` illumination in `0:255`, `NaN` where undefined.
"""
hillshade(ras, grid::AbstractGridSpec; azimuth=315.0, zenith=45.0,
    method=defaultmethod(hillshade, grid), kw...) =
    _hillshade(method, ras, grid; azimuth, zenith, kw...)

"""
    multihillshade(dem; azimuth=[225, 270, 315, 360], zenith=45.0)

[`hillshade`](@ref) combined over several light sources, weighted by the square
of the sine of the aspect–azimuth difference (Mark, 1992).
"""
multihillshade(ras, grid::AbstractGridSpec; azimuth=(225.0, 270.0, 315.0, 360.0),
    zenith=45.0, method=defaultmethod(hillshade, grid), kw...) =
    _multihillshade(method, ras, grid; azimuth, zenith, kw...)

defaultmethod(::typeof(hillshade), ::RectilinearGrid) = Horn()
defaultmethod(::typeof(hillshade), ::CellGrid) = PlaneFit()

# Geomorphometry's illumination formula verbatim, taking the gradients as
# arguments so the windowed and record estimators share it. `a` is the aspect as
# an *unsigned* angle in radians measured the way that code measures it; it is
# kept in that form so the two implementations remain comparable term by term.
@inline function _illumangle(dzdx, dzdy)
    if dzdx != 0
        a = atan(dzdx, dzdy)
        a < 0 && (a += 2π)
    else
        a = π / 2
        dzdy < 0 && (a += 2π)
    end
    return a
end

# The angle the illumination formula wants is the *downslope* bearing, so the
# gradients are negated on the way in. Geomorphometry gets there by accident:
# `LocalFilters` addresses its kernel as `A[i - j]`, so the sums its `horn`
# helper labels "north" and "east" are actually the south and west ones, and its
# `δzδx`/`δzδy` are the negated true gradients. `aspect` then undoes the flip in
# `compass`, and `hillshade` — which has no `compass` — relies on it.
@inline function _shade(dzdx, dzdy, zenithr, azimuthr)
    a = _illumangle(-dzdx, -dzdy)
    s = atan(hypot(dzdx, dzdy))
    return max(0.0, 255 * (cos(zenithr) * cos(s) + sin(zenithr) * sin(s) * cos(azimuthr - a)))
end

@inline function _multishade(dzdx, dzdy, zenithr, azimuth)
    a = _illumangle(-dzdx, -dzdy)
    s = atan(hypot(dzdx, dzdy))
    α = cos(zenithr) * cos(s)
    β = sin(zenithr) * sin(s)
    total = 0.0
    weights = 0.0
    for az in azimuth
        weight = sin(a - deg2rad(az - 90))^2
        weights += weight
        total += weight * (α + β * cos(deg2rad(az) - a))
    end
    return max(0.0, 255 * total / weights)
end

function _hillshade(method::Union{Horn,ZevenbergenThorne}, ras, grid::RectilinearGrid;
        azimuth, zenith)
    zenithr, azimuthr = deg2rad(zenith), deg2rad(azimuth)
    return mapwindow(ras, grid) do w, sp
        _isnodata(Stencils.center(w)) && return NaN
        _shade(_gradients(method, w, sp)..., zenithr, azimuthr)
    end
end
_hillshade(method::Union{Horn,ZevenbergenThorne}, ras, ::CellGrid; kw...) =
    _windowedonly(method)

function _hillshade(::PlaneFit, ras, grid::AbstractGridSpec; azimuth, zenith, kw...)
    zenithr, azimuthr = deg2rad(zenith), deg2rad(azimuth)
    return mapneighbors(ras, grid, NeighborRings(); needs=PLANEFIT_NEEDS, kw...) do _, v, nbs
        _isnodata(v) && return NaN
        dzdx, dzdy = _planefitgradients(v, nbs)
        isnan(dzdx) ? NaN : _shade(dzdx, dzdy, zenithr, azimuthr)
    end
end

function _multihillshade(method::Union{Horn,ZevenbergenThorne}, ras, grid::RectilinearGrid;
        azimuth, zenith)
    zenithr = deg2rad(zenith)
    az = Tuple(azimuth)
    return mapwindow(ras, grid) do w, sp
        _isnodata(Stencils.center(w)) && return NaN
        _multishade(_gradients(method, w, sp)..., zenithr, az)
    end
end
_multihillshade(method::Union{Horn,ZevenbergenThorne}, ras, ::CellGrid; kw...) =
    _windowedonly(method)

function _multihillshade(::PlaneFit, ras, grid::AbstractGridSpec; azimuth, zenith, kw...)
    zenithr = deg2rad(zenith)
    az = Tuple(azimuth)
    return mapneighbors(ras, grid, NeighborRings(); needs=PLANEFIT_NEEDS, kw...) do _, v, nbs
        _isnodata(v) && return NaN
        dzdx, dzdy = _planefitgradients(v, nbs)
        isnan(dzdx) ? NaN : _multishade(dzdx, dzdy, zenithr, az)
    end
end

"""
    pssm(dem; exaggeration=2.3, method=defaultmethod(slope, grid))

Perceptually Shaded Slope Map (Pingel and Clarke, 2014): exaggerated
[`slope`](@ref), suitable for display as a grayscale image.
"""
pssm(ras, grid::AbstractGridSpec; exaggeration=2.3,
    method=defaultmethod(slope, grid), kw...) =
    _slope(method, ras, grid; exaggeration, kw...)
