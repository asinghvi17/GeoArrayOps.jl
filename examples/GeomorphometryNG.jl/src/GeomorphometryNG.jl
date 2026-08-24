"""
    GeomorphometryNG

A grid-interface experiment for Geomorphometry: one set of terrain algorithms
that runs unchanged over rectilinear rasters and DGG cell grids.

A grid spec ([`RectilinearGrid`](@ref), [`CellGrid`](@ref)) carries the topology
and the manifold; a *request* — a tuple of [`NeighborField`](@ref) singletons
such as `(Value(), Distance())` — states which per-neighbor quantities a kernel
reads, and the drivers ([`mapneighbors`](@ref), [`eachneighbor`](@ref)) build
exactly those and nothing else. [`neighbortable`](@ref) is the traversal
analogue, serving the sweep family ([`floodsweep`](@ref), [`settle`](@ref),
[`flowdirection`](@ref), [`flowaccumulation`](@ref)) with a slot-stable,
position-keyed adjacency.

This is the packaged form of `examples/grid_interface_poc_v3.jl`; that file
remains the historical record and the narrative introduction.

The cell-grid backend lives in the `DiscreteGlobalGrids` extension: load
`DiscreteGlobalGrids` to get it.
"""
module GeomorphometryNG

using Rasters
import Stencils
import GeoFormatTypes
using GeoFormatTypes: EPSG
using GeometryOpsCore: Manifold, Planar, Spherical, Geodesic
# The priority flood uses QuickHeaps' array-backed queue directly rather than
# through Geomorphometry's re-export.
using QuickHeaps: FastPriorityQueue, enqueue!, dequeue!
# The sweep family reuses Geomorphometry's direction encodings and its
# flow-method singletons. The import is qualified because both namespaces define
# `slope`, `roughness`, `flowaccumulation`, `Horn`, ....
import Geomorphometry as GM
using Geomorphometry: D8, DInf, FD8, FlowDirection, FlowDirectionConvention,
    FlowDirectionMethod, D8D, LDD, ispit

# Grid specs and the façade adapter that builds them.
export AbstractGridSpec, RectilinearGrid, CellGrid, spatialparts
# The neighborhood vocabulary and the request API.
export NeighborRings, NeighborField, Index, Value, Distance, Bearing
# Reusable geometry and the drivers a user kernel is written against.
export neighborgeometry, precompute, mapneighbors, mapwindow, eachneighbor
# Areas and the traversal primitive.
export cellarea, neighbortable
# Named algorithms.
export steepest_slope, flow_direction, slope, aspect, defaultmethod
export Horn, PlaneFit, ZevenbergenThorne, MaximumDownwardGradient, MDG
export laplacian, plan_curvature, profile_curvature, tangential_curvature
export hillshade, multihillshade, pssm
export topographic_position_index, TPI, terrain_ruggedness_index, TRI, roughness
export bathymetric_position_index, BPI, roughness_index_elevation, rugosity, entropy
export prominence, percentile_elevation, pitremoval
# The sweep family, and what it carries.
export floodsweep, settle, flowdirection, flowaccumulation, flowpartition,
    RingSlot, RingMask, downstreamposition, slotgeometry
export depression_depth, depression_volume, drainage_potential
export topographic_wetness_index, TWI, stream_power_index, SPI
export height_above_nearest_drainage
export skewness_balancing
export spread, Tomlin

# The grid interface: the specs, the request API, the reusable geometry, the
# drivers a kernel is written against, the two hoisted structures
# (`neighborgeometry`, `neighbortable`) and the façade adapters. Dependency
# order.
include("interface/gridspecs.jl")
include("interface/requestapi.jl")
include("interface/geometry.jl")
include("interface/drivers.jl")
include("interface/cellarea.jl")
include("interface/neighbortable.jl")
include("interface/facade.jl")

# The algorithms, filed the way Geomorphometry files them.
include("algorithms/utils.jl")
include("algorithms/terrain.jl")
include("algorithms/relative.jl")
include("algorithms/hydrology.jl")
include("algorithms/flowdir.jl")
include("algorithms/skew.jl")
include("algorithms/spread.jl")

# The façade methods are generated last, over every algorithm defined above.
for f in (:steepest_slope, :flow_direction, :slope, :aspect, :settle,
          :flowdirection, :flowaccumulation, :topographic_position_index,
          :terrain_ruggedness_index, :roughness, :laplacian, :plan_curvature,
          :profile_curvature, :tangential_curvature, :hillshade, :multihillshade,
          :pssm, :bathymetric_position_index, :roughness_index_elevation,
          :rugosity, :entropy, :prominence, :percentile_elevation, :pitremoval,
          :depression_depth, :depression_volume, :drainage_potential,
          :topographic_wetness_index, :stream_power_index,
          :height_above_nearest_drainage, :skewness_balancing)
    names = get(OUTPUTNAMES, f, (f,))
    @eval function $f(input; kwargs...)
        spatial, rest = splitspatial(; kwargs...)
        ras, grid = spatialparts(input; spatial...)
        data = $f(ras, grid; rest...)
        return rebuildoutput(input, grid, data, $(QuoteNode(names)))
    end
end

# Two public entry points the loop above cannot generate. `flowpartition` takes
# its method first, and `spread` takes three arrays of which only `friction`
# carries the grid.

function flowpartition(method::FlowDirectionMethod, input; kwargs...)
    spatial, rest = splitspatial(; kwargs...)
    ras, grid = spatialparts(input; spatial...)
    return flowpartition(method, ras, grid; rest...)
end

function spread(points, initial, friction; kwargs...)
    spatial, rest = splitspatial(; kwargs...)
    ras, grid = spatialparts(friction; spatial...)
    data = spread(points, initial, ras, grid; rest...)
    return rebuildoutput(friction, grid, data; name=:spread)
end

end # module GeomorphometryNG
