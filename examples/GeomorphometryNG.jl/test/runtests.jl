# The behavior checks from `examples/grid_interface_poc_v3.jl`, as a test suite.
#
# Run it from the repo root with
#
#     julia --project=examples/GeomorphometryNG.jl -e 'using Pkg; Pkg.test()'
#
# DiscreteGlobalGrids is a weakdep of the package and a test dependency, so the
# suite needs the environment `Pkg.test` builds from `[extras]`/`[targets]`;
# `--project=examples/GeomorphometryNG.jl` on its own cannot load it.
#
# Each file below is one section of that PoC's `# ## Behavior checks`, with
# `@assert` replaced by `@test` and the two deliberate warnings asserted with
# `@test_logs` rather than through a hand-rolled collecting logger. Nothing else
# changed: the fixtures, the expressions and their order are the PoC's.
#
# The files are `include`d from inside the testsets, so they evaluate at module
# scope and share their fixtures the way the PoC's linear script does, while the
# `@test`s land in the enclosing testset.

using Test

using Rasters
import Stencils
import DiscreteGlobalGrids as DGG
using GeoFormatTypes: EPSG
using GeometryOpsCore: Planar, Spherical, Geodesic
import Geomorphometry as GM
using Geomorphometry: D8, DInf, FD8, FlowDirection, D8D, LDD, ispit

using GeomorphometryNG
# Rasters exports a `cellarea` of its own, so the one under test has to be named
# explicitly rather than picked up from the `using` above.
using GeomorphometryNG: cellarea
# The checks reach below the public surface on purpose: they are the PoC's
# checks, and several of them are about internal representation choices — which
# geometry payload was built, whether a request can stream, how many slots a
# neighbor table has.
using GeomorphometryNG: AUTHALIC_RADIUS_M, Absent, CellAreas, NoCentroids,
    NoGeometry, OnDemandCentroids, Requested, RowGeometry, StoredCentroids,
    UniformGeometry, NORTH_UP_NEIGHBORS,
    LOCAL_NEEDS, PLANEFIT_NEEDS, STEEPEST_NEEDS, IGEO7_TO_LDD,
    axismap, cellindices, cellkey, geometryat, isboundary, manifold,
    neighborrecord, northupstencil, northupwindow, nslots, partitionrange,
    requestfields, slots, storageposition,
    _accumulatedown!, _accumulateweighted!, _boundaryseeds, _cellencode!, _data,
    _fillcellareas!, _planefit_slope, _steepest_slope, _streamable, _tpi_kernel

# The one place the PoC's collecting logger is still wanted: silencing the
# deliberate no-outlet warning while an allocation gate is measured.
silently(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

@testset verbose = true "GeomorphometryNG" begin
    @testset "grid objects and adapters" begin include("gridobjects.jl") end
    @testset "the request API" begin include("requestapi.jl") end
    @testset "local statistics" begin include("localstats.jl") end
    @testset "nodata semantics" begin include("nodata.jl") end
    @testset "spherical rectilinear grids" begin include("spherical.jl") end
    @testset "cell grids" begin include("cellgrid.jl") end
    @testset "cell geometry is requested, not cached" begin include("cellgeometry.jl") end
    @testset "precompute is a performance choice" begin include("precompute.jl") end

    @testset verbose = true "the sweep family" begin
        @testset "helpers and area indexing" begin include("sweep_setup.jl") end
        @testset "rectilinear fixtures" begin include("sweep_rect.jl") end
        @testset "nodata and closed cells" begin include("sweep_nodata.jl") end
        @testset "degenerate grids, methods, façade" begin include("sweep_edge.jl") end
        @testset "cell grids" begin include("sweep_cells.jl") end
    end

    @testset verbose = true "the ported algorithms" begin
        @testset "terrain derivatives" begin include("terrain_derivatives.jl") end
        @testset "neighborhood statistics" begin include("relative_more.jl") end
        @testset "hydrology" begin include("hydrology_more.jl") end
        @testset "multi-direction flow" begin include("multidirection.jl") end
        @testset "spread and skewness balancing" begin include("spread_skew.jl") end
    end
end
