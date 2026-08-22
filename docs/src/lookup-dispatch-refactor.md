# Grid-driven spatial dispatch in Geomorphometry

Status: design discussion; no package implementation has been started. A
standalone executable proof of concept is available in
[`examples/grid_interface_poc.jl`](../../examples/grid_interface_poc.jl).

## Summary

Geomorphometry currently supports ordinary matrices, Rasters.jl rasters,
GeoArrays, and—in the DiscreteGlobalGrids extension—rasters indexed by discrete
global grid cells. These inputs expose their spatial structure in different
ways. Most of the core package nevertheless operates as though the first array
axis is X/east and the second is Y/north.

The proposed refactor introduces a common boundary between public APIs and
numerical algorithms. At that boundary, an input is resolved into two parts:

- `ras`: the array or wrapper on which the calculation operates;
- `grid`: a Geomorphometry-owned object describing the topology, geometry, and
  relationship between logical spatial axes and storage axes.

Internal methods dispatch on the grid. This makes grid capabilities explicit
instead of using dimensionality, concrete Raster types, or DimensionalData
lookup types as indirect proxies.

The adaptation layer may retain DimensionalData dimensions and lookups inside
the grid, but numerical algorithms do not inspect them. They speak a smaller
vocabulary: cells, neighborhoods, stencils, spacing, area, distance, bearing,
direction, storage access, and output allocation.

Rectilinear axis order is represented directly in type space. For example,
`RectilinearGrid{(2, 1)}` means that logical X is stored on array axis 2 and
logical Y is stored on array axis 1. A `CellGrid` has no axis-order parameter:
for now, a Cells dimension must be the only array axis.

## Goals

The refactor should:

- provide a common spatial-input protocol for Raster and `AbstractMatrix`;
- make grid capabilities available for dispatch without exposing lookup details
  to algorithms;
- accept `spatialdims=nothing` by default and autodetect known spatial
  dimensions;
- allow `spatialdims=...` to select spatial dimensions explicitly;
- handle Raster dimensions stored as either X,Y or Y,X;
- allow extensions to declare additional spatial dimension types;
- recognize `DGG.Cells` as spatial in the DiscreteGlobalGrids extension;
- support both two-dimensional rectilinear surfaces and one-dimensional DGG
  cell collections;
- preserve input wrappers, metadata, dimension order, and output lookup
  information;
- make unsupported grid/algorithm combinations fail clearly;
- preserve the current behavior of ordinary matrices unless an existing
  inconsistency is being deliberately corrected.

## Non-goals for the first stage

The first stage does not need to make every algorithm work on every grid.
Lookup-driven dispatch should express capabilities rather than imply universal
support. For example, a Horn finite-difference stencil is meaningful on a
regular rectilinear grid but has no automatic interpretation on an arbitrary
DGG cell topology.

The initial refactor does not need to solve general batching across Time, Band,
or other nonspatial dimensions. Rectilinear inputs are two-dimensional
surfaces, and a DGG Cells dimension must be the sole axis. Support for spatial
axes embedded in higher-dimensional arrays can be designed later.

## Current behavior and its limitations

### Raster support is incidental

There is currently no central Raster adapter. A two-dimensional Raster works
mainly because it is an `AbstractMatrix`. Core matrix methods receive the Raster
wrapper directly, while the Rasters extension supplies cell size and a special
three-dimensional allocator for horizon outputs.

The relevant glue is in
[`GeomorphometryRastersExt.jl`](../../ext/GeomorphometryRastersExt.jl). It
currently provides:

- Raster-specific horizon allocation;
- Raster-specific cell-size derivation;
- forwarding of rectangular outlets and neighbors to the parent matrix.

This arrangement preserves wrappers in many cases because `similar`, `copy`,
and broadcasting happen to preserve them. It does not provide an explicit place
to resolve spatial dimensions, validate lookups, align multiple inputs, or
rebuild outputs consistently.

### X,Y and Y,X are not equivalent today

The flow-direction convention explicitly defines array axis 1 as X/east and
array axis 2 as Y/north. Terrain stencils use the same convention. See
[`flowdir.jl`](../../src/flowdir.jl) and
[`terrain.jl`](../../src/terrain.jl).

Raster cell size, however, is obtained by querying dimensions in X,Y order:

- DimensionalData returns queried dimensions in the requested query order;
- this query order is independent of the Raster's storage order;
- the numerical data is not permuted to match the query.

For a Raster stored as Y,X, the result is therefore an X,Y cell-size tuple
applied to data whose physical axes are Y,X. With unequal spacing this gives
wrong distances and derivative scaling. Even with equal spacing, directional
results such as aspect and flow-direction codes can be wrong because the axes
have different semantic roles.

Current DimensionalData documentation confirms that `dims` supports selection
by type, instance, symbol, tuple, and predicate, and that a selector tuple can
control result order. `dimnum` is the appropriate way to recover physical array
positions. See the
[DimensionalData information API](https://github.com/rafaqz/DimensionalData.jl/blob/main/docs/src/get_info.md).

### Horizon uses a different axis convention

The horizon implementation describes cell size as row size followed by column
size. It moves north/south along axis 1 and east/west along axis 2. That is the
ordinary matrix-image convention, but it conflicts with the X-first convention
used by terrain derivatives and flow directions.

The existing Raster horizon test uses square cells and equal X/Y steps, so it
checks wrapper preservation without exposing the directional mismatch. See
[`horizon.jl`](../../src/horizon.jl) and its tests in
[`test/horizon.jl`](../../test/horizon.jl).

### DGG support is tied to the whole Raster type

The current DiscreteGlobalGrids extension defines a `CellsRaster` alias for a
one-dimensional Raster whose only dimension is `DGG.Cells` carrying a
`DGG.CellLookup`. It then adds cell-ID-based `eachindex`, `getindex`,
`setindex!`, and `checkbounds` behavior to that Raster type.

Neighborhoods, flow accumulation, HAND, cell area, distance, and bearing are
implemented through methods on the complete Raster type. See
[`GeomorphometryDiscreteGlobalGridsExt.jl`](../../ext/GeomorphometryDiscreteGlobalGridsExt.jl).

This works, but it couples topology and algorithms to the complete nested Raster
type. A `CellGrid` adapter can retain the `CellLookup` privately while exposing
only Geomorphometry's topology and geometry vocabulary. Although DGG itself can
support a Cells axis among additional dimensions, this refactor deliberately
keeps the current one-axis restriction.

## Proposed architecture

### 1. A public normalization boundary

Every public operation that semantically acts on a spatial surface should
resolve its input once. Resolution performs these tasks:

1. discover or resolve the selected spatial dimensions;
2. validate dimension roles and lookup capabilities;
3. construct a Geomorphometry grid that encapsulates the dimensions, lookups,
   topology, geometry, and storage-axis mapping;
4. produce `(ras, grid)` for an internal function barrier;
5. retain enough boundary context to rebuild results correctly.

Composite operations should not call other public operations with the original
input. They should resolve once and call private grid-aware implementations.
This applies to combinations such as PSSM calling slope, or TWI combining slope
and flow accumulation.

The new resolver should not be named `decompose`, because that name is already
exported for flow-direction decomposition.

### 2. A dimension-type oracle

The package should expose or at least document a small extension point that
answers whether a dimension type is spatial.

The intended registrations are:

- the core fallback reports nonspatial;
- the Rasters/DimensionalData integration reports spatial for `XDim` and
  `YDim`;
- the DGG extension reports spatial for `DGG.Cells`.

The oracle should operate on the dimension type, with an instance convenience
method if useful. Type-based classification is stable, easy for extensions to
add, and amenable to specialization.

The oracle belongs to the adapter layer. Numerical algorithms do not call it and
do not depend on DimensionalData dimension types.

### 3. Spatial-dimension resolution

When `spatialdims=nothing`, the resolver should inspect the input dimensions and
select those whose types are classified as spatial.

When `spatialdims` is supplied, it should resolve DimensionalData selectors to
actual dimensions. Types, dimension instances, and symbols are natural forms
because DimensionalData already supports them.

Explicit selection should override autodetection, but it should not silently
confuse query order with storage order. The resolver must use actual dimension
positions to construct the correct grid type.

Central validation should report useful errors for:

- no detected spatial dimensions;
- a selected dimension that does not exist;
- duplicate X or Y roles;
- an unsupported number or combination of spatial dimensions;
- a spatial dimension whose lookup lacks the capability required by the
  selected algorithm;
- unsupported nonspatial dimensions.

### 4. Geomorphometry grid types

The grid is an owned semantic adapter, not merely a dispatch tag. It wraps any
external dimensions and lookups required to implement the grid interface, but
those objects remain private implementation details.

#### Rectilinear grids

The type parameter of `RectilinearGrid` is a tuple mapping logical coordinate
roles to storage axes:

- `RectilinearGrid{(1, 2)}` means X is storage axis 1 and Y is storage axis 2;
- `RectilinearGrid{(2, 1)}` means X is storage axis 2 and Y is storage axis 1.

The parameter always has the meaning `(xaxis, yaxis)`. It is not the storage
order of a dimension tuple and does not include lookup direction.

The grid value may privately hold:

- the resolved X and Y dimensions;
- their underlying lookups;
- coordinate orientation;
- regular spacing or coordinate access;
- CRS or metric context where needed;
- information used to validate and rebuild aligned outputs.

Only facts that change compiled behavior need to appear in type space. The
axis-order tuple is one such fact. Actual coordinate and dimension values remain
ordinary fields. Forward/reverse lookup order may already specialize through
field types; it does not need to be added to the main tuple parameter.

Keeping data in storage order avoids implicit permutation of disk-backed arrays,
GPU arrays, Cartesian indices, masks, and auxiliary inputs. The grid translates
between logical X/Y directions and physical array axes.

#### Cell grids

`CellGrid` represents a grid whose only array axis is a Cells dimension. It does
not carry an axis-order parameter. Its adapter validates the one-dimensional
restriction and privately retains the cell dimension and compatible
`CellLookup`.

The initial design intentionally does not support Time × Cells or another
placement of the Cells axis. This keeps cell storage addressing and traversal
unambiguous while the core protocol is established.

### 5. The algorithm vocabulary

Numerical algorithms should not inspect dimensions or lookups. They should use
a compact grid interface built around five concepts.

#### Cell traversal

- iterate storage-addressable cells;
- obtain boundary or outlet cells;
- choose a suitable traversal order where relevant.

#### Neighborhood topology

- obtain neighboring cells;
- map a function over neighborhoods;
- interpret a neighborhood request such as adjacency, Moore radius, an annulus,
  or a named directional stencil.

#### Geometry

- constant spacing where available;
- cell area;
- distance and bearing between cells;
- logical direction and flow-direction encoding.

#### Storage mapping

- read and write a cell without requiring logical cell identity to equal its
  raw array index;
- translate logical X/Y offsets into storage-axis offsets;
- translate topology-native neighbors into storage positions.

#### Output construction

- allocate a spatially compatible result;
- rebuild the public wrapper with its original order and metadata;
- add output dimensions such as horizon direction.

Rectilinear implementations derive these operations from their privately held
X/Y dimensions and lookups. Cell implementations use cell adjacency, polygon
area, centroids, and relative-cell arithmetic. Algorithms dispatch on the grid
type and use only this vocabulary.

### 6. `mapneighbors` and semantic stencils

Geomorphometry's existing `mapneighbors` should become a central execution
interface over `(ras, grid, neighborhood)`. It is currently a convenience
fallback built on `eachindex(dem)` and `neighbors(dem, cell)`, with a DGG method
on the complete Raster type. The refactor moves that backend decision into the
grid.

For `RectilinearGrid`, a stencil is expressed in logical X/Y coordinates. The
grid transforms semantic directions such as east, west, north, and south into
storage offsets using its `(xaxis, yaxis)` parameter and coordinate orientation.
The same Horn or Moore stencil therefore has the same geographic meaning for
X,Y and Y,X storage.

For `CellGrid`, the initial neighborhood is adjacency. Mapping can use the
lookup's preferred storage order and optimized DGG traversal without changing
Base indexing behavior on Raster.

This supports two useful levels of algorithm:

- topology-only reductions such as roughness, TPI, TRI, and prominence can use
  adjacent-neighbor mapping on either grid;
- directional finite-difference algorithms request a rectilinear stencil and
  dispatch only on `RectilinearGrid`.

`mapneighbors` is not the only primitive. Stateful algorithms such as priority
flood and flow accumulation need direct cell and neighbor iteration because
their traversal depends on previously computed state. They should use the same
grid vocabulary at a lower level.

Iteration policy is distinct from axis order. Axis order belongs to the grid
type. Policies such as storage order, reverse order, threaded mapping, or
priority order are arguments or separate policy tokens supported by the grid.

### 7. AbstractMatrix behavior

An ordinary `AbstractMatrix` has no named dimensions. Its fallback should use
the existing Geomorphometry convention:

- axis 1 has the synthetic X/east role;
- axis 2 has the synthetic Y/north role;
- its grid retains `axes(A)` as the index-coordinate source;
- default cell spacing is `(1, 1)`.

Using `axes(A)` is important. Constructing new one-based ranges from `size(A)`
would discard offset or custom axes and further entrench dense, one-based array
assumptions.

If explicit `spatialdims` is supported for plain matrices, integer axis
selectors are the most natural meaning. Silently ignoring the keyword would be
misleading.

### 8. DGG behavior

A DGG input has one spatial dimension and produces a `CellGrid`. The adapter
requires that Cells be the sole axis and that it contain a compatible
`CellLookup`.

The dimension oracle should recognize `DGG.Cells`, while algorithm dispatch
must separately require a compatible lookup. A `Cells` dimension containing a
generic categorical or no-lookup value is spatial by name but does not provide
DGG neighborhood geometry.

This distinction allows:

- topology-based operations to work on DGG grids;
- rectilinear stencil algorithms to fail clearly;
- DGG-specific optimized traversal to remain available;
- removal of much of the current full-Raster coupling and Base indexing
  overrides once internal algorithms consume `(ras, grid)` directly.

Support still depends on DGG capabilities. For example, current relative-cell
implementations of DInf, FD8, and HAND are available only for IGeo7, while D8
can be represented more generally.

## Cell size, orientation, and geographic meaning

The public meaning of an explicit `cellsize` tuple should be consistent across
the package. The recommended contract is semantic `(x, y)`, independent of
storage order. This agrees with the current `cellsize` documentation and with
the geographic meaning of flow-direction codes.

The rectilinear grid maps semantic X/Y spacing to physical array axes. Lookup
ordering and coordinate direction must be handled centrally:

- forward and reverse ordered lookups affect the sign of coordinate motion;
- distances use absolute physical lengths;
- bearings, aspects, and direction codes retain orientation;
- the same orientation rules must be used by terrain, hydrology, horizon, and
  spread.

Regular rectilinear lookups can provide a constant cell size. Irregular lookups
cannot generally be reduced to one global tuple. Algorithms that require a
constant finite-difference spacing should reject irregular sampling or provide
a different implementation. Neighbor-based algorithms can instead use
per-cell distances and bearings.

Geographic degree coordinates also require care. A degree step is not a
constant metric distance over latitude. The Raster adapter must retain the
required CRS or coordinate-system context inside the grid when metric results
are needed. Reducing `ras` to its raw parent without preserving that context is
insufficient.

The current geographic Raster cell-size implementation returns `(1, 1)` before
its intended conversion logic, making the remaining code unreachable. The
refactor should expose this as an explicit metric-policy decision rather than
carrying the hardcoded behavior forward accidentally.

## Algorithm capability groups

### Rectilinear metric or stencil algorithms

These algorithms require two ordered coordinate axes and, in many cases,
regular spacing:

- Horn and Zevenbergen–Thorne slope and aspect;
- Laplacian and curvature methods;
- hillshade and multihillshade;
- rugosity;
- progressive and simple morphological filters;
- matrix spread algorithms;
- horizon angle, sky-view factor, total viewshed, and viewshed;
- fixed-window entropy, percentile, and similar stencil operations.

They should dispatch on `RectilinearGrid` and produce a clear unsupported-grid
error for a Cells lookup unless a meaningful specialized implementation is
added.

### Topology-based neighborhood algorithms

These primarily need cells and neighbors:

- roughness;
- topographic position index;
- terrain ruggedness index;
- prominence;
- generic neighbor mapping;
- priority-flood depression filling.

They can share high-level logic across rectangular and DGG topologies, with
lookup-specific neighbor traversal underneath.

### Topology plus geometry algorithms

Hydrological algorithms combine traversal with metric and directional
semantics:

- cell area initialization;
- D8, DInf, and FD8 accumulation;
- HAND;
- depression volume;
- topographic wetness and stream power indices;
- drainage potential.

These should resolve topology, cell area, distance, bearing, and direction
encoding through grid-aware helpers. Composite indices should reuse already
resolved spatial parts when calling slope or accumulation.

### Lookup-transparent operations

Some operations use array values but not neighborhood geometry, for example
skewness balancing and depression-depth subtraction. They may carry lookup
context only so outputs are rebuilt consistently. They should not accept and
ignore spatial keywords unless the common boundary provides meaningful
behavior.

## Multiple aligned inputs

Operations taking more than one array need a single reference grid. Every
aligned input must be resolved and validated against that reference.

Affected APIs include:

- spread, with source points, initial values, and friction;
- flow accumulation, with accumulation and closed arrays;
- HAND, with a stream mask;
- in-place source/destination methods;
- morphological filters with array-valued parameters.

The implementation should validate more than size equality. Dimensions,
lookups, order, and coordinates must either match or be deliberately aligned.
Independently autodetecting and normalizing every operand without a common
reference can silently combine differently ordered grids.

## Output reconstruction

The observable API includes wrapper and metadata preservation. Decomposing an
input to its raw parent loses information required to reconstruct:

- Raster dimensions and their original order;
- name, metadata, CRS, mapped CRS, reference dimensions, and missing value;
- GeoArray affine transform and CRS;
- DGG lookup identity;
- tuple results such as accumulation plus flow directions;
- additional dimensions such as horizon direction or Band.

The public boundary must therefore retain the original input or another rebuild
context. A centralized rebuild hook should handle array results, tuple results,
mutation, and outputs with newly introduced dimensions.

Horizon output is an important test case: it adds a direction dimension while
preserving the source's spatial dimension order. The current Raster allocator
assumes exactly two original dimensions and appends Band after both; it does not
cover spatial axes embedded among Time or Band dimensions.

## Nonspatial dimensions

The initial grid types deliberately describe one surface at a time:

- `RectilinearGrid` requires a two-dimensional input containing X and Y;
- `CellGrid` requires a one-dimensional input whose only dimension is Cells.

Existing behavior such as PMF's singleton third dimension can remain as an
explicit public convenience that selects a two-dimensional surface before grid
construction. General Time/Band batching and Time × Cells layouts are deferred
rather than partially encoded in the first grid protocol.

## Extension implications

### Rasters and ArchGDAL

The current Rasters extension requires both Rasters and ArchGDAL, even though
dimension and lookup resolution does not require GDAL. The lookup adapter should
ideally load with Rasters alone, with ArchGDAL-dependent CRS behavior isolated
in a narrower extension.

### DiscreteGlobalGrids

The DGG extension should register `Cells` with the dimension oracle and define
capabilities for `(CellLookup,)`. Optimized adjacency and traversal methods can
remain specialized by lookup type. Normal package tests should load and exercise
the extension; DGG is currently a weak dependency but is absent from the test
target in [`Project.toml`](../../Project.toml).

### GeoArrays

GeoArray has no DimensionalData dimension tuple. Its geometry lives in an affine
transform and CRS. Its adapter should construct a `RectilinearGrid{(1, 2)}`
directly from that affine geometry.
Treating it as an ordinary matrix would preserve values but erase the geometry
needed for cell size and rebuilding.

### Eikonal

The Eikonal extension currently strips only the friction wrapper and assumes
that other inputs share its positional layout. It should receive already
aligned spatial parts, or explicitly validate all participating inputs before
handing parent storage to the solver.

## Dispatch and inference considerations

Rectilinear axis order is a tuple-valued type parameter. The invariant is
`(xaxis, yaxis)`, giving concrete types such as `RectilinearGrid{(1, 2)}` and
`RectilinearGrid{(2, 1)}`. `CellGrid` needs no axis parameter while Cells is
restricted to the only axis.

Avoid building runtime `Vector{Any}` collections of dimensions or selectors.
The adapter should construct a concrete grid at a function barrier. Runtime
dimension selection is acceptable at the boundary as long as the resulting
grid type is concrete for the numerical implementation.

Public methods should resolve keywords before internal positional dispatch.
Julia does not dispatch on keywords, so adding `spatialdims` alone cannot resolve
the existing overlap between `AbstractMatrix` and generic indexed-grid methods.
The method tables will need to be reorganized around the grid as an internal
positional argument.

Internal extensions should specialize consistently on the grid argument.
Mixing methods that are more specific in the Raster argument with methods that
are more specific in the grid argument can create crossed-specificity and
ambiguities.

## Migration plan

### Phase 1: establish the protocol

- define the spatial-dimension oracle and adapter layer;
- define spatial-dimension resolution and validation;
- define `RectilinearGrid{(xaxis, yaxis)}` and `CellGrid`;
- define the cell, neighborhood, geometry, storage, and allocation interfaces;
- add AbstractMatrix, Raster, and DGG adapters;
- add output-rebuild and input-alignment hooks;
- write parity and error-behavior tests without changing numerical kernels.

### Phase 2: move geometric primitives

- resolve cell size through the grid rather than directly from the original DEM;
- centralize axis orientation and lookup direction;
- move cells, neighbors, cell area, distance, and bearing behind grid-aware
  methods;
- preserve topology-specific optimized implementations.

### Phase 3: reorganize public APIs

- convert public methods into one-time normalization façades;
- move numerical work into private `(ras, grid, ...)` methods;
- preserve grid-specific default methods through grid dispatch;
- make composite functions call private resolved implementations.

### Phase 4: migrate algorithm families

- topology-based relative metrics and priority flood;
- terrain derivatives and plotting helpers;
- hydrology and DGG specializations;
- morphology and spread;
- horizon and extra-dimension allocation;
- lookup-transparent operations where common rebuilding adds value.

### Phase 5: remove obsolete coupling

- remove Raster-type aliases used only to expose lookup capabilities;
- remove DGG Base indexing overrides that are no longer required internally;
- separate Rasters-only lookup integration from ArchGDAL-dependent behavior;
- document the extension protocol for third-party dimension and lookup types.

## Test strategy

The most important invariant test is representation independence: the same
physical surface stored as X,Y and Y,X must produce equivalent coordinate-aware
results after aligning output dimensions.

Required test groups include:

- X,Y and Y,X rasters with unequal axis lengths and unequal cell steps;
- forward and reverse ordered X and Y lookups;
- explicit `spatialdims` and autodetection producing equivalent grids;
- informative errors for missing, duplicate, or unsupported spatial axes;
- regular versus irregular lookup capability checks;
- wrapper, metadata, CRS, missing-value, and dimension-order preservation;
- horizon output with an added direction dimension;
- observer and source-index behavior under alternate dimension order;
- multi-input alignment for spread, masks, and accumulation arrays;
- plain dense matrices retaining current behavior;
- views, custom axes, and offset-indexed `AbstractMatrix` inputs;
- DGG Cells autodetection and explicit selection;
- rejection of DGG Cells with any additional array axis;
- invalid `Cells` dimensions without a `CellLookup`;
- D8 on general DGG cells and IGeo7-specific DInf/FD8/HAND behavior;
- Raster lookup integration when ArchGDAL is not loaded;
- inference checks for default and explicitly selected dimensions;
- GPU/backend smoke tests for the selected order-handling strategy.

## Decisions to confirm before implementation

1. **Explicit cell-size semantics:** should a tuple always mean semantic X,Y,
   or should it follow storage/selection order?

   Recommendation: always semantic X,Y for rectilinear grids.

2. **Public oracle:** should the spatial-dimension oracle be exported and
   documented for downstream extensions?

   Recommendation: yes. It is the intended extension point for dimensions such
   as `DGG.Cells`.

3. **Meaning of `ras`:** should it remain wrapper-preserving, or be raw parent
   storage?

   Recommendation: keep enough wrapper/context at the boundary to support CRS,
   metadata, backend selection, and reliable rebuilding. Raw storage may still
   be extracted inside backend-specific kernels.

## Expected outcome

After this refactor, exported algorithms should no longer decide grid behavior
primarily from `AbstractMatrix` versus a concrete Raster alias. They receive
`(ras, grid)` and dispatch on Geomorphometry's own grid vocabulary.

This should make ordinary matrices remain simple, make X,Y and Y,X rasters
semantically equivalent, let DGG expose topology through `CellLookup`, and
provide a controlled path for future lookup types without replicating whole
algorithm entry points in every extension.
