# Migrating from PCMS 0.3 to 0.4

PCMS 0.4 reorganizes the library so that downstream codes migrate once: the
public API is now a small set of interfaces and factory functions, everything
else is internal. Behavior and numerical results do not change.

## Building and linking

- **Requirements:** CMake 3.23 or newer.
- **Link `pcms::pcms`.** It links every PCMS library that was built and is the
  supported way to use PCMS. `pcms::core` still works in 0.4 as a deprecated
  alias of `pcms::pcms` (CMake prints a deprecation warning) and will be
  removed. The other targets are `pcms::utility`, `pcms::discretization`,
  `pcms::localization`, `pcms::field`, `pcms::transfer`, `pcms::coupler`,
  `pcms::transient`, `pcms::capi`, `pcms::fortranapi`.
- **Optional capabilities are off by default.** A plain
  `cmake -S . -B build` needs only Kokkos and MPI and builds the core: the
  uniform-grid backend (Lagrange and spline spaces), the point-cloud geometry
  and its support search, coordinate systems and
  maps, QueryPoints, interpolation and copy, derivative sources, the
  timestepper, and the Python module with `-DPCMS_ENABLE_PYTHON=ON`. Everything
  else is a flag named after what it builds (the backend namespace, the Python
  module `pcms.<name>`, the key of `pcms.features`, and the macro in
  `pcms/configuration.hpp` all share the name):

  | Flag (`PCMS_ENABLE_...`) | Builds | Needs |
  |---|---|---|
  | `OMEGA_H` | the Omega_h backend | Omega_h |
  | `OMEGA_H_MESHFIELDS` | MeshFields spaces on Omega_h meshes | flag `OMEGA_H`, MeshFields |
  | `CONSERVATIVE_PROJECTION` | conservative and control-variate projections | flag `OMEGA_H_MESHFIELDS`, PETSc |
  | `POLYNOMIAL_RECONSTRUCTION` | polynomial-reconstruction (MLS) spaces on point clouds | KokkosKernels |
  | `XGC` | the XGC backend | flag `OMEGA_H` |
  | `MFEM` | the MFEM backend | MFEM |
  | `COUPLER` | the coupler | redev, ADIOS2 |
  | `C` | the C client API with `COUPLER` and `XGC`, the C transfer API with `OMEGA_H` and `POLYNOMIAL_RECONSTRUCTION` | one of the two pairs of flags |
  | `FORTRAN` | the Fortran API (the parts `C` builds) | flag `C` |
  | `PYTHON` | the Python module | pybind11 |
  | `EXAMPLES` | the in-tree examples | |

  A flag whose prerequisite flag is off stops the configuration with an error
  that names both; nothing is switched on or off silently. To get the previous
  full build, pass `-DPCMS_ENABLE_OMEGA_H=ON -DPCMS_ENABLE_OMEGA_H_MESHFIELDS=ON
  -DPCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION=ON -DPCMS_ENABLE_XGC=ON -DPCMS_ENABLE_COUPLER=ON
  -DPCMS_ENABLE_C=ON -DPCMS_ENABLE_FORTRAN=ON -DPCMS_ENABLE_EXAMPLES=ON` (and
  `-DPCMS_ENABLE_CONSERVATIVE_PROJECTION=ON` for what `PCMS_ENABLE_PETSC=ON`
  built). `PCMS_ENABLE_PRINT` stays `ON`.
- **Renamed options** (the old names still work in 0.4, with a deprecation
  warning): `PCMS_ENABLE_Python` → `PCMS_ENABLE_PYTHON`,
  `PCMS_ENABLE_Fortran` → `PCMS_ENABLE_FORTRAN`,
  `PCMS_ENABLE_MESHFIELDS` → `PCMS_ENABLE_OMEGA_H_MESHFIELDS`,
  `PCMS_ENABLE_PETSC` → `PCMS_ENABLE_CONSERVATIVE_PROJECTION`,
  `PETSC_LINK_STATIC` → `PCMS_PETSC_LINK_STATIC`. The macros
  `PCMS_ENABLE_MESHFIELDS` and `PCMS_ENABLE_PETSC` are still defined in
  `pcms/configuration.hpp` for 0.4.
- **Removed options:** `PCMS_ENABLE_SERVER` and `PCMS_ENABLE_CLIENT` (they
  gated nothing since the coupler communication merged) and
  `PCMS_ENABLE_TRANSIENT` (the timestepper is always built). Setting one
  prints a deprecation message and has no effect; remove it. The macros
  `PCMS_ENABLE_SERVER` and `PCMS_ENABLE_CLIENT` are no longer defined
  (`PCMS_ENABLE_TRANSIENT` still is, always).
- **Python:** the module no longer needs the coupler or Omega_h. `pcms.omega_h`
  (with `PCMS_ENABLE_OMEGA_H`) and the polynomial-reconstruction and
  projection factories exist only when built; `pcms.features` is a read-only dict of the
  capability flags (`pcms.features["omega_h"]`, ...) to test for them.
- **`CreatePolynomialReconstructionSpace`** (and `pcms.MLSOptions`,
  `pcms.RadialBasisFunction` and `pcms.create_polynomial_reconstruction_space`
  in Python) exist only with `PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION`; the
  point-cloud geometries (`point_cloud::CreateGeometry`,
  `omega_h::CreatePointCloudGeometry`) are always built.

- **Finding PCMS:** the package files are in `<prefix>/lib/cmake/pcms`. The
  package no longer bakes in the paths of its dependencies: make them findable
  as usual (`CMAKE_PREFIX_PATH`, `<Package>_ROOT` or `<Package>_DIR`) for
  Kokkos, KokkosKernels, Omega_h, and, if PCMS was built with them,
  redev, ADIOS2, MeshFields, MFEM and PETSc.
- **Dependencies of consumers:** Omega_h, KokkosKernels, MeshFields, PETSc,
  MFEM and MPI are private dependencies of the core libraries: a program that
  uses only the core API needs Kokkos and PCMS (and, for a static PCMS build,
  the private dependencies at link time; the package config finds them
  then). Code that includes `pcms/backends/omega_h/omega_h.hpp` or Omega_h
  itself links `Omega_h::omega_h` (`find_package(pcms)` finds Omega_h when
  the Omega_h backend is installed); code that calls MPI links `MPI::MPI_CXX`.
  `pcms::coupler` and the C/Fortran APIs still bring redev, ADIOS2 and MPI.
- **Versions:** `find_package(pcms 0.4)` accepts only 0.4.x
  (`SameMinorVersion`); shared libraries have SOVERSION `0.4`.

## Headers

- All C++ headers end in `.hpp` (the C API headers in `pcms/capi/` stay `.h`).
- Only the public headers are installed, under `include/pcms/<library>/`.
  `pcms/pcms.hpp` includes the whole public API; backend headers are included
  separately (`pcms/backends/omega_h/omega_h.hpp`,
  `pcms/backends/mfem/mfem.hpp`, `pcms/backends/xgc/xgc.hpp`).
- The top-level `pcms.h` is replaced by `pcms/pcms.hpp`.
- The geometry factories of the uniform-grid and point-cloud backends are in
  `pcms/backends/uniform_grid/geometry.hpp` and
  `pcms/backends/point_cloud/geometry.hpp`; `MLSOptions` is in
  `pcms/field/mls_options.hpp`.
- Installed headers below a `detail/` folder (`pcms/<library>/detail/`,
  `pcms/backends/<name>/detail/`) and the names in `detail` namespaces are not
  API: they may change in any release.
- Headers that are no longer installed are private; their classes must be
  replaced by the public API below.
- Core public headers no longer include Omega_h, MFEM, MeshFields, PETSc,
  ADIOS2 or MPI. Code that got Omega_h or MPI declarations through a PCMS
  header must include the headers it uses (`pcms/utility/assert.hpp` no
  longer includes `<mpi.h>`; `PCMS_ALWAYS_ASSERT(cond, comm[, msg])` needs it
  at the call site).
- `Discretization` is removed (it is not installed and no class refers to it).
  A layout holds the `Geometry` its DOF holders live on: `FieldLayout::GetGeometry()`
  (public, `std::shared_ptr<const Geometry>`), and `FunctionSpace::GetDiscretization()`,
  `FieldLayout::GetDiscretization()`, `Geometry::GetDiscretization()` and
  `EmptyDiscretization` are removed; use `GetGeometry()`. Custom layouts store the
  geometry they were built on and return it from `GetGeometry()`; custom `Geometry`
  kinds drop their `GetDiscretization` override.
- `SupportResults` (`pcms/geometry/supports.hpp`, returned by
  `Geometry::Neighbourhood`) holds Kokkos views in `DeviceMemorySpace`
  instead of Omega_h arrays: `supports_ptr` and `supports_idx` are
  `Kokkos::View<const LO*, DeviceMemorySpace>`, `radii2` is
  `Kokkos::View<Real*, DeviceMemorySpace>`. Indexing (`supports_ptr[i]`),
  `size()` and `data()` work as before; to build one from Omega_h arrays pass
  their `view()`.

## Creating spaces, fields and transfers

The concrete classes are internal. A space is created on a `Geometry`, which
owns the discretization (the entities the space lives on), the coordinate
system and the point search; use the factories, which return the interfaces
`Geometry`, `FunctionSpace`, `FieldFactory`, `FieldLayout` and
`TransferOperator`:

| 0.3 | 0.4 |
|---|---|
| `LagrangeFunctionSpace::FromMesh(mesh, order, n, csys, gid_name, backend, layout_name)` | `CreateLagrangeSpace(omega_h::CreateGeometry(mesh, csys, {.global_id_name = gid_name}), order, n, {.backend = backend, .layout_name = layout_name})` |
| `... FromMesh(mesh, order, n, csys, owned_mask, ...)` | `CreateLagrangeSpace(omega_h::CreateGeometry(mesh, csys, {.global_id_name = gid_name, .owned_mask = owned_mask}), ...)` |
| `LagrangeFunctionSpace::FromUniformGrid(grid, n, csys, order)` | `CreateLagrangeSpace(uniform_grid::CreateGeometry(grid, csys), order, n)` |
| `LagrangeFunctionSpace::Backend::{MeshFields,OmegaH}` | `LagrangeBackend::{MeshFields,OmegaH}` (in `LagrangeOptions`) |
| `SplineFunctionSpace::FromUniformGrid(grid, csys)` | `CreateSplineSpace(uniform_grid::CreateGeometry(grid, csys))` |
| `PolynomialReconstructionFunctionSpace::Create(points, opts)` | `CreatePolynomialReconstructionSpace(point_cloud::CreateGeometry(points, neighbourhood), opts)` |
| `PolynomialReconstructionFunctionSpace::FromMesh(mesh, dim, csys, opts)` | `CreatePolynomialReconstructionSpace(omega_h::CreatePointCloudGeometry(omega_h::CreateGeometry(mesh, csys), dim, neighbourhood), opts)` |
| `EmptyFieldLayout(csys, n)` | `CreateEmptyLayout(csys, n)` |
| `OmegaHConservativeProjection(src, tgt, mass)` | `CreateConservativeProjection(src, tgt, {mass})` |
| `OmegaHControlVariateProjection(src, tgt, samples, sampling, seed, mass)` | `CreateControlVariateProjection(src, tgt, {samples, seed, mass, to_source})` |
| `XGCFieldFactory(rc, n, layout_name)` | `xgc::CreateFieldFactory(rc, n, layout_name)` |
| `factory.CreateField<T>(name, make_unique<XGCFieldData<T>>(factory.GetXGCLayout(), basis, data, comm))` | `xgc::CreateField<T>(factory, name, data, comm, basis)` |
| `MFEMFieldFactory(pmesh, pfes, gf, csys, layout_name)` | `mfem::CreateFieldFactory(pmesh, pfes, gf, csys, layout_name)` |
| `MFEMLayout::OverlapMaskFromAttribute(pmesh, attr)` | `mfem::OverlapMaskFromAttribute(pmesh, attr)` |
| `OverlapMask(size, omega_h_mask)` | `OverlapMask(size, omega_h::OverlapMaskFromArray(omega_h_mask, size))` |
| `MakeConstRank2View(omega_h_array, dim)` | `omega_h::MakeConstRank2View(omega_h_array, dim)` |
| `CreateUniformGridFromMesh<dim>(mesh, ...)` | `omega_h::CreateUniformGridFromMesh<dim>(mesh, ...)` |
| `CreateUniformGridBinaryField<2>(mesh, ...).second` | `omega_h::CreateUniformGridBinaryField(mesh, ...)` |
| `Interpolator<T> interp(src, tgt[, policy \| bound])` | `auto interp = CreateInterpolator<T>(src, tgt[, policy \| bound])` |
| `Copy<T> copy(src, tgt)` | `auto copy = CreateCopy<T>(src, tgt)` |
| `csys::Cartesian::Create(dim)` | `csys::CreateCartesian(dim)` |
| `csys::CylindricalRThetaZ::Create()`, `csys::CylindricalRZ::Create()` | `csys::CreateCylindricalRThetaZ()`, `csys::CreateCylindricalRZ()` |
| `CartesianToCylindrical()`, `CylindricalToCartesian()` | `CreateCartesianToCylindrical()`, `CreateCylindricalToCartesian()` |
| `PoloidalPlaneEmbedding(theta0)` | `CreatePoloidalPlaneEmbedding(theta0)` |
| `FixedTimestepper(dt)` | `CreateFixedTimestepper(dt)` |

(`omega_h::`, `uniform_grid::`, `point_cloud::`, `xgc::` and `mfem::` are
`pcms::omega_h::` and so on; all other names are in `pcms::`.) Factories
for a combination PCMS does not support (e.g. a Lagrange space on a
point-cloud geometry) throw `pcms_error`. The projection factories throw in a
build without PETSc and MeshFields.

- Factories return the interface: `std::shared_ptr<const Geometry>`,
  `std::shared_ptr<const CoordinateSystem>`, `std::shared_ptr<const
  CoordinateMap>`, `std::unique_ptr<TransferOperator<T>>`,
  `std::unique_ptr<Timestepper>`. The concrete classes (the geometry kinds,
  `Interpolator`, `Copy`, the coordinate systems and maps,
  `FixedTimestepper`) are internal. A function that needs a particular kind
  checks it at runtime and throws `pcms_error` (e.g.
  `omega_h::CreatePointCloudGeometry` for a geometry that is not an Omega_h
  geometry).
- Objects hold what they need: the transfer factories take
  `std::shared_ptr<const PointEvaluatorFactory>` and
  `std::shared_ptr<const FunctionSpace>` (pass the pointer, not `*space`)
  and the operators keep them alive.
- The `FunctionSpace` constructor (for spaces implemented outside PCMS)
  takes the geometry as a required argument; it throws `pcms_error` if the
  geometry is null or is not the evaluator source's.
- A geometry must be owned by a `std::shared_ptr` (the factories return one);
  spaces and BoundPoints hold it, and what it refers to (the Omega_h mesh)
  must outlive it. Creating a geometry builds the mesh's containing-element
  search: create one per mesh and create every space on it, so that the
  spaces share the search and accept each other's BoundPoints.
- There are no discretization factories (`Discretization` is removed; layouts
  hold their geometry), and a geometry owns its
  global-id tag and ownership mask (`OmegaHGeometryOptions`).
  PCMS does not hand a mesh back: keep the Omega_h mesh you created the
  geometry from.
- The coordinate system belongs to the geometry; the space factories no
  longer take one.
- `MLSOptions` keeps `degree`, `lambda`, `tol`, `decay_factor` and `basis`;
  `radius`, `min_req_supports` and `adapt_radius` are `NeighbourhoodOptions`
  (`radius`, `min_supports`, `adapt_radius`) of the point-cloud geometry,
  whose support search is the geometry's.
- `MonteCarloSampling` is removed: the control-variate projection samples
  uniformly at random.

EQDSK moved to `pcms/eqdsk/eqdsk.hpp` (no longer included by
`pcms/pcms.hpp`; it is expected to move to FusionIO 2.0) together with
`EQDSKData` (from `pcms/utility/eqdsk.hpp`). `MakeEQDSKField(data)` is
`CreateEQDSKFunction(data)`, which returns the psi `Function<Real>` (it owns
its spline space: `GetSpace()`), and `MakeEQDSKField(filename)` is
`ReadEQDSK(filename)`, which returns `EQDSKFunction {data, psi}`.
`EQDSKField` and `EQDSKFieldWithData` are removed.

`UniformGrid::ClosestCellID` and `IsPointInBounds` accept any point type
indexable with `[i]` (`Kokkos::Array`, a subview, `Omega_h::Vector`).

## Evaluating sources at points

A source (`PointEvaluatorFactory`: a function space, or a gradient,
divergence or curl source) is bound to its query points once by its
geometry, and evaluators are created from the BoundPoints without searching
again. `EvaluationRequest` and its `From*` factories are removed, and so is
`TransformedSource`: a coordinate map is an argument of the bind.

| 0.3 | 0.4 |
|---|---|
| `source.CreatePointEvaluator<T>(EvaluationRequest::FromCoordinates(coords, policy))` | `source.CreatePointEvaluator<T>(source.GetGeometry()->Bind(coords, policy))` |
| `source.CreatePointEvaluator<T>(EvaluationRequest::FromFunctionSpace(space, policy))` | `source.CreatePointEvaluator<T>(source.GetGeometry()->Bind(space.GetLayout()->GetDOFHolderCoordinates(), policy))`, or `CreateInterpolator` |
| `TransformedSource(source, map)` | `geometry->Bind(coords, policy, map)`: the points are in the map's source coordinate system, are located at their images, and the evaluators re-express the values through the map |
| (none) | `source.CreatePointEvaluator<T>(bound, DifferentialOperator::Gradient)`; `Gradient(source)`, `Divergence(source)`, `Curl(source)` are sources whose evaluators write the operator applied to the field (`pcms/field/gradient_source.hpp`, public); `Supports(op)` tells which a source implements |
| (none) | `geometry->Bind(coords, OwnerHints{...}, policy, map)`: points whose owners are known (generated in the elements of the geometry) skip the search; `geometry->Rebind(previous, coords)` recycles the storage of a BoundPoints nothing else holds |
| (none) | `QueryPoints query(geometry, policy, map); query.Locate(coords); query.Evaluate(function, out, op)`: evaluation at points that move between calls; `Evaluate(function, points, out, op, policy, map)` (`pcms/field/query_points.hpp`) does one evaluation |
| (none) | `geometry->Neighbourhood(coords, NeighbourhoodOptions, owners)`: the supports of reconstruction spaces; `geometry->Differential(bound)`: the Jacobian at the bound points |

A `PointEvaluatorFactory` implemented outside PCMS implements
`GetGeometry()`, `Supports(op)` and `CreatePointEvaluatorImpl(type, bound,
op)` (`pcms/field/point_evaluator_factory.hpp`), and a `Geometry` implemented
outside PCMS implements `GetCoordinateSystem()` and
`BindImpl(request)` (`pcms/geometry/geometry.hpp` and `bound_points.hpp`:
`BoundPoints` is the per-point data of a bind). A `PointEvaluator`
implements `NumPoints()` and `EvaluateImpl` instead of overriding
`Evaluate`, which checks the output's shape and basis and fills the points
the BoundPoints did not find.

`TransformedSource`, `InterpolationBase`, `MLSMeshInterpolation`,
`MLSPointCloudInterpolation`, `LocalizationFactory` and its two
implementations, `FindOwningElements` (removed),
`GetOwningElementId(s)`, the COO assembly utilities and `MonteCarloSampling`
are removed; the replacements are the bind above, `CreateInterpolator` (also
for reconstruction spaces) and `Geometry::Neighbourhood`.

## Namespaces

- The public API stays in `pcms::` (coordinate systems in `pcms::csys::`).
- Each backend has one flat namespace, `pcms::<name>` (`omega_h`, `mfem`,
  `xgc`, `point_cloud`, `uniform_grid`, `omega_h_meshfields`); its public
  entry points are there (`omega_h::CreateGeometry`,
  `uniform_grid::CreateGeometry`, `point_cloud::CreateGeometry`, ...).
  There is no `pcms::backend::` and no Python `pcms.backend`.
- Internal code is in one namespace per library (`pcms::utility`,
  `pcms::localization`, `pcms::field`, `pcms::transfer`, `pcms::coupler`)
  or in its backend's namespace. Installed code that is not API is in
  `pcms::detail` (or `pcms::<backend>::detail`), in headers below `detail/`.
- MFEM's namespace `mfem` collides with `pcms::mfem`: inside `namespace pcms`,
  or after `using namespace pcms;`, write MFEM's as `::mfem::`.
- The MFEM and XGC backends support send and receive only so far (field
  factories for the coupler); they have no geometry, spaces or evaluators
  yet.

## Coupler

`pcms::Mode` is now a PCMS enum (`Deferred`, `Synchronous`) instead of an
alias of `redev::Mode`. `SendField`, `ReceiveField`, `SendData`,
`ReceiveData` and the handles' `Send`/`Receive` take `pcms::Mode`: replace
`redev::Mode::X` with `pcms::Mode::X` in those calls. The `OverlapMask`
constructors taking Omega_h arrays are gone (see the table).

## Macros

`MAX_DIM` → `PCMS_MAX_DIM`, `ACTIVE_GPU_EXECUTION` →
`PCMS_ACTIVE_GPU_EXECUTION`, and the configuration macro
`PCMS_ENABLE_Fortran` → `PCMS_ENABLE_FORTRAN`.

## Python

The Python module mirrors the C++ factories; the class-based constructors
are removed:

| 0.3 | 0.4 |
|---|---|
| `pcms.LagrangeFunctionSpace.from_mesh(mesh, order, n, csys, backend=...)` | `pcms.create_lagrange_space(pcms.omega_h.create_geometry(mesh, csys), order, n, backend=pcms.LagrangeBackend.X)` |
| `pcms.LagrangeFunctionSpace.from_uniform_grid(grid, n, csys, order)` | `pcms.create_lagrange_space(pcms.uniform_grid.create_geometry(grid, csys), order, n)` |
| `pcms.PolynomialReconstructionFunctionSpace.create(coords, opts)` | `pcms.create_polynomial_reconstruction_space(pcms.point_cloud.create_geometry(coords), opts)` |
| `pcms.PolynomialReconstructionFunctionSpace.from_mesh(mesh, dim, csys, opts)` | `pcms.create_polynomial_reconstruction_space(pcms.omega_h.create_point_cloud_geometry(pcms.omega_h.create_geometry(mesh, system=csys), dim), opts)` |
| `pcms.OmegaHConservativeProjection(src, tgt, ...)` | `pcms.create_conservative_projection(src, tgt, ...)` |
| `pcms.OmegaHControlVariateProjection(src, tgt, ...)` | `pcms.create_control_variate_projection(src, tgt, ...)` |
| `pcms.create_uniform_grid_binary_field(mesh, divs)` | `pcms.omega_h.create_uniform_grid_binary_field(mesh, divs)` |
| `pcms.create_uniform_grid_from_mesh(mesh, divs)` | `pcms.omega_h.create_uniform_grid_from_mesh(mesh, divs)` |
| `pcms.map_entity_field_to_vertices_average(mesh, values, dim)` | read the values as a P0 field (`pcms.create_lagrange_space(geometry, 0, n)`) and interpolate it with `pcms.create_interpolator` |
| `pcms.map_entity_field_to_vertices_mls(mesh, values, dim)` | `pcms.create_polynomial_reconstruction_space(pcms.omega_h.create_point_cloud_geometry(geometry, dim), opts)` and `pcms.create_interpolator` |
| `pcms.compute_entity_centroids(mesh, dim)` | the DOF-holder coordinates of a P0 space (`space.get_layout().get_dof_holder_coordinates()`) |
| `pcms.Interpolator(src, tgt[, policy \| bound])` | `pcms.create_interpolator(src, tgt[, policy \| bound])` |
| `pcms.Copy(src, tgt)` | `pcms.create_copy(src, tgt)` |
| `pcms.csys.Cartesian.create(dim)`, `pcms.csys.CylindricalRThetaZ.create()`, `pcms.csys.CylindricalRZ.create()` | `pcms.csys.create_cartesian(dim)`, `pcms.csys.create_cylindrical_r_theta_z()`, `pcms.csys.create_cylindrical_rz()` |
| `pcms.CartesianToCylindrical()`, `pcms.CylindricalToCartesian()` | `pcms.create_cartesian_to_cylindrical()`, `pcms.create_cylindrical_to_cartesian()` |
| `pcms.same_coordinate_system(a, b)`, `pcms.same_value_basis(a, b)` | `a == b` |
| `pcms.has_identity_metric(s)` | `s.has_orthogonal_basis() and s.has_unit_scale_factors()` |
| `pcms.patch_size_for_degree(d)` | removed (set `NeighbourhoodOptions.min_supports`) |

Geometries are opaque `pcms.Geometry` objects (the `pcms.OmegaHGeometry`
class is removed) and coordinate systems are `pcms.CoordinateSystem` objects
(compare them with `==`, or use `kind()`). An object keeps alive what it was
created from: a geometry its mesh, a space its geometry, an evaluator its
source and bound points, a transfer operator its source and target spaces.
PCMS does not hand a mesh back (`pcms.omega_h.get_mesh` is removed).
`pcms.Geometry.bind`, `pcms.QueryPoints`, `pcms.gradient` and
`pcms.evaluate` are the Python forms of the bind and the evaluation above.
The MLSOptions fields `radius`, `min_req_supports` and `adapt_radius` are
`pcms.NeighbourhoodOptions` fields (`radius`, `min_supports`,
`adapt_radius`) of the point-cloud geometry.

## Porting code that used internals

Notes for code written against 0.3 internals, such as open pull requests:

- CMake source lists: keep the reorganization's list and add the branch's
  new files.
- New headers: name them `.hpp`, decide public (`include/pcms/<lib>/`) or
  private (next to the sources), and list them accordingly.
- Files the reorganization rewrote beyond rename detection (git reports a
  modify/delete conflict): re-apply the branch's hunks to the new file, or
  delete it if the branch deleted it.
- Names that are no longer public: 0.3 declared a number of helpers in public
  headers; in 0.4 they are in `pcms::detail` (array copy and view helpers of
  `utility/arrays.hpp`, `TypeEnumFromType`, `ForEach`, the coordinate
  system and value-basis helpers, the backends' BoundPoints and evaluators,
  the coupler's communicators and planners). They are not API and may
  change; code that used one should use the public equivalent or copy the few
  lines.
- Extending PCMS: the types a public signature or a protected hook names stay
  public (`PointEvaluator`, `PointEvaluatorFactory`, `BindRequest`,
  `DifferentialView`, `SupportResults`, `OwnerHints`, `DOFHolderOwners`,
  `FieldData`, `DOFAccess`, `FieldFactory`, `CreateEmptyLayout`, the
  `FunctionSpace` constructor). The concrete helpers moved to `pcms::detail`
  (`SimpleFieldData`, `PointCloudLayout`, the backends' BoundPoints and
  evaluators, the coupler's communicators and planners, the MLS helpers).
- Code that calls the point search directly (pull requests #345, #416 and
  #417): the search `Result` now reports the entity the point lies on
  (`dimensionality`, `entity`; `OutsideId`/`IsOutside` for a point outside;
the search is now private, in `src/pcms/backends/omega_h/point_search.hpp`)
  and its `parametric` coordinates on that entity, not an element and its
  barycentric coordinates. The search does not choose an element: a consumer
  that needs one takes the lowest adjacent element for a point on a shared
  entity (the Omega_h geometry's binding does, or the element of a `OwnerHints`
  hint) from the binding's per-point `Results()`. `FindOwningElements` and
  `GetOwningElementId(s)` are removed. `UniformGridPointSearch` is the search
  of the uniform-grid geometry (the cell and the position in it).
- Code that uses the (internal) integrators:
  `LinearFormIntegrator::GetIntegrationPoints()` returns global coordinates
  in the mesh's coordinate system, and the Omega_h integrators'
  `IntegrationPointOwners()` returns the element of each point as
  `OwnerHints`, to bind them without a search. `intersectTargets` takes the
  source `OmegaHGeometry`.


## C and Fortran APIs

Unchanged, except that the C headers are installed in `include/pcms/capi/`
and the installed include directory is `include` (include them as
`pcms/capi/<name>.h`, as the build tree already required).
