# Changelog

## 0.4.0

A one-time breaking reorganization so that downstream codes migrate once.
The reorganization does not change numerical behavior. See [MIGRATION.md](MIGRATION.md).

### Changed (breaking)

- The geometry kinds `omega_h::OmegaHGeometry`, `uniform_grid::UniformGridGeometry<Dim>` and
  `point_cloud::PointCloudGeometry` are public types and their factories return them
  (`std::shared_ptr<const Kind>`); `omega_h::GetMesh(geometry)` is `OmegaHGeometry::GetMesh()`.
- Public-API trim: helpers that 0.3 declared in public headers are in
  `pcms::detail` (not API); the types a public signature or a protected hook
  names stay public (see MIGRATION.md).
- Public API: only interfaces and factory functions are public. Spaces are
  created on a `Geometry` (`CreateLagrangeSpace`, `CreateSplineSpace`,
  `CreatePolynomialReconstructionSpace`; `uniform_grid::CreateGeometry`,
  `point_cloud::CreateGeometry`, `omega_h::CreateGeometry`), plus
  `CreateEmptyLayout`, `CreateConservativeProjection`,
  `CreateControlVariateProjection` and the backend entry points in
  `pcms/backends/<name>/<name>.hpp`. There is one concrete `FunctionSpace`;
  the geometry kinds, layouts, discretizations and projections are
  internal.
- Headers: all C++ headers are `.hpp`; only public headers are installed, in
  `include/pcms/<library>/`; `pcms.h` is replaced by the umbrella
  `pcms/pcms.hpp`; core public headers no longer include Omega_h, MFEM,
  MeshFields, PETSc or ADIOS2.
- Namespaces: backend code is in `pcms::<name>`, internal code in
  `pcms::<library>`; the public API stays in `pcms::`.
- Omega_h types left the core public signatures (`MakeConstRank2View`,
  `CreateUniformGridFromMesh` and the Omega_h `OverlapMask` constructors moved
  to the Omega_h backend; `UniformGrid` point queries take any indexable
  point).
- `pcms::Mode` is a PCMS enum instead of an alias of `redev::Mode`.
- Point evaluation binds a source to its points once, on its geometry
  (`Geometry::Bind`, `CreatePointEvaluator(bound, op)`); `EvaluationRequest`,
  `TransformedSource` (the coordinate map is an argument of the bind),
  `InterpolationBase` and the `MLS*Interpolation` classes,
  `LocalizationFactory`, `FindOwningElements`, `GetOwningElementId(s)` and
  `MonteCarloSampling` are removed. The point search reports the entity a
  point lies on and its parametric coordinates on it.
- Python API mirrors the C++ factories; the class-based constructors are
  removed.
- Macros `MAX_DIM` and `ACTIVE_GPU_EXECUTION` are `PCMS_`-prefixed.
- Build: CMake 3.23 required; options renamed to `PCMS_ENABLE_<UPPERCASE>` /
  `PCMS_` form (old names still accepted in 0.4, deprecated); `pcms::core` is
  a deprecated alias of `pcms::pcms`; `pcms::transfer` no longer depends on
  the coupler, redev or ADIOS2; package files installed to
  `lib/cmake/pcms`; no baked dependency paths; `SameMinorVersion`
  compatibility; libraries have `VERSION`/`SOVERSION`; the SWIG interface
  files are no longer installed.
- Dependencies: Omega_h (and MPI, KokkosKernels, MeshFields, PETSc, MFEM) are
  private dependencies of the core libraries; a core-only consumer needs only
  Kokkos and PCMS. Users of the Omega_h backend header link
  `Omega_h::omega_h`. The package config finds private dependencies only for
  static builds. `pcms/utility/assert.hpp` no longer includes `<mpi.h>`.
- `SupportResults` holds Kokkos views instead of Omega_h arrays.
- `Discretization` is removed: layouts hold their `Geometry` (`FieldLayout::GetGeometry()`);
  `GetDiscretization()` and `EmptyDiscretization` are gone.

### Added

- `PCMS_ENABLE_COUPLER` to build without redev and ADIOS2.
- A build without Omega_h (`PCMS_ENABLE_OMEGA_H=OFF`, with the server, the
  coupler, XGC and MeshFields off): the uniform-grid Lagrange and spline
  spaces and point-cloud MLS (`point_cloud::CreateGeometry`,
  `CreatePolynomialReconstructionSpace`) need only Kokkos, KokkosKernels and
  MPI.
- `Geometry` (`Bind` with out-of-bounds policy, coordinate map and owner
  hints, `Rebind`, `Neighbourhood`, `Differential`), `BoundPoints`,
  `OwnerHints`, `QueryPoints` and the one-shot `Evaluate`
  (`pcms/field/query_points.hpp`).
- Differential operators: `CreatePointEvaluator(bound, op)` with
  `DifferentialOperator`, and `Gradient(source)`, `Divergence(source)`,
  `Curl(source)` (in `pcms::detail`).
- Factory functions and backend headers listed above; the Omega_h backend in
  Python (`pcms.omega_h`).
- Tests that every public header compiles on its own and that core public
  headers include no backend library; CI builds and runs
  `examples/external-usage-example` against the install tree.
