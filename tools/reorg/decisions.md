# PCMS breaking reorganization: decisions and amendments

This file records the decisions made with the maintainer before and during the
reorganization, including where they deliberately differ from the original task
description (`refactor.md`). Reviewers should read this first.

## Process

| Item | Decision |
|---|---|
| Base | `reorg/breaking-block` from `upstream/develop` (c33f816c) |
| PR #379, #401 | Still open; judged not to conflict materially; not blocking |
| PR #410 (CMake 3.23) | Will not be merged; the first commit of this branch bumps `cmake_minimum_required` to 3.23 (required by `FILE_SET`, already in use) |
| Version | Minor bump: 0.3.0 → 0.4.0 |
| Push / PR | Commits are local; the maintainer pushes and opens the draft PR |
| Verification per commit | Minimal-static and Full-shared build + `ctest` |
| Verification per phase boundary and at tip | Additionally Full-static, `cmake --install`, and `examples/external-usage-example` built against the install tree |
| Open branches (#400 dof-access, external-eval-infra, fusion-io, nimrod) | Reorg lands first; those branches migrate afterwards with `tools/migrate-pcms.sh` |
| Old CMake target names (`pcms::core`, ...) | Kept for one release as deprecated `INTERFACE` targets linking `pcms::pcms` |
| Old header paths | No forwarding stubs; `tools/reorg/header-map.txt` + `tools/migrate-pcms.sh` cover them |
| Fortran SWIG wrappers | Not regenerated: they include `"pcms/capi/*.h"`, which keeps its path, and the C API does not change |
| Commit granularity | Deviation from the plan's one-commit-per-sub-library rule, decided at review: the `.h`→`.hpp` renames, the `pcms/`-rooted includes, the moves to `include/` and the install rules are one commit each (each mechanical commit lists its commands in order; every one reproduces its diff). Fixes to the reorg scripts are folded into the commits that add the scripts, and the decision/inventory updates into the two docs commits. Backend moves and sub-library namespaces stay one commit per backend / library. |
| Verification | Per commit (Minimal-static for CMake changes, Full-shared for all) through Phase 5; after that, the end of each phase, the commits with temporary Python aliases and the tip, bisecting only on failure (decided at review) |

### Baseline limitation

On develop the coupler is built unconditionally (`find_package(redev REQUIRED)`,
`find_package(ADIOS2 REQUIRED)` at top level, `smoke_test` links `pcms::core`).
The "Minimal: no coupler" configuration therefore cannot be built until Phase 2.
Until then "Minimal" means: Omega_h only; meshfields, PETSc, MFEM, XGC, Python,
C, Fortran and transient off; static.

## Public API (amendment: the task originally excluded API redesign)

The change is limited to *how objects are constructed*; semantics of the
existing interfaces do not change.

### Principles

- Public headers contain **interfaces** and **free-function factories** only.
  All concrete implementation classes are private (stay in `src/`, not
  installed). In-tree targets and tests may still include private headers.
- Function spaces are created **by space type, from a `Discretization`**. Mesh
  sources produce discretizations.
- Codes with their own representation (XGC, MFEM in place, later NIMROD,
  Fusion-IO, BOUT++) get **per-code factories** returning the same interfaces.
- Only backend headers (`include/pcms/backends/<name>/`) may include or name
  backend-library types.
- Factories return `std::shared_ptr<Interface>`; transfer factories return
  `std::unique_ptr<TransferOperator<T>>`.

Rationale for interface-only (vs. keeping concrete classes public with factories
moved out): an audit of every non-override public method on concrete classes
found no external user that needs one. The only uses outside the owning
sub-library are the C API, the Python bindings and tests (XGC
`GetXGCLayout`/`SetLayoutName`, MFEM `OverlapMaskFromAttribute`,
`OmegaHDiscretization::GetMesh`); these become factory arguments or free
functions. Making a class public later is non-breaking; hiding it later is not.

### Public interfaces and types (`namespace pcms`)

Interfaces: `Discretization`, `FieldLayout`, `FunctionSpace`, `Field<T>` /
`Function<T>`, `FieldData` / `DOFAccess`, `PointEvaluator<T>` /
`PointEvaluatorFactory`, `TransferOperator<T>`, `CoordinateSystem` /
`CoordinateMap`, coupler types, `Timestepper`.

Backend-free types users construct directly: `Interpolator<T>`, `Copy<T>`,
`TransformedSource`, coordinate systems and maps, `FixedTimestepper`.

### Factories

```cpp
// Discretizations (core)
std::shared_ptr<const Discretization> CreateUniformGridDiscretization(const UniformGrid<Dim>&);
std::shared_ptr<const Discretization> CreatePointCloudDiscretization(...);
std::shared_ptr<const Discretization> CreateEmptyDiscretization(...);
std::shared_ptr<const FieldLayout>    CreateEmptyLayout(csys, ncomp);   // for out-of-tree evaluation-only spaces

// Discretizations (backend)
namespace pcms::backend::omega_h {
std::shared_ptr<const Discretization>
CreateDiscretization(::Omega_h::Mesh& mesh, std::string global_id_name = "global", /* owned mask */ ...);  // borrows mesh, as today
::Omega_h::Mesh& GetMesh(const Discretization&);
}

// Function spaces, by type (core)
std::shared_ptr<FunctionSpace> CreateLagrangeSpace(disc, order, ncomp, csys, LagrangeOptions = {});
        // LagrangeOptions: today's MeshFields/OmegaH evaluator choice, layout_name
std::shared_ptr<FunctionSpace> CreateSplineSpace(disc, ..., SplineOptions = {});
std::shared_ptr<FunctionSpace> CreatePolynomialReconstructionSpace(disc, ..., PolynomialReconstructionOptions = {});

// Adopted application spaces (backend headers)
namespace pcms::backend::xgc  { CreateFunctionSpace(..., layout_name); CreateField<T>(space, name, data, plane_comm); }
namespace pcms::backend::mfem { CreateFunctionSpace(::mfem::ParFiniteElementSpace&, ..., layout_name);
                                OverlapMaskFromAttribute(::mfem::ParMesh&, int); }

// Transfers (core header include/pcms/transfer/, always declared;
// throw "not available in this build" if the implementing backend is disabled)
std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(src_space, tgt_space, opts = {});
std::unique_ptr<TransferOperator<Real>> CreateControlVariateProjection(src_space, tgt_space, opts = {});
```

Space and transfer factories dispatch on the concrete discretization type and
call internal implementation functions directly (e.g. `detail::FromOmegaH`,
`detail::FromUniformGrid`). `Discretization` does not know about function
spaces (virtual hooks on `Discretization` were rejected as coupling the two).

Unsupported combinations (e.g. a Lagrange space on an XGC discretization) throw
with a message naming both types.

### Coupler

`pcms::Mode` (today `using Mode = redev::Mode`) becomes a real enum
`{Deferred, Synchronous}`, converted to `redev::Mode` inside the coupler. Public
signatures use `pcms::Mode`; redev no longer appears in public headers.

### Python

The Python API is changed to match the C++ shape in this release
(`pcms.backend.omega_h.create_discretization(mesh)`,
`pcms.create_lagrange_space(disc, ...)`); the old names are removed.

## Backends and build layout (amendment to Phase 3)

- Namespace `pcms::backend::<name>` (instead of `pcms::<name>_backend`); bridges
  `pcms::backend::omega_h_meshfields`. Inside `pcms::backend`, `mfem::` would
  find `pcms::backend::mfem`, so backend code writes `::mfem::`. Other library
  namespaces (Kokkos, Omega_h, ...) are not shadowed and are written without
  a leading `::` (amended at review: qualifying them too changed about 1000
  lines and caused most conflicts with open PRs). Nothing enforces this; an
  unqualified `mfem::` that resolves to `pcms::backend::mfem` fails to
  compile, so the mistake cannot go unnoticed.
- Backends are **directories, namespaces and headers**, not libraries:
  sources in `src/pcms/backends/<name>/`, public headers in
  `include/pcms/backends/<name>/`, compiled into the core target that uses them
  (`pcms::field`, `pcms::transfer`, `pcms::coupler`) behind
  `PCMS_ENABLE_<BE>` / `PCMS_HAS_<BE>`. A bridge is compiled when both backends
  are enabled.
- No `pcms::<be>` targets and no `pcms_backend_part()` helper. Reason: separate
  backend libraries only pay off for static-link dependency hygiene and
  `COMPONENTS` support, neither of which is in scope; they would force a
  dispatch layer above `pcms::field`.
- Downstream links **`pcms::pcms`**, the documented, stable entry point. Splitting
  backends into libraries later (e.g. when Omega_h becomes optional) is then an
  internal change.

## Deferred (follow-ups, non-breaking)

- `ReadDiscretization(comm, path, opts)` with file-type dispatch (needs a
  partitioning option; owns the mesh it reads). Further mesh sources (MFEM
  mesh → shared infrastructure, gmsh, BOUT++).
- **Multi-field sources** (EQDSK, NIMROD, Fusion-IO, VMEC): a `FieldSource`
  layer on top of spaces (`Quantities()`, `CreateFunction(quantity)`,
  `QuantityInfo`, `OpenEQDSK`, `OpenFieldSource(path)`), designed separately.
  In this reorg `MakeEQDSKField` only changes its struct member to
  `std::shared_ptr<FunctionSpace>`.
- Batch evaluation of several quantities in one call; EQDSK flux functions
  (needs a composition evaluator).
- Making Omega_h optional for core internals.

## Include guards (amendment to Phase 1.4)

`#pragma once` is not standard C++ and is **not used**. Every header keeps an
`#ifndef`/`#define`/`#endif` guard. Instead of switching to `#pragma once`,
guard names are normalized to a path-derived scheme
(`pcms/field/field.hpp` → `PCMS_FIELD_FIELD_HPP`) in one mechanical pass
**after Phase 5**, when header paths are final. Before that, only headers that
had no guard (or `#pragma once`) are fixed.

## Phase 3 notes (backends)

- Layout: `src/pcms/backends/<name>/` with a `CMakeLists.txt` that adds the
  backend's sources and headers to the core target that uses them (see
  "Backends and build layout"); namespaces `pcms::backend::<name>`. Done for
  `xgc`, `mfem`, `point_cloud`, `uniform_grid`, `omega_h` and the bridge
  `omega_h_meshfields`.
- **MeshFields is a bridge.** Every MeshFields file needs the Omega_h adapter
  (`MeshFieldsAdapterLayout` is built from an `OmegaHDiscretization`), so the
  whole MeshFields adapter is in `backends/omega_h_meshfields/`, built only
  when both Omega_h and MeshFields are enabled. There is no `backends/meshfields`.
- **`function_space/lagrange` stays in `field`.** The plan listed it under the
  Omega_h adapter, but it is the generic Lagrange space (it also serves uniform
  grids) behind `CreateLagrangeSpace`.
- **Omega_h-based transfer algorithms stay in `transfer`** (conservative and
  control-variate projection, the integrators, `mesh_intersection`); the plan
  only listed field/localization parts for the adapter. They are private.
- `point_cloud` and `uniform_grid` have no external dependency and are used by
  other backends (the MFEM layout builds a `PointCloudDiscretization`); this is
  not treated as a bridge.
- `SupportResults`, the result type of the `LocalizationFactory` interface,
  moved out of `adj_search.hpp` into `localization/support_results.hpp` before
  the adjacency search moved into the Omega_h adapter.
- `transfer/linear_interpolant.hpp` and `transfer/multidimarray.hpp` moved to
  `backends/uniform_grid/` (only the uniform-grid evaluator uses them; field
  must not include transfer headers).
- In a backend with its own nested `detail` namespace, references to
  `pcms::detail` names are written `pcms::detail::X`.

## Phase 0 checkpoint (approved)

The classification in `inventory.md` is approved with these decisions:

1. Coupler communicator templates (`field_communicator`, `field_layout_communicator`,
   `field_exchange_planner`, `global_communicator`, `field_serializer`) become
   installed-detail (`include/pcms/coupler/detail/`); redev remains a public
   dependency of `pcms::coupler` only.
2. XGC and MFEM factories are `pcms::backend::xgc::CreateFieldFactory` and
   `pcms::backend::mfem::CreateFieldFactory`, returning
   `std::shared_ptr<FieldFactory>` (they are not evaluable function spaces).
3. `SimpleFieldData` stays a public class.
4. `spline_interpolator.hpp` becomes private.
5. `CreateUniformGridBinaryField` is a public free function in the Omega_h
   backend header.
6. `utility/print.hpp` is public (`setStdout`/`setStderr` let applications
   redirect PCMS output); `bounding_box`, `array_mask`, `mpi_type` are private.
7. Localization stays private (its public API is deferred).

Further decisions: the external usage example is rewritten in Phase 9 to use
only the public API; Python tests/examples are updated to the new Python names;
`UniformGrid::ClosestCellID` takes any indexable point type (template) instead of
`Omega_h::Vector<dim>` (`std::array` is not usable in device code without
relaxed constexpr); `FindPETSc.cmake` stays installed because the package
config's `find_dependency(PETSc)` needs it.

## Phases 5-9 notes

- **Factories (Phase 5).** Space factories take a `Discretization` and call the
  existing builders with the same arguments (the space still builds its own
  layout and internal discretization from the mesh or grid). The Omega_h
  discretization carries the global-id tag name and the optional owned mask,
  which were `LagrangeFunctionSpace::FromMesh` arguments. Polynomial
  reconstruction over points takes the points directly; there is no public
  point-cloud discretization factory (nothing public would consume it).
  `CreateUniformGridBinaryField` returns only the `Field` (its layout type is
  private).
- **Reclassified while splitting:** `utility/bounding_box.hpp` is public
  (`UniformGrid::GetCellBBOX` returns `AABBox`); `transfer/monte_carlo_sampling.hpp`
  is public (enum used in the projection options).
- **Install rules:** every library's installed FILE_SET contains only its
  headers under `include/pcms`; private headers are private sources.
  `pcms::localization` switched from `PUBLIC_HEADER` (installed flat) to a
  FILE_SET. The C API install include directory is `include`.
- **PUBLIC dependencies that remain:** Kokkos; redev (coupler and its
  installed-detail templates); perfstubs (`utility/detail/profile.hpp`); MFEM
  (MFEM backend header); Omega_h (Omega_h and XGC backend headers, and core
  internals). MeshFields is PRIVATE.
- **Omega_h out of core public signatures (Phase 6):** no new vector type;
  `UniformGrid` point queries are templates on the point type.
  The Python bindings follow: `create_uniform_grid_from_mesh`,
  `compute_entity_centroids` and `map_entity_field_to_vertices_{average,mls}`
  take an Omega_h mesh and moved to `pcms.backend.omega_h`; the
  `compute_entity_centroids_few` alias is removed; `patch_size_for_degree`
  (no mesh) stays top-level (maintainer decision during review).
- **Sub-library namespaces (Phase 7):** private code of utility,
  discretization, localization, field, transfer and the coupler's internal
  classes (including its installed-detail templates) moved to `pcms::<lib>`.
  Transient has no internal code. `field/evaluator/interpolator_aliases.hpp`
  keeps its generic aliases (`member_type`, `ScratchMatView`, ...) in `pcms::`
  because other files define unrelated local aliases with those names.
- **Packaging (Phase 8):** package files in `lib/cmake/pcms`
  (`PCMS_INSTALL_CMAKEDIR`); `SameMinorVersion`; version 0.4.0, SOVERSION 0.4;
  SWIG `.i` files no longer installed; `FindPETSc.cmake` kept (needed by the
  config's `find_dependency(PETSc)`); baked `HINTS` removed; baked
  `PETSC_DIR`/`PETSC_ARCH` defaults kept (used by the installed FindPETSc).
- **Guards and migration (Phase 9):** see `tools/check-public-includes.py`,
  `test/public_headers/`, `tools/migrate-pcms.sh`, `MIGRATION.md`.
