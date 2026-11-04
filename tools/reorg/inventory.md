# PCMS reorganization inventory (Phase 0)

Baseline: `upstream/develop` c33f816c. Decisions referenced below are in
[`decisions.md`](decisions.md).

## 1. Baseline build and test results

Environment: container image derived from the PCMS dev image (Spack: GCC 13.3.1,
MPICH 5.0.2, Kokkos 4.7.04 (Serial+OpenMP), Omega_h 11.4.0-scorec, meshfields
main, PETSc 3.26.0, redev main, ADIOS2 2.12.1, MFEM 4.10.0, Python 3.14,
pybind11 3.0.4). `RelWithDebInfo`, `ctest -j1`.

| Configuration | Options | Result |
|---|---|---|
| Minimal-static | `BUILD_SHARED_LIBS=OFF`, Omega_h on; meshfields, PETSc, MFEM, XGC, Python, C, Fortran, transient off | build OK, **187/187** tests pass |
| Full-shared | `BUILD_SHARED_LIBS=ON`, everything on (Omega_h, meshfields, PETSc, XGC, MFEM, C, Fortran, Python, transient, server, client) | build OK, **257/257** tests pass |
| Install | `cmake --install` of Full-shared; `examples/external-usage-example` built against the install tree | OK |
| Full-static | as Full-shared, `BUILD_SHARED_LIBS=OFF`, Python off | build OK, **257/257** tests pass |

Limitation: the coupler cannot be disabled on develop (`redev` and `ADIOS2` are
`REQUIRED` at top level, `smoke_test` links `pcms::core`), so the "Minimal: no
coupler" configuration does not exist yet. It becomes buildable in Phase 2.

## 2. Target dependency graph (before)

Generated with `cmake --graphviz` on the Full-shared build (executables and
tests excluded); the raw file is [`deps-before.dot`](deps-before.dot).
Solid = PUBLIC, dotted = PRIVATE, dashed = INTERFACE.

```
pcms                         -> Kokkos::kokkos (private), MPI::MPI_C (private), PETSc::PETSc (private), Python::Module (private), meshfields::meshfields (private), pcms_core (private), pcms_transfer (private), pybind11::headers (private), pybind11::module (private)
pcms_capi                    -> pcms_capi_core (interface), pcms_capi_transfer (interface)
pcms_capi_core               -> MPI::MPI_C, pcms_core (private)
pcms_capi_transfer           -> MPI::MPI_C, pcms_transfer (private)
pcms_core                    -> Kokkos::kokkos, MPI::MPI_CXX, pcms_coupler, pcms_field, pcms_localization, pcms_utility, perfstubs, redev::redev
pcms_coupler                 -> Kokkos::kokkos, MPI::MPI_CXX, Omega_h::omega_h, meshfields::meshfields, pcms_field, pcms_utility, perfstubs, redev::redev
pcms_discretization          -> Kokkos::kokkos, Omega_h::omega_h (private), pcms_utility
pcms_field                   -> Kokkos::kokkoskernels (private), Omega_h::omega_h, meshfields::meshfields, mfem, pcms_discretization, pcms_localization, pcms_utility
pcms_fortranapi              -> pcms_fortranapi_core (interface), pcms_fortranapi_transfer (interface)
pcms_fortranapi_core         -> Kokkos::kokkos, pcms_capi_core (private)
pcms_fortranapi_transfer     -> Kokkos::kokkos, pcms_capi_transfer (private)
pcms_localization            -> Kokkos::kokkos, Omega_h::omega_h, pcms_discretization, pcms_utility
pcms_pcms                    -> pcms_capi (interface), pcms_core (interface), pcms_fortranapi (interface), pcms_transfer (interface), pcms_transient (interface)
pcms_transfer                -> Omega_h::omega_h, PETSc::PETSc (private), meshfields::meshfields (private), pcms_core, pcms_field, pcms_localization
pcms_transient               -> pcms_utility
pcms_utility                 -> Kokkos::kokkos, Omega_h::omega_h, perfstubs
```

Layering problems visible in the graph:

- `pcms_transfer` links `pcms_core` (= coupler + redev + ADIOS2).
- `pcms_coupler` links Omega_h and meshfields; `pcms_utility` links Omega_h.
- `pcms::core` is "coupler plus everything"; the C API's `capi_core` links it.

## 2b. Target dependency graph after Phase 2

Raw file: [`deps-after-phase2.dot`](deps-after-phase2.dot) (Full configuration).

```
pcms                         -> Kokkos::kokkos (private), MPI::MPI_C (private), PETSc::PETSc (private), Python::Module (private), meshfields::meshfields (private), pcms_coupler (private), pcms_transfer (private), pybind11::headers (private), pybind11::module (private)
pcms_capi                    -> pcms_capi_core (interface), pcms_capi_transfer (interface)
pcms_capi_core               -> MPI::MPI_C, pcms_coupler (private)
pcms_capi_transfer           -> MPI::MPI_C, pcms_transfer (private)
pcms_core                    -> pcms_pcms (interface)
pcms_coupler                 -> Kokkos::kokkos, MPI::MPI_CXX, Omega_h::omega_h, meshfields::meshfields, pcms_field, pcms_utility, perfstubs, redev::redev
pcms_discretization          -> Kokkos::kokkos, Omega_h::omega_h (private), pcms_utility
pcms_field                   -> Kokkos::kokkoskernels (private), Omega_h::omega_h, meshfields::meshfields, mfem, pcms_discretization, pcms_localization, pcms_utility
pcms_fortranapi              -> pcms_fortranapi_core (interface), pcms_fortranapi_transfer (interface)
pcms_fortranapi_core         -> Kokkos::kokkos, pcms_capi_core (private)
pcms_fortranapi_transfer     -> Kokkos::kokkos, pcms_capi_transfer (private)
pcms_localization            -> Kokkos::kokkos, Omega_h::omega_h, pcms_discretization, pcms_utility
pcms_pcms                    -> pcms_capi (interface), pcms_coupler (interface), pcms_discretization (interface), pcms_field (interface), pcms_fortranapi (interface), pcms_localization (interface), pcms_transfer (interface), pcms_transient (interface), pcms_utility (interface)
pcms_transfer                -> Omega_h::omega_h, PETSc::PETSc (private), meshfields::meshfields (private), pcms_field, pcms_localization
pcms_transient               -> pcms_utility
pcms_utility                 -> Kokkos::kokkos, Omega_h::omega_h, perfstubs
```

Changes:

- `pcms_transfer` no longer links `pcms_core`; transfer and everything below
  it is free of redev and ADIOS2.
- `pcms::core` is a deprecated interface alias of `pcms::pcms`; `pcms::pcms`
  links every library that was built.
- `pcms::localization` does not use any field header: the coordinate types moved
  to utility, the field-layout adapters moved to field.
- `PCMS_ENABLE_COUPLER=OFF` builds without redev and ADIOS2 (no
  `pcms_coupler`, `pcms_capi_core`, `pcms_fortranapi_core`). Minimal
  configuration (static, coupler off): 172/172 tests pass.

Remaining in the graph, handled in later phases: `pcms_utility` and
`pcms_localization` link Omega_h publicly (Phases 3, 5, 6);
`pcms_coupler` links Omega_h and meshfields publicly (Phase 5).

## 2c. Target dependency graph at the end of the reorganization

Raw file: [`deps-after.dot`](deps-after.dot) (Full configuration, shared).

```
pcms                     -> Kokkos::kokkos (private), MPI::MPI_C (private), PETSc::PETSc (private), Python::Module (private), meshfields::meshfields (private), pcms_coupler (private), pcms_transfer (private), pybind11::headers (private), pybind11::module (private)
pcms_capi                -> pcms_capi_core (interface), pcms_capi_transfer (interface)
pcms_capi_core           -> MPI::MPI_C, pcms_coupler (private)
pcms_capi_transfer       -> MPI::MPI_C, pcms_transfer (private)
pcms_core                -> pcms_pcms (interface)
pcms_coupler             -> Kokkos::kokkos, MPI::MPI_CXX, Omega_h::omega_h, meshfields::meshfields (private), pcms_field, pcms_utility, perfstubs, redev::redev
pcms_discretization      -> Kokkos::kokkos, Omega_h::omega_h (private), pcms_utility
pcms_field               -> Kokkos::kokkoskernels (private), Omega_h::omega_h, meshfields::meshfields (private), mfem, pcms_discretization, pcms_localization, pcms_utility
pcms_fortranapi          -> pcms_fortranapi_core (interface), pcms_fortranapi_transfer (interface)
pcms_fortranapi_core     -> Kokkos::kokkos, pcms_capi_core (private)
pcms_fortranapi_transfer -> Kokkos::kokkos, pcms_capi_transfer (private)
pcms_localization        -> Kokkos::kokkos, Omega_h::omega_h, pcms_discretization, pcms_utility
pcms_pcms                -> pcms_capi (interface), pcms_coupler (interface), pcms_discretization (interface), pcms_field (interface), pcms_fortranapi (interface), pcms_localization (interface), pcms_transfer (interface), pcms_transient (interface), pcms_utility (interface)
pcms_public_header_check -> mfem (private), pcms_pcms (private)
pcms_transfer            -> Omega_h::omega_h, PETSc::PETSc (private), meshfields::meshfields (private), pcms_field, pcms_localization
pcms_transient           -> pcms_utility
pcms_utility             -> Kokkos::kokkos, Omega_h::omega_h, perfstubs
```

Changes since 2b: `pcms_coupler` and `pcms_field` link MeshFields
privately (Phase 5). Omega_h stays a public dependency of `pcms_utility`,
`pcms_localization`, `pcms_field`, `pcms_transfer` and `pcms_coupler`:
the Omega_h and XGC backend headers are public and compiled into those
targets, and Omega_h is still a required dependency. No core public header
includes Omega_h (checked by `tools/check-public-includes.py`).
`pcms_public_header_check` is the test target of
`test/public_headers/`.

### Final public headers (57 files under `include/pcms/`, plus the generated `pcms/version.hpp` and `pcms/configuration.hpp`)

- `(top level)`: `pcms.hpp`
- `backends/mfem`: `mfem.hpp`
- `backends/omega_h`: `omega_h.hpp`
- `backends/xgc`: `reverse_classification.hpp`, `xgc.hpp`
- `capi`: `client.h`, `interpolator.h`, `kokkos.h`, `mesh.h`
- `coupler`: `coupler.hpp`, `coupler_types.hpp`, `overlap_mask.hpp`
- `coupler/detail`: `field_communicator.hpp`, `field_exchange_planner.hpp`, `field_layout_communicator.hpp`, `field_serializer.hpp`, `global_communicator.hpp`
- `discretization`: `discretization.hpp`
- `field`: `coordinate_map.hpp`, `eqdsk_field.hpp`, `evaluation_request.hpp`, `field.hpp`, `field_data.hpp`, `field_factory.hpp`, `field_layout.hpp`, `function_space.hpp`, `function_spaces.hpp`, `out_of_bounds_policy.hpp`, `point_evaluator.hpp`, `point_evaluator_factory.hpp`, `simple_field_data.hpp`, `transformed_source.hpp`, `value_view.hpp`
- `localization`: `mls_options.hpp`
- `transfer`: `copy.hpp`, `interpolator.hpp`, `mass_matrix_type.hpp`, `monte_carlo_sampling.hpp`, `projection.hpp`, `transfer_operator.hpp`
- `transient`: `timestepper.hpp`
- `utility`: `arrays.hpp`, `assert.hpp`, `bounding_box.hpp`, `coordinate_system.hpp`, `coordinate_view.hpp`, `entity_types.hpp`, `eqdsk.hpp`, `memory_spaces.hpp`, `print.hpp`, `types.hpp`, `uniform_grid.hpp`
- `utility/coordinate_systems`: `cartesian.hpp`, `cylindrical.hpp`
- `utility/detail`: `common.hpp`, `inclusive_scan.hpp`, `profile.hpp`

## 3. Installed files

127 files under `include/pcms` (122 headers, 5 files of
`external/flibhpc`: `LICENSE`, `README.md`, `mpi.i`, `openacc.i`, `thrust.i`),
plus `lib64/cmake/pcms/FindPETSc.cmake`. Every C++ header in `src/pcms`
(outside the Python and Fortran bindings) is installed except
`localization/localization_path_selection.h`. The per-header list with
backend includes and the proposed classification is in section 5.

Headers that include a backend library directly: **44**. Headers that reach
one only through other PCMS headers: **50**, almost all through
`utility/arrays.h` → `<Omega_h_array.hpp>`. Removing that one include
(`MakeConstRank2View(Omega_h::Read<T>, int)` moves to the Omega_h backend)
removes Omega_h from the transitive closure of most headers.

## 4. Omega_h types in public signatures

All installed headers using `Omega_h::` are listed in section 5 (column
"Backend includes"). Under the proposed classification, only these
declarations remain in headers that stay installed:

| Header | Line | Declaration |
|---|---|---|
| `pcms/coupler/overlap_mask.h` | 42 | `OverlapMask(size_t size, Omega_h::HostRead<Omega_h::I8> is_overlap_host)` |
| `pcms/coupler/overlap_mask.h` | 51 | `OverlapMask(size_t size, Omega_h::Read<Omega_h::I8> is_overlap_device)` |
| `pcms/discretization/discretization/xgc_reverse_classification.h` | 111 | `template <typename T = Omega_h::LO>` |
| `pcms/discretization/discretization/xgc_reverse_classification.h` | 113 | `Omega_h::Mesh& mesh, std::string numbering = "simNumbering",` |
| `pcms/discretization/discretization/xgc_reverse_classification.h` | 116 | `auto class_ids_h = Omega_h::HostRead<Omega_h::ClassId>(` |
| `pcms/discretization/discretization/xgc_reverse_classification.h` | 120 | `auto vert_id = Omega_h::HostRead<T>(mesh.get_array<T>(0, numbering));` |
| `pcms/utility/arrays.h` | 248 | `auto MakeConstRank2View(Omega_h::Read<T> array, int dim)` |
| `pcms/utility/uniform_grid.h` | 36 | `ClosestCellID(const Omega_h::Vector<dim>& point) const` |
| `pcms/utility/uniform_grid.h` | 178 | `UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,` |
| `pcms/utility/uniform_grid.h` | 182 | `auto bbox = Omega_h::get_bounding_box<dim>(&mesh);` |
| `pcms/utility/uniform_grid.h` | 212 | `UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,` |

Handling:

- `OverlapMask`: add a constructor from a pcms `Rank1View`; the Omega_h
  constructors become free functions in the Omega_h backend header.
- `xgc_reverse_classification.h`: `ConstructRCFromOmegaHMesh` moves to an
  Omega_h+XGC bridge header; the rest has no Omega_h.
- `MakeConstRank2View(Omega_h::Read<T>, int)`: moves to the Omega_h backend.
- `UniformGrid::ClosestCellID(Omega_h::Vector<dim>)`: takes the new pcms
  fixed-size vector type (Phase 6). `CreateUniformGridFromMesh(Omega_h::Mesh&)`
  moves to the Omega_h backend.

The other uses of `Omega_h::` (188 declaration lines in 30 headers) are in
headers proposed as private.

## 5. Proposed header classification

- **public**: `include/pcms/<lib>/`, no backend library includes.
- **backend public**: `include/pcms/backends/<name>/`, may include its library.
- **installed-detail**: `include/pcms/<lib>/detail/`, required by public
  templates.
- **private**: stays under `src/` (or moves to `src/pcms/backends/<name>/`),
  not installed. In-tree targets and tests may still include it.

New public headers created for the factories (see `decisions.md`):
`include/pcms/backends/omega_h/omega_h.hpp`, `include/pcms/backends/mfem/mfem.hpp`,
`include/pcms/backends/xgc/xgc.hpp`, plus the space and projection factory
declarations in `include/pcms/field/` and `include/pcms/transfer/`.

| Class | Count |
|---|---|
| public | 43 |
| backend public | 1 |
| installed-detail | 8 |
| private | 70 |

| Current installed header | Class | Proposed location | Backend includes (direct / via) | Notes |
|---|---|---|---|---|
| `pcms/capi/client.h` | public | `include/pcms/capi/client.h` | – / – | C API, stays .h |
| `pcms/capi/interpolator.h` | public | `include/pcms/capi/interpolator.h` | – / – | C API, stays .h |
| `pcms/capi/kokkos.h` | public | `include/pcms/capi/kokkos.h` | – / – | C API, stays .h |
| `pcms/capi/mesh.h` | public | `include/pcms/capi/mesh.h` | – / – | C API, stays .h |
| `pcms/configuration.h` | public | `include/pcms/configuration.hpp` | – / – | generated |
| `pcms/coupler/coupler.hpp` | public | `include/pcms/coupler/coupler.hpp` | – / Omega_h/redev | ApplicationComm, handles; includes the communicator templates below (redev) |
| `pcms/coupler/coupler_types.h` | public | `include/pcms/coupler/coupler_types.hpp` | redev / – | Mode becomes a pcms enum |
| `pcms/coupler/field_communicator.hpp` | installed-detail | `include/pcms/coupler/detail/field_communicator.hpp` | redev / Omega_h | template used by coupler.hpp; uses redev (see open question) |
| `pcms/coupler/field_exchange_planner.h` | installed-detail | `include/pcms/coupler/detail/field_exchange_planner.hpp` | redev / Omega_h | used by coupler.hpp; redev |
| `pcms/coupler/field_layout_communicator.h` | installed-detail | `include/pcms/coupler/detail/field_layout_communicator.hpp` | redev / Omega_h | used by coupler.hpp; redev |
| `pcms/coupler/field_serializer.h` | installed-detail | `include/pcms/coupler/detail/field_serializer.hpp` | – / Omega_h | used by coupler.hpp |
| `pcms/coupler/global_communicator.h` | installed-detail | `include/pcms/coupler/detail/global_communicator.hpp` | redev / – | used by coupler.hpp; redev |
| `pcms/coupler/overlap_mask.h` | public | `include/pcms/coupler/overlap_mask.hpp` | Omega_h / – | argument of ApplicationComm::SetLayoutOverlapMask; Omega_h ctors move to backend |
| `pcms/coupler/partition.h` | private | `src/pcms/coupler/partition.hpp` | redev / – | internal |
| `pcms/discretization/discretization.h` | public | `include/pcms/discretization/discretization.hpp` | – / Omega_h | interface + new Create*Discretization factories |
| `pcms/discretization/discretization/empty.hpp` | private | `src/pcms/discretization/empty.hpp` | – / Omega_h | impl; factory CreateEmptyDiscretization |
| `pcms/discretization/discretization/omega_h.hpp` | private | `src/pcms/backends/omega_h/discretization.hpp` | Omega_h / – | impl; factory backend::omega_h::CreateDiscretization |
| `pcms/discretization/discretization/point_cloud.hpp` | private | `src/pcms/backends/point_cloud/discretization.hpp` | – / Omega_h | impl; factory CreatePointCloudDiscretization |
| `pcms/discretization/discretization/uniform_grid.hpp` | private | `src/pcms/backends/uniform_grid/discretization.hpp` | – / Omega_h | impl; factory CreateUniformGridDiscretization |
| `pcms/discretization/discretization/xgc.hpp` | private | `src/pcms/backends/xgc/discretization.hpp` | – / Omega_h | impl |
| `pcms/discretization/discretization/xgc_reverse_classification.h` | backend public | `include/pcms/backends/xgc/reverse_classification.hpp` | Omega_h / – | ReverseClassificationVertex is an XGC input (C API); ConstructRCFromOmegaHMesh moves to an omega_h+xgc bridge |
| `pcms/field/coordinate_map.hpp` | public | `include/pcms/field/coordinate_map.hpp` | – / Omega_h |  |
| `pcms/field/coordinate_system.hpp` | public | `include/pcms/field/coordinate_system.hpp` | – / – |  |
| `pcms/field/coordinate_systems/cartesian.hpp` | public | `include/pcms/field/coordinate_systems/cartesian.hpp` | – / – |  |
| `pcms/field/coordinate_systems/cylindrical.hpp` | public | `include/pcms/field/coordinate_systems/cylindrical.hpp` | – / – |  |
| `pcms/field/coordinate_view.hpp` | public | `include/pcms/field/coordinate_view.hpp` | – / Omega_h |  |
| `pcms/field/data/mesh_fields.h` | private | `src/pcms/backends/meshfields/field_data.hpp` | – / Omega_h/meshfields |  |
| `pcms/field/data/mfem.hpp` | private | `src/pcms/backends/mfem/field_data.hpp` | MFEM / Omega_h |  |
| `pcms/field/data/point_cloud.h` | private | `src/pcms/backends/point_cloud/field_data.hpp` | – / Omega_h |  |
| `pcms/field/data/simple.h` | public | `include/pcms/field/simple_field_data.hpp` | – / Omega_h | generic FieldData over user views; used by C API and tests to wrap user arrays (question: public or replace with a factory) |
| `pcms/field/data/xgc.h` | private | `src/pcms/backends/xgc/field_data.hpp` | – / Omega_h | replaced by backend::xgc::CreateField |
| `pcms/field/element_dispatch.h` | private | `src/pcms/field/element_dispatch.hpp` | – / – |  |
| `pcms/field/eqdsk_field.h` | public | `include/pcms/field/eqdsk_field.hpp` | – / Omega_h | MakeEQDSKField; struct member becomes shared_ptr<FunctionSpace> |
| `pcms/field/evaluation_request.h` | public | `include/pcms/field/evaluation_request.hpp` | – / Omega_h |  |
| `pcms/field/evaluator/mesh_fields.h` | private | `src/pcms/backends/meshfields/evaluator.hpp` | – / Omega_h/meshfields |  |
| `pcms/field/evaluator/mesh_fields_backend.h` | private | `src/pcms/backends/omega_h_meshfields/evaluator_backend.hpp` | meshfields / Omega_h | needs Omega_h + meshfields: bridge |
| `pcms/field/evaluator/mls_interpolation.hpp` | private | `src/pcms/field/evaluator/mls_interpolation.hpp` | Omega_h / – | RadialBasisFunction enum moves to mls_options |
| `pcms/field/evaluator/mls_interpolation_impl.hpp` | private | `src/pcms/field/evaluator/mls_interpolation_impl.hpp` | Omega_h / – |  |
| `pcms/field/evaluator/mls_options.h` | public | `include/pcms/field/mls_options.hpp` | – / Omega_h | options for polynomial reconstruction (Python binds it); stop including mls_interpolation.hpp |
| `pcms/field/evaluator/mls_point_cloud.h` | private | `src/pcms/backends/point_cloud/mls_evaluator.hpp` | Omega_h / – |  |
| `pcms/field/evaluator/omega_h_lagrange.h` | private | `src/pcms/backends/omega_h/lagrange_evaluator.hpp` | Omega_h / – |  |
| `pcms/field/evaluator/pcms_interpolator_aliases.hpp` | private | `src/pcms/field/evaluator/interpolator_aliases.hpp` | – / Omega_h | stutter dropped (Phase 4) |
| `pcms/field/evaluator/pcms_interpolator_logger.hpp` | private | `src/pcms/field/evaluator/interpolator_logger.hpp` | – / Omega_h | stutter dropped (Phase 4) |
| `pcms/field/evaluator/pcms_interpolator_view_utils.hpp` | private | `src/pcms/field/evaluator/interpolator_view_utils.hpp` | Omega_h / – | stutter dropped (Phase 4) |
| `pcms/field/evaluator/point_cloud.h` | private | `src/pcms/backends/point_cloud/evaluator.hpp` | Omega_h / – |  |
| `pcms/field/evaluator/spline_interpolator.hpp` | private | `src/pcms/backends/uniform_grid/spline_interpolator.hpp` | – / Omega_h | installed on purpose in an earlier PR ("install spline_interpolator.hpp"); question: any external user? |
| `pcms/field/evaluator/uniform_grid.h` | private | `src/pcms/backends/uniform_grid/evaluator.hpp` | – / Omega_h |  |
| `pcms/field/evaluator/uniform_grid_spline.h` | private | `src/pcms/backends/uniform_grid/spline_evaluator.hpp` | – / Omega_h |  |
| `pcms/field/field.h` | public | `include/pcms/field/field.hpp` | – / Omega_h |  |
| `pcms/field/field_data.h` | public | `include/pcms/field/field_data.hpp` | – / Omega_h | FieldData / DOFAccess interfaces |
| `pcms/field/field_evaluator_factory.h` | private | `src/pcms/field/field_evaluator_factory.hpp` | – / Omega_h | internal evaluator-factory interface (expert access kept internal) |
| `pcms/field/field_factory.h` | public | `include/pcms/field/field_factory.hpp` | – / Omega_h | interface; XGC/MFEM adapters are FieldFactory, not FunctionSpace |
| `pcms/field/field_layout.h` | public | `include/pcms/field/field_layout.hpp` | – / Omega_h |  |
| `pcms/field/function_space.h` | public | `include/pcms/field/function_space.hpp` | – / Omega_h | interface + new space factories (or a sibling function_spaces.hpp) |
| `pcms/field/function_space/lagrange.h` | private | `src/pcms/field/function_space/lagrange.hpp` | Omega_h / – | impl; factory CreateLagrangeSpace |
| `pcms/field/function_space/mfem.hpp` | private | `src/pcms/backends/mfem/field_factory.hpp` | MFEM / Omega_h | impl; factory backend::mfem::CreateFieldFactory |
| `pcms/field/function_space/polynomial_reconstruction.hpp` | private | `src/pcms/field/function_space/polynomial_reconstruction.hpp` | Omega_h / – | impl; factory CreatePolynomialReconstructionSpace |
| `pcms/field/function_space/spline.h` | private | `src/pcms/field/function_space/spline.hpp` | – / Omega_h | impl; factory CreateSplineSpace |
| `pcms/field/function_space/xgc.h` | private | `src/pcms/backends/xgc/field_factory.hpp` | – / Omega_h | impl; factory backend::xgc::CreateFieldFactory |
| `pcms/field/gid_permutation.hpp` | private | `src/pcms/field/gid_permutation.hpp` | – / Omega_h |  |
| `pcms/field/layout/empty.h` | private | `src/pcms/field/layout/empty.hpp` | – / Omega_h | impl; factory CreateEmptyLayout |
| `pcms/field/layout/mesh_fields.h` | private | `src/pcms/backends/omega_h_meshfields/layout.hpp` | Omega_h / – | Omega_h + meshfields |
| `pcms/field/layout/mfem.hpp` | private | `src/pcms/backends/mfem/layout.hpp` | MFEM / Omega_h | OverlapMaskFromAttribute becomes a backend free function |
| `pcms/field/layout/omega_h_entity.h` | private | `src/pcms/backends/omega_h/entity_layout.hpp` | Omega_h / – |  |
| `pcms/field/layout/omega_h_lagrange.h` | private | `src/pcms/backends/omega_h/lagrange_layout.hpp` | Omega_h / – |  |
| `pcms/field/layout/point_cloud.h` | private | `src/pcms/backends/point_cloud/layout.hpp` | – / Omega_h |  |
| `pcms/field/layout/uniform_grid.h` | private | `src/pcms/backends/uniform_grid/layout.hpp` | – / Omega_h |  |
| `pcms/field/layout/xgc.h` | private | `src/pcms/backends/xgc/layout.hpp` | – / Omega_h |  |
| `pcms/field/out_of_bounds_policy.h` | public | `include/pcms/field/out_of_bounds_policy.hpp` | – / – |  |
| `pcms/field/point_evaluator.h` | public | `include/pcms/field/point_evaluator.hpp` | – / Omega_h |  |
| `pcms/field/point_evaluator_factory.hpp` | public | `include/pcms/field/point_evaluator_factory.hpp` | – / Omega_h |  |
| `pcms/field/transformed_source.hpp` | public | `include/pcms/field/transformed_source.hpp` | – / Omega_h |  |
| `pcms/field/uniform_grid_binary_field.h` | private | `src/pcms/backends/omega_h/uniform_grid_binary_field.hpp` | Omega_h / – | builds a uniform-grid field from an Omega_h mesh; Python binds it (question: public backend function?) |
| `pcms/field/value_view.hpp` | public | `include/pcms/field/value_view.hpp` | – / Omega_h |  |
| `pcms/localization/adj_search.hpp` | private | `src/pcms/backends/omega_h/adj_search.hpp` | – / Omega_h | Omega_h adapter (Phase 3.6) |
| `pcms/localization/localization_factory.h` | private | `src/pcms/localization/localization_factory.hpp` | – / Omega_h |  |
| `pcms/localization/mesh_localization.h` | private | `src/pcms/localization/mesh_localization.hpp` | Omega_h / – |  |
| `pcms/localization/mls_support_helpers.h` | private | `src/pcms/localization/mls_support_helpers.hpp` | Omega_h / – |  |
| `pcms/localization/point_cloud_localization.h` | private | `src/pcms/localization/point_cloud_localization.hpp` | Omega_h / – |  |
| `pcms/localization/point_search.h` | private | `src/pcms/localization/point_search.hpp` | Omega_h / – |  |
| `pcms/localization/queue_visited.hpp` | private | `src/pcms/backends/omega_h/queue_visited.hpp` | Omega_h / – | Omega_h adapter (Phase 3.6) |
| `pcms/pcms.h` | public | `include/pcms/pcms.hpp` | – / – | becomes umbrella header (Phase 4) |
| `pcms/transfer/bilinear_form_integrator.hpp` | private | `src/pcms/transfer/bilinear_form_integrator.hpp` | PETSc / – |  |
| `pcms/transfer/conservative_projection_solver.hpp` | private | `src/pcms/transfer/conservative_projection_solver.hpp` | Omega_h/PETSc / – |  |
| `pcms/transfer/coo_assembly_utils.hpp` | private | `src/pcms/transfer/coo_assembly_utils.hpp` | Omega_h/PETSc / – |  |
| `pcms/transfer/copy.h` | public | `include/pcms/transfer/copy.hpp` | – / Omega_h |  |
| `pcms/transfer/interpolation_base.h` | private | `src/pcms/transfer/interpolation_base.hpp` | Omega_h / – | legacy MLS interpolation used by the C API |
| `pcms/transfer/interpolation_helpers.h` | private | `src/pcms/transfer/interpolation_helpers.hpp` | Omega_h / – |  |
| `pcms/transfer/interpolator.h` | public | `include/pcms/transfer/interpolator.hpp` | – / Omega_h |  |
| `pcms/transfer/linear_form_integrator.hpp` | private | `src/pcms/transfer/linear_form_integrator.hpp` | PETSc / Omega_h |  |
| `pcms/transfer/linear_interpolant.hpp` | private | `src/pcms/transfer/linear_interpolant.hpp` | – / – | used by uniform-grid evaluator + a test |
| `pcms/transfer/mass_matrix_integrator.hpp` | private | `src/pcms/transfer/mass_matrix_integrator.hpp` | Omega_h/meshfields / – |  |
| `pcms/transfer/mass_matrix_type.hpp` | public | `include/pcms/transfer/mass_matrix_type.hpp` | – / – | enum in projection options (Python binds it) |
| `pcms/transfer/mesh_intersection.hpp` | private | `src/pcms/backends/omega_h/mesh_intersection.hpp` | Omega_h / – |  |
| `pcms/transfer/monte_carlo_sampling.hpp` | private | `src/pcms/transfer/monte_carlo_sampling.hpp` | – / – |  |
| `pcms/transfer/multidimarray.hpp` | private | `src/pcms/transfer/multidimarray.hpp` | – / – |  |
| `pcms/transfer/omega_h_conservative_projection.hpp` | private | `src/pcms/backends/omega_h/conservative_projection.hpp` | – / Omega_h | factory CreateConservativeProjection |
| `pcms/transfer/omega_h_control_variate_projection.hpp` | private | `src/pcms/backends/omega_h/control_variate_projection.hpp` | – / Omega_h | factory CreateControlVariateProjection |
| `pcms/transfer/omega_h_form_integrator_utils.hpp` | private | `src/pcms/backends/omega_h_meshfields/form_integrator_utils.hpp` | Omega_h/meshfields / – | Omega_h + meshfields |
| `pcms/transfer/omega_h_intersection_rhs_integrator.hpp` | private | `src/pcms/backends/omega_h/intersection_rhs_integrator.hpp` | – / Omega_h/PETSc |  |
| `pcms/transfer/omega_h_mass_integrator.hpp` | private | `src/pcms/backends/omega_h/mass_integrator.hpp` | – / Omega_h/PETSc |  |
| `pcms/transfer/omega_h_mc_rhs_integrator.hpp` | private | `src/pcms/backends/omega_h/mc_rhs_integrator.hpp` | – / Omega_h/PETSc |  |
| `pcms/transfer/petsc_utils.hpp` | private | `src/pcms/transfer/petsc_utils.hpp` | PETSc / – |  |
| `pcms/transfer/transfer_operator.hpp` | public | `include/pcms/transfer/transfer_operator.hpp` | – / – |  |
| `pcms/transient/timestepper.hpp` | public | `include/pcms/transient/timestepper.hpp` | – / – |  |
| `pcms/utility/array_mask.h` | private | `src/pcms/utility/array_mask.hpp` | – / Omega_h | no external user |
| `pcms/utility/arrays.h` | public | `include/pcms/utility/arrays.hpp` | Omega_h / – | View aliases; Omega_h overload MakeConstRank2View moves to the Omega_h backend |
| `pcms/utility/assert.h` | public | `include/pcms/utility/assert.hpp` | – / – | pcms_error, PCMS_ALWAYS_ASSERT |
| `pcms/utility/bounding_box.h` | private | `src/pcms/utility/bounding_box.hpp` | – / – | only a test uses it |
| `pcms/utility/common.h` | installed-detail | `include/pcms/utility/detail/common.hpp` | – / – | find_or_error helpers used by coupler.hpp templates |
| `pcms/utility/entity_types.h` | public | `include/pcms/utility/entity_types.hpp` | – / – | entity-dimension constants used by public interfaces |
| `pcms/utility/eqdsk.h` | public | `include/pcms/utility/eqdsk.hpp` | – / Omega_h | EQDSKData (MakeEQDSKField input) |
| `pcms/utility/inclusive_scan.h` | installed-detail | `include/pcms/utility/detail/inclusive_scan.hpp` | – / – | used by coupler.hpp templates |
| `pcms/utility/memory_spaces.h` | public | `include/pcms/utility/memory_spaces.hpp` | – / – |  |
| `pcms/utility/mesh_geometry.h` | private | `src/pcms/backends/omega_h/mesh_geometry.hpp` | Omega_h / – | Omega_h geometry kernels |
| `pcms/utility/mpi_type.h` | private | `src/pcms/utility/mpi_type.hpp` | – / – |  |
| `pcms/utility/omega_h_array_utils.h` | private | `src/pcms/backends/omega_h/array_utils.hpp` | Omega_h / – |  |
| `pcms/utility/print.h` | public | `include/pcms/utility/print.hpp` | – / – | used by examples/external-usage-example |
| `pcms/utility/profile.h` | installed-detail | `include/pcms/utility/detail/profile.hpp` | – / – | PCMS_FUNCTION_TIMER used inside public templates; perfstubs |
| `pcms/utility/types.h` | public | `include/pcms/utility/types.hpp` | – / – | LO, GO, Real, Type |
| `pcms/utility/uniform_grid.h` | public | `include/pcms/utility/uniform_grid.hpp` | Omega_h / – | UniformGrid; ClosestCellID takes a pcms vector (Phase 6); CreateUniformGridFromMesh moves to the Omega_h backend |
| `pcms/version.h` | public | `include/pcms/version.hpp` | – / – | generated |

Not installed today and stays private: `localization/localization_path_selection.h`.

Installed non-header files: `pcms/external/flibhpc/*` (SWIG `.i` files) and
`FindPETSc.cmake` are removal candidates in Phase 8.

## 6. Questions for the checkpoint

1. **Coupler templates and redev.** `coupler.hpp` (public) includes
   `field_communicator.hpp`, `field_exchange_planner.h`,
   `field_layout_communicator.h`, `global_communicator.h` and
   `field_serializer.h`, which are templates over the field value type and
   use redev in their bodies. Moving them to `detail/` keeps redev as a public
   dependency of the coupler. Removing it needs a hand-written change: a
   non-template (pimpl) boundary plus explicit instantiation for the closed set
   of value types in `FieldVariant`. That is behavior-preserving but is more
   than a move. Do it in this release, or list redev as a remaining public
   dependency of `pcms::coupler`?
2. **XGC and MFEM factory names.** `XGCFieldFactory` and `MFEMFieldFactory`
   derive from `FieldFactory`, not `FunctionSpace` (their fields are not
   evaluable). `decisions.md` says `CreateFunctionSpace`; the correct name is
   `backend::xgc::CreateFieldFactory` / `backend::mfem::CreateFieldFactory`
   returning `std::shared_ptr<FieldFactory>`. OK?
3. **`SimpleFieldData`** (`field/data/simple.h`): the generic `FieldData` over
   user-provided views, used by the C API and tests. Proposed public. Keep the
   class public, or hide it behind a `CreateFieldData(...)` factory?
4. **`spline_interpolator.hpp`** was installed on purpose by an earlier PR. The
   only user outside its sub-library is a test. Proposed private. Is there an
   external user (e.g. a PSPLINE replacement)?
5. **`uniform_grid_binary_field.h`** (builds a uniform-grid field from an
   Omega_h mesh) is bound in Python. Proposed: private implementation plus a
   free function in the Omega_h backend header. OK?
6. **`utility/print.h`, `bounding_box.h`, `array_mask.h`, `mpi_type.h`.**
   `print.h` stays public because the external usage example uses it; the
   others have no external user and become private.
