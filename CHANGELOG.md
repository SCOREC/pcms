# Changelog

## 0.4.0

A one-time breaking reorganization so that downstream codes migrate once.
Numerical behavior is unchanged. See [MIGRATION.md](MIGRATION.md) and
`tools/migrate-pcms.sh`; the decisions behind the changes are in
`tools/reorg/decisions.md`.

### Changed (breaking)

- Public API: only interfaces and factory functions are public
  (`CreateLagrangeSpace`, `CreateSplineSpace`,
  `CreatePolynomialReconstructionSpace`, `CreateUniformGridDiscretization`,
  `CreateEmptyDiscretization`, `CreateEmptyLayout`,
  `CreateConservativeProjection`, `CreateControlVariateProjection`, and the
  backend entry points in `pcms/backends/<name>/<name>.hpp`). The concrete
  spaces, layouts, discretizations and projections are internal.
- Headers: all C++ headers are `.hpp`; only public headers are installed, in
  `include/pcms/<library>/`; `pcms.h` is replaced by the umbrella
  `pcms/pcms.hpp`; core public headers no longer include Omega_h, MFEM,
  MeshFields, PETSc or ADIOS2.
- Namespaces: backend code is in `pcms::backend::<name>`, internal code in
  `pcms::<library>`; the public API stays in `pcms::`.
- Omega_h types left the core public signatures (`MakeConstRank2View`,
  `CreateUniformGridFromMesh` and the Omega_h `OverlapMask` constructors moved
  to the Omega_h backend; `UniformGrid` point queries take any indexable
  point).
- `pcms::Mode` is a PCMS enum instead of an alias of `redev::Mode`.
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

### Added

- `PCMS_ENABLE_COUPLER` to build without redev and ADIOS2.
- Factory functions and backend headers listed above; `FunctionSpace`
  discretizations in Python (`pcms.backend.omega_h`).
- Tests that every public header compiles on its own and that core public
  headers include no backend library; CI builds and runs
  `examples/external-usage-example` against the install tree.
- `tools/migrate-pcms.sh` with `tools/reorg/header-map.txt` and
  `tools/reorg/symbol-map.txt`.
