# Migrating from PCMS 0.3 to 0.4

PCMS 0.4 reorganizes the library so that downstream codes migrate once: the
public API is now a small set of interfaces and factory functions, everything
else is internal. Behavior and numerical results do not change.

Most of the mechanical work is done by `tools/migrate-pcms.sh`:

```sh
tools/migrate-pcms.sh --dry-run path/to/your/code   # show what would change
tools/migrate-pcms.sh path/to/your/code             # apply it
```

It rewrites includes, renamed names, macros, CMake options and targets, and
Python code, and lists every place that needs a manual change (with the
replacement). The rest of this document explains those changes.

## Building and linking

- **Requirements:** CMake 3.23 or newer.
- **Link `pcms::pcms`.** It links every PCMS library that was built and is the
  supported way to use PCMS. `pcms::core` still works in 0.4 as a deprecated
  alias of `pcms::pcms` (CMake prints a deprecation warning) and will be
  removed. The other targets are `pcms::utility`, `pcms::discretization`,
  `pcms::localization`, `pcms::field`, `pcms::transfer`, `pcms::coupler`,
  `pcms::transient`, `pcms::capi`, `pcms::fortranapi`.
- **Renamed options** (the old names still work in 0.4, with a deprecation
  warning): `PCMS_ENABLE_Python` → `PCMS_ENABLE_PYTHON`,
  `PCMS_ENABLE_Fortran` → `PCMS_ENABLE_FORTRAN`,
  `PETSC_LINK_STATIC` → `PCMS_PETSC_LINK_STATIC`.
- **New option** `PCMS_ENABLE_COUPLER` (default `ON`): with `OFF`, PCMS builds
  without redev and ADIOS2 (no coupler, no C/Fortran client API, no Python).
- **Finding PCMS:** the package files are in `<prefix>/lib/cmake/pcms`. The
  package no longer bakes in the paths of its dependencies: make them findable
  as usual (`CMAKE_PREFIX_PATH`, `<Package>_ROOT` or `<Package>_DIR`) for
  Kokkos, KokkosKernels, perfstubs, Omega_h, and, if PCMS was built with them,
  redev, ADIOS2, MeshFields, MFEM and PETSc.
- **Versions:** `find_package(pcms 0.4)` accepts only 0.4.x
  (`SameMinorVersion`); shared libraries have SOVERSION `0.4`.

## Headers

- All C++ headers end in `.hpp` (the C API headers in `pcms/capi/` stay `.h`).
- Only the public headers are installed, under `include/pcms/<library>/`.
  `pcms/pcms.hpp` includes the whole public API; backend headers are included
  separately (`pcms/backends/omega_h/omega_h.hpp`,
  `pcms/backends/mfem/mfem.hpp`, `pcms/backends/xgc/xgc.hpp`).
- The top-level `pcms.h` is replaced by `pcms/pcms.hpp`.
- `tools/reorg/header-map.txt` lists every old path and its replacement;
  headers marked "private" are no longer installed and their classes must be
  replaced by the public API below.
- Core public headers no longer include Omega_h, MFEM, MeshFields, PETSc or
  ADIOS2. Code that got Omega_h declarations through a PCMS header must
  include the Omega_h headers it uses.

## Creating spaces, fields and transfers

The concrete classes are internal. Use the factories, which return the
interfaces `Discretization`, `FunctionSpace`, `FieldFactory`, `FieldLayout`
and `TransferOperator`:

| 0.3 | 0.4 |
|---|---|
| `LagrangeFunctionSpace::FromMesh(mesh, order, n, csys, gid_name, backend, layout_name)` | `CreateLagrangeSpace(backend::omega_h::CreateDiscretization(mesh, gid_name), order, n, csys, {implementation, layout_name})` |
| `... FromMesh(mesh, order, n, csys, owned_mask, ...)` | `CreateLagrangeSpace(backend::omega_h::CreateDiscretization(mesh, owned_mask, gid_name), ...)` |
| `LagrangeFunctionSpace::FromUniformGrid(grid, n, csys, order)` | `CreateLagrangeSpace(CreateUniformGridDiscretization(grid), order, n, csys)` |
| `LagrangeFunctionSpace::Backend::{MeshFields,OmegaH}` | `LagrangeImplementation::{MeshFields,OmegaH}` (in `LagrangeOptions`) |
| `SplineFunctionSpace::FromUniformGrid(grid, csys)` | `CreateSplineSpace(CreateUniformGridDiscretization(grid), csys)` |
| `PolynomialReconstructionFunctionSpace::Create(points, opts)` | `CreatePolynomialReconstructionSpace(points, opts)` |
| `PolynomialReconstructionFunctionSpace::FromMesh(mesh, dim, csys, opts)` | `CreatePolynomialReconstructionSpace(backend::omega_h::CreateDiscretization(mesh), dim, csys, opts)` |
| `EmptyFieldLayout(csys, n)` | `CreateEmptyLayout(csys, n)` |
| `OmegaHConservativeProjection(src, tgt, mass)` | `CreateConservativeProjection(src, tgt, {mass})` |
| `OmegaHControlVariateProjection(src, tgt, samples, sampling, seed, mass)` | `CreateControlVariateProjection(src, tgt, {samples, sampling, seed, mass})` |
| `XGCFieldFactory(rc, n, layout_name)` | `backend::xgc::CreateFieldFactory(rc, n, layout_name)` |
| `factory.CreateField<T>(name, make_unique<XGCFieldData<T>>(factory.GetXGCLayout(), basis, data, comm))` | `backend::xgc::CreateField<T>(factory, name, data, comm, basis)` |
| `MFEMFieldFactory(pmesh, pfes, gf, csys, layout_name)` | `backend::mfem::CreateFieldFactory(pmesh, pfes, gf, csys, layout_name)` |
| `MFEMLayout::OverlapMaskFromAttribute(pmesh, attr)` | `backend::mfem::OverlapMaskFromAttribute(pmesh, attr)` |
| `OverlapMask(size, omega_h_mask)` | `OverlapMask(size, backend::omega_h::OverlapMaskFromArray(omega_h_mask))` |
| `MakeConstRank2View(omega_h_array, dim)` | `backend::omega_h::MakeConstRank2View(omega_h_array, dim)` |
| `CreateUniformGridFromMesh<dim>(mesh, ...)` | `backend::omega_h::CreateUniformGridFromMesh<dim>(mesh, ...)` |
| `CreateUniformGridBinaryField<2>(mesh, ...).second` | `backend::omega_h::CreateUniformGridBinaryField(mesh, ...)` |

(`backend::` is `pcms::backend::`; all other names are in `pcms::`.) Factories
for a combination PCMS does not support (e.g. a Lagrange space on an XGC
discretization) throw `pcms_error`. The projection factories throw in a build
without PETSc and MeshFields.

`MakeEQDSKField` is unchanged except that `EQDSKField::space` is a
`std::shared_ptr<FunctionSpace>`.

`UniformGrid::ClosestCellID` and `IsPointInBounds` accept any point type
indexable with `[i]` (`Kokkos::Array`, a subview, `Omega_h::Vector`).

## Namespaces

- The public API stays in `pcms::`.
- Backend code is in `pcms::backend::<name>` (`omega_h`, `mfem`, `xgc`,
  `point_cloud`, `uniform_grid`, `omega_h_meshfields`).
- Internal code is in `pcms::utility`, `pcms::discretization`,
  `pcms::localization`, `pcms::field`, `pcms::transfer` and `pcms::coupler`.
  `tools/reorg/symbol-map.txt` lists every name that moved.

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

The Python module mirrors the C++ factories; the class-based constructors are
removed:

| 0.3 | 0.4 |
|---|---|
| `pcms.LagrangeFunctionSpace.from_mesh(mesh, order, n, csys, backend=...)` | `pcms.create_lagrange_space(pcms.backend.omega_h.create_discretization(mesh), order, n, csys, implementation=pcms.LagrangeImplementation.X)` |
| `pcms.LagrangeFunctionSpace.from_uniform_grid(grid, n, csys, order)` | `pcms.create_lagrange_space(pcms.create_uniform_grid_discretization(grid), order, n, csys)` |
| `pcms.PolynomialReconstructionFunctionSpace.create(coords, opts)` | `pcms.create_polynomial_reconstruction_space(coords, opts)` |
| `pcms.PolynomialReconstructionFunctionSpace.from_mesh(mesh, dim, csys, opts)` | `pcms.create_polynomial_reconstruction_space(pcms.backend.omega_h.create_discretization(mesh), dim, csys, opts)` |
| `pcms.OmegaHConservativeProjection(src, tgt, ...)` | `pcms.create_conservative_projection(src, tgt, ...)` |
| `pcms.OmegaHControlVariateProjection(src, tgt, ...)` | `pcms.create_control_variate_projection(src, tgt, ...)` |
| `space.mesh()` | `pcms.backend.omega_h.get_mesh(space.get_discretization())` |
| `pcms.create_uniform_grid_binary_field(mesh, divs)` | `pcms.backend.omega_h.create_uniform_grid_binary_field(mesh, divs)` |
| `pcms.create_uniform_grid_from_mesh(mesh, divs)` | `pcms.backend.omega_h.create_uniform_grid_from_mesh(mesh, divs)` |
| `pcms.map_entity_field_to_vertices_{average,mls}(mesh, ...)` | `pcms.backend.omega_h.map_entity_field_to_vertices_{average,mls}(mesh, ...)` |
| `pcms.compute_entity_centroids(mesh, dim)`, `pcms.compute_entity_centroids_few(mesh, dim)` | `pcms.backend.omega_h.compute_entity_centroids(mesh, dim)` |

`tools/reorg/python_api_names.py` (run by `migrate-pcms.sh`) rewrites these
calls.

## Porting open pull requests

A pull request written against 0.3 is ported with a merge and the scripts
that made the reorganization (all are idempotent):

1. Bring the branch up to date with the last `develop` commit before the
   reorganization (c33f816c) and resolve those conflicts as usual.
2. Merge the reorganization. In conflicting C++, Python and Markdown files
   keep the branch's side of each conflicting hunk (steps 3 and 4 re-apply
   the reorganization's changes to them); git applies the branch's edits to
   renamed and moved files:

   ```sh
   git merge -X find-renames <develop with the reorganization>
   for f in $(git diff --name-only --diff-filter=U -- '*.cpp' '*.hpp' '*.h' '*.py' '*.md'); do
     # only files changed on both sides; a file one side deleted or moved
     # (modify/delete) is left for the by-hand step
     [ "$(git ls-files -u -- "$f" | wc -l)" -eq 3 ] || continue
     git show :2:$f > $f.ours; git show :1:$f > $f.base; git show :3:$f > $f.theirs
     git merge-file --ours $f.ours $f.base $f.theirs && mv $f.ours $f && rm $f.base $f.theirs && git add $f
   done
   ```

   Resolve the remaining conflicts by hand (CMake lists and modify/delete
   conflicts, see below), then
   commit the merge. Do not use `git merge -X ours` for the whole merge: in
   CMake files it silently drops the reorganization's changes.

3. After committing the merge, re-apply the namespace changes and keep them
   only in the files that carry the branch's changes (those that differ from
   the merged `develop`, `HEAD^2`):

   ```sh
   for b in xgc mfem omega_h_meshfields point_cloud uniform_grid omega_h; do
     tools/reorg/backend_namespace.py $b tools/reorg/namespaces/$b.txt; done
   for l in utility discretization localization field transfer coupler; do
     tools/reorg/sublib_namespace.py $l tools/reorg/namespaces/sublib-$l.txt \
       tools/reorg/namespaces/sublib-$l-files.txt; done
   git diff --name-only | grep -vxF -f <(git diff --name-only HEAD^2 HEAD) | xargs -r git checkout --
   ```

4. Migrate includes and names in the branch's files, and format them:

   ```sh
   tools/migrate-pcms.sh --in-tree $(git diff --name-only --diff-filter=d HEAD^2)
   clang-format -i $(git diff --name-only --diff-filter=d HEAD^2 | grep -E '\.(cpp|hpp|h)$')
   ```

   `--in-tree` rewrites private headers and internal names to their new
   locations instead of reporting them.

What is left by hand:

- CMake source lists: keep the reorganization's list and add the branch's
  new files.
- New headers: name them `.hpp`, decide public (`include/pcms/<lib>/`) or
  private (next to the sources), and list them accordingly.
- Files the reorganization rewrote beyond rename detection (git reports a
  modify/delete conflict): re-apply the branch's hunks to the new file, or
  delete it if the branch deleted it.
- New or renamed files in a sub-library or backend directory: before step 3,
  add them to `tools/reorg/namespaces/sublib-<lib>-files.txt` (or the
  backend's list), remove the names they replace, and add the namespace-scope
  names they declare to `sublib-<lib>.txt`, so that step 3 moves them into
  `pcms::<lib>` with the rest of the library.

## C and Fortran APIs

Unchanged, except that the C headers are installed in `include/pcms/capi/`
and the installed include directory is `include` (include them as
`pcms/capi/<name>.h`, as the build tree already required).
