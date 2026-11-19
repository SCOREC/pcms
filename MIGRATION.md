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
- **Without Omega_h:** `PCMS_ENABLE_OMEGA_H=OFF` (with
  `PCMS_ENABLE_SERVER=OFF`, which forces Omega_h on, and the coupler, XGC
  and MeshFields off) builds the core: uniform-grid and point-cloud
  geometries and their spaces, QueryPoints, MLS. Mesh geometries, the mesh
  search, the projections and the C and Fortran APIs need Omega_h: they are
  switched off (with a status message) when it is off.
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
- `tools/reorg/header-map.txt` lists every old path and its replacement;
  headers marked "private" are no longer installed and their classes must be
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

(`` is `pcms::`; all other names are in `pcms::`.) Factories
for a combination PCMS does not support (e.g. a Lagrange space on a
point-cloud geometry) throw `pcms_error`. The projection factories throw in a
build without PETSc and MeshFields.

- A geometry must be owned by a `std::shared_ptr` (the factories return one);
  spaces and BoundPoints hold it, and what it refers to (the Omega_h mesh)
  must outlive it. Creating a geometry builds the mesh's containing-element
  search: create one per mesh and create every space on it, so that the
  spaces share the search and accept each other's BoundPoints.
- There are no discretization factories (`Discretization` is removed; layouts
  hold their geometry), and a geometry owns its
  global-id tag and ownership mask (`OmegaHGeometryOptions`).
  `OmegaHGeometry::GetMesh()` returns the mesh of an Omega_h geometry; the
  geometry kinds are public and their factories return them.
- The coordinate system belongs to the geometry; the space factories no
  longer take one.
- `MLSOptions` keeps `degree`, `lambda`, `tol`, `decay_factor` and `basis`;
  `radius`, `min_req_supports` and `adapt_radius` are `NeighbourhoodOptions`
  (`radius`, `min_supports`, `adapt_radius`) of the point-cloud geometry,
  whose support search is the geometry's.
- `MonteCarloSampling` is removed: the control-variate projection samples
  uniformly at random.

`MakeEQDSKField` is unchanged except that `EQDSKField::space` is a
`std::shared_ptr<FunctionSpace>`.

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
| `source.CreatePointEvaluator<T>(EvaluationRequest::FromFunctionSpace(space, policy))` | `source.CreatePointEvaluator<T>(source.GetGeometry()->Bind(space.GetLayout()->GetDOFHolderCoordinates(), policy))`, or an `Interpolator` |
| `TransformedSource(source, map)` | `geometry->Bind(coords, policy, map)`: the points are in the map's source coordinate system, are located at their images, and the evaluators re-express the values through the map |
| (none) | `source.CreatePointEvaluator<T>(bound, DifferentialOperator::Gradient)`; `Gradient(source)`, `Divergence(source)`, `Curl(source)` are sources whose evaluators write the operator applied to the field (`pcms/field/gradient_source.hpp`); `Supports(op)` tells which a source implements |
| (none) | `geometry->Bind(coords, OwnerHints{...}, policy, map)`: points whose owners are known (generated in the elements of the geometry) skip the search; `geometry->Rebind(previous, coords)` recycles the storage of a BoundPoints nothing else holds |
| (none) | `QueryPoints query(geometry, policy, map); query.Locate(coords); query.Evaluate(function, out, op)`: evaluation at points that move between calls; `Evaluate(function, points, out, op, policy, map)` (`pcms/field/query_points.hpp`) does one evaluation |
| (none) | `geometry->Neighbourhood(coords, NeighbourhoodOptions, owners)`: the supports of reconstruction spaces; `geometry->Differential(bound)`: the Jacobian at the bound points |

A `PointEvaluatorFactory` implemented outside PCMS implements
`GetGeometry()`, `Supports(op)` and `CreatePointEvaluatorImpl(type, bound,
op)` (`pcms/field/point_evaluator_factory.hpp`; deriving from
`TypedPointEvaluatorFactory<T, Bound>` and implementing `CreateEvaluator`
gets the value-type and operator checks), and a `Geometry` implemented
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
are removed; the replacements are the bind above, `Interpolator` (also for
reconstruction spaces) and `Geometry::Neighbourhood`.

## Namespaces

- The public API stays in `pcms::`.
- Backend code is in `pcms::<name>` (`omega_h`, `mfem`, `xgc`,
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
| `space.mesh()` | `pcms.omega_h.get_mesh(space.get_geometry())` |
| `pcms.create_uniform_grid_binary_field(mesh, divs)` | `pcms.omega_h.create_uniform_grid_binary_field(mesh, divs)` |
| `pcms.create_uniform_grid_from_mesh(mesh, divs)` | `pcms.omega_h.create_uniform_grid_from_mesh(mesh, divs)` |
| `pcms.map_entity_field_to_vertices_{average,mls}(mesh, ...)` | `pcms.omega_h.map_entity_field_to_vertices_{average,mls}(mesh, ...)` |
| `pcms.compute_entity_centroids(mesh, dim)`, `pcms.compute_entity_centroids_few(mesh, dim)` | `pcms.omega_h.compute_entity_centroids(mesh, dim)` |

A geometry keeps its mesh alive and a space keeps its geometry alive;
`pcms.Geometry.bind`, `pcms.QueryPoints`, `pcms.gradient` and
`pcms.evaluate` are the Python forms of the bind and the evaluation above.
The MLSOptions fields `radius`, `min_req_supports` and `adapt_radius` are
`pcms.NeighbourhoodOptions` fields (`radius`, `min_supports`,
`adapt_radius`) of the point-cloud geometry.

`tools/reorg/python_api_names.py` (run by `migrate-pcms.sh`) rewrites these
calls; it reports the MLSOptions fields that moved.

## Porting open pull requests

A pull request written against 0.3 is ported by merging the reorganization in
two parts, before and after its namespace commit ("refactor: move internal
code into the sub-library and backend namespaces", `$N` below), and applying
that commit's spec to the merged tree in between, so that the pull request's
code moves to the new namespaces with the rest. This needs a pull request
that builds (in the Full configuration); for one that does not, see the end
of this section.

`<develop>` is the reorganized develop. The first part does not contain
`tools/migrate-pcms.sh` yet, so use it from a separate checkout:

```sh
git worktree add ../pcms-0.4 <develop>
N=$(git log --format=%H -1 --grep='move internal code into the sub-library and backend namespaces' <develop>)
```

1. Bring the branch up to date with the last `develop` commit before the
   reorganization (the parent of the first reorganization commit, "build:
   require CMake 3.23") and resolve those conflicts as usual.
2. Merge the reorganization up to the namespace commit. In conflicting C++,
   Python and Markdown files keep the branch's side of each conflicting hunk;
   git applies the branch's edits to renamed and moved files:

   ```sh
   git merge -X find-renames $N~1
   for f in $(git diff --name-only --diff-filter=U -- '*.cpp' '*.hpp' '*.h' '*.py' '*.md'); do
     # only files changed on both sides; a file one side deleted or moved
     # (modify/delete) is left for the by-hand step
     [ "$(git ls-files -u -- "$f" | wc -l)" -eq 3 ] || continue
     git show :2:$f > $f.ours; git show :1:$f > $f.base; git show :3:$f > $f.theirs
     git merge-file --ours $f.ours $f.base $f.theirs && mv $f.ours $f && rm $f.base $f.theirs && git add $f
   done
   ```

   Resolve the remaining conflicts by hand (CMake lists and modify/delete
   conflicts, see below). Then review what the branch changes relative to the
   reorganization, `git diff $N~1 -- <the branch's files>`: a kept hunk that
   only edited code the reorganization removed or replaced brings that code
   back (for example, comment edits on `OverlapMask`'s Omega_h constructors,
   which became `pcms::omega_h::OverlapMaskFromArray`); take the
   reorganization's side there. Commit the merge. Do not use
   `git merge -X ours` for the whole merge: in CMake files it silently drops
   the reorganization's changes.
3. Migrate includes, macros, CMake options and Python in the branch's files,
   but not C++ names (step 4 moves those), and make it build:

   ```sh
   ../pcms-0.4/tools/migrate-pcms.sh --in-tree --no-names $(git diff --name-only --diff-filter=d HEAD^2)
   ```

   Commit.
4. Apply the namespace spec (`tools/reorg/refactor/README.md` lists what the
   tool needs). A new file in a sub-library or backend directory that should
   move with it is added to `tools/reorg/refactor/namespaces.spec` first.

   ```sh
   tools/reorg/refactor/build.sh <bindir>
   cmake -S . -B <build> -DCMAKE_EXPORT_COMPILE_COMMANDS=ON <every PCMS_ENABLE_* option ON> ...
   PCMS_REFACTOR=<bindir>/pcms-refactor \
     tools/reorg/refactor/move_namespaces.sh . <build> tools/reorg/refactor/namespaces.spec
   clang-format -i $(git diff --name-only -- '*.cpp' '*.hpp' '*.h')   # clang-format 18.1.3
   ```

   Commit. Outside the branch's files the tree is now the namespace commit's
   (`git diff --name-only $N` lists only the branch's files).
5. Record the namespace commit as merged and merge the rest, keeping the
   branch's side of conflicting hunks as in step 2; then migrate the
   remaining names and format the branch's files:

   ```sh
   git merge -s ours $N
   git merge <develop>
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
  evaluators, the coupler's communicators and planners, the MLS helpers); the
  symbol map lists the 0.3 names as `cpp-private`.
- Uses of APIs that changed shape: `migrate-pcms.sh` reports them (the
  `manual` entries of `tools/reorg/symbol-map.txt`).
- Code that calls the point search directly (pull requests #345, #416 and
  #417): the search `Result` now reports the entity the point lies on
  (`dimensionality`, `entity`; `OutsideId`/`IsOutside` for a point outside)
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

A pull request that does not build before the reorganization is merged in one
step instead: step 2 with `<develop>` in place of `$N~1`, then step 5 without
the two merges. `migrate-pcms.sh --in-tree` then writes the moved names in
full (`pcms::localization::GridPointSearch2D`) where the tool would have
written them as briefly as lookup allows.

## C and Fortran APIs

Unchanged, except that the C headers are installed in `include/pcms/capi/`
and the installed include directory is `include` (include them as
`pcms/capi/<name>.h`, as the build tree already required).
