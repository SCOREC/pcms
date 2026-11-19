## Creating a Release

0. select the new version number following https://semver.org/:

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
>    - MAJOR version when you make incompatible API changes
>    - MINOR version when you add functionality in a backward compatible manner
>    - PATCH version when you make backward compatible bug fixes


1. create a github issue documenting significant release changes; review the commit log and closed issues to find them

```
This issue is to document functionality and features added to pcms since the #.#.# release (SHA1 of prior release):

New functionality or feature support:

- <feature> (SHA1,issueNumber)
- ...

Bug Fixes:

- <feature> (SHA1,issueNumber)
- ...

Other Updates and Improvements:

- <feature> (SHA1,issueNumber)
- ...
```

2. apply the issue/PR label 'v#.#.#' to significant issues and PR that are part of the release
3. increase the pcms version # in [CMakeLists.txt](https://github.com/Sichao25/pcms/blob/yus/format_check_test/CMakeLists.txt#L4) in the `develop` branch
4. commit; include the issue # in the commit message

```
pcms version #.#.#                                                                                                                                                        
see issue #<###>
```

5. push
6. create the tag `git tag -a v#.#.# -m "pcms version #.#.#"`
7. push the tag `git push origin v#.#.#`

## Branch Management

The `develop` branch is the main development branch where all new features and bug fixes are merged. The `master` branch is the stable release branch that reflects the latest released version of pcms. All changes should be made in the `develop` branch and then merged into `master` when ready for release. 

It is generally recommended to use forks for development to keep the main repository clean and to facilitate code review. Contributors can create a fork of the repository, make their changes in a new branch, and then submit a pull request to the `develop` branch of the main repository. Ensure that the pull request includes a clear description of the changes made and references any relevant issues.

When creating a new feature or fixing a bug in main repo, create a new branch from `develop` with a descriptive name (e.g., `feature/new-feature`, `bugfix/issue-123`).

## Maintain CI/CD

Current CI/CD contains:
- Cmake build and test using GitHub Actions triggered on pushes and pull requests to the `develop`,`main` branch. The workflow file is [`.github/workflows/cmake-test.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/cmake-test.yml).
- Format check using clang-format and cmake-format. Workflow file is [`.github/workflows/format-check.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/format-check.yml).
- Static code analysis using clang-tidy. Combined in the main cmake workflow at [`.github/workflows/cmake-test.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/cmake-test.yml).
- Weekly globus build/test in supercomputer Perlmutter. Workflow file is in [`.github/workflows/perlmutter/`](https://github.com/SCOREC/pcms/tree/develop/.github/workflows/perlmutter) directory.
- Nightly build in [CDash dashboard](https://my.cdash.org/index.php?project=SCOREC). The build script is in [`cdash/`](https://github.com/SCOREC/pcms/tree/develop/cdash) directory. The script is triggered by a cron job on a SCOREC workstation daily at early morning.

To maintain the CI/CD pipeline, ensure that any changes to the build or test processes are reflected in those workflow file. Regularly review the workflow for updates to dependencies or tools used in the pipeline. When adding new steps to the CI/CD process, ensure they are properly integrated and tested to avoid disruptions. Try to avoid third party actions for security and maintenance reasons.

## Conventions

How the public API is organized; keep new code to these rules.

1. **Minimal public API.** If in doubt, make a name private: it can be made
   public later, but not taken back. Every public name is needed by an
   example, a factory signature or a Python binding.
2. **Python mirrors C++.** The module path is the namespace path
   (`pcms::omega_h` is `pcms.omega_h`), and there are no Python-only
   functions.
3. **Backends.** Each backend has one flat namespace (`pcms::omega_h`,
   `pcms::mfem`, `pcms::xgc`, `pcms::omega_h_meshfields`,
   `pcms::uniform_grid`, `pcms::point_cloud`); its headers are in
   `include/pcms/backends/<name>/` and its sources in
   `src/pcms/backends/<name>/`. A backend folder holds everything that
   supports the backend: geometry, layouts, evaluators and special-purpose
   search. Code shared by several backends lives in core (`UniformGrid<dim>`
   is public; the MLS support search is private, in `src/pcms/localization/`).
4. **Internal namespaces.** Inside the library there is one namespace per
   library (`pcms::utility`, `pcms::field`, `pcms::transfer`,
   `pcms::coupler`, ...) and `pcms::<backend>` for backends. Anything
   installed that is not public is in a header below a `detail/` folder and
   in a `detail` namespace: `include/pcms/<lib>/detail/` with `pcms::detail`,
   `include/pcms/backends/<b>/detail/` with `pcms::<b>::detail`. Public
   headers have no `detail` blocks of their own; they include the `detail/`
   header. `tools/check-public-includes.py` checks this.
5. **Factories return the base type.** Concrete classes of polymorphic
   families (geometries, transfers, coordinate systems, coordinate maps,
   timesteppers) live in `src/` and are never installed. Factories return
   `std::shared_ptr<const Geometry>`, `std::unique_ptr<TransferOperator<T>>`,
   `std::shared_ptr<const CoordinateSystem>`,
   `std::shared_ptr<const CoordinateMap>` or `std::unique_ptr<Timestepper>`.
   A function that needs a particular kind checks it at runtime and throws
   `pcms_error`. Non-polymorphic value types (`Field`, `Function`,
   `QueryPoints`, `OverlapMask`, `UniformGrid`, options structs) stay public
   classes. Every factory validates its inputs.
6. **No pimpl.** Keeping the concrete classes in `src/` makes pimpl and
   forward-declaration workarounds unnecessary.
7. **Naming.** `Create*` builds a PCMS object, `Make*` wraps existing memory
   in a non-owning view, `Read*` loads from a file.
8. **Ownership is explicit.** An object that needs another to stay alive
   takes and stores a `std::shared_ptr` to it. Python bindings never use
   `py::keep_alive`.
9. **`::mfem`.** MFEM's namespace collides with `pcms::mfem`: inside
   `namespace pcms`, write `::mfem::`.

**Extension points.** External code may implement `Geometry`,
`PointEvaluatorFactory`, `BoundPoints` and `FunctionSpace` (and `FieldLayout`,
`FieldData`, `CoordinateMap`, `Timestepper`); their protected and extension
APIs are public and stay so.
