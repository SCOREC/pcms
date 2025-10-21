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


## Maintain CI/CD

Current CI/CD contains:
- Cmake build and test using GitHub Actions triggered on pushes and pull requests to the `develop`,`main` branch. The workflow file is [`.github/workflows/cmake-test.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/cmake-test.yml).
- Format check using clang-format and cmake-format. Workflow file is [`.github/workflows/format-check.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/format-check.yml).
- Static code analysis using clang-tidy. Combined in the main cmake workflow at [`.github/workflows/cmake-test.yml`](https://github.com/SCOREC/pcms/blob/develop/.github/workflows/cmake-test.yml).
- Weekly globus build/test in supercomputer Perlmutter. Workflow file is in [`.github/workflows/perlmutter/`](https://github.com/SCOREC/pcms/tree/develop/.github/workflows/perlmutter) directory.
- Nightly build in [CDash dashboard](https://my.cdash.org/index.php?project=SCOREC). The build script is in [`cdash/`](https://github.com/SCOREC/pcms/tree/develop/cdash) directory. The script is triggered by a cron job on a SCOREC workstation daily at early morning.

To maintain the CI/CD pipeline, ensure that any changes to the build or test processes are reflected in those workflow file. Regularly review the workflow for updates to dependencies or tools used in the pipeline. When adding new steps to the CI/CD process, ensure they are properly integrated and tested to avoid disruptions. Try to avoid third party actions for security and maintenance reasons.
