# PCMS: Parallel Coupler For Multimodel Simulations
 PCMS provides infrastructure to couple multiple simulation codes at scale on exascale systems. It uses scalable methodologies for transferring data between distributed applications. PCMS takes a hierarchical approach to field transfer using as much information as possible, such as mesh topology, to accommodate conservation and other physical constraints. In the context of limited information, alternative strategies, such as radial basis functions, are used to support field transfer. PCMS has been demonstrated on coupled fusion simulations that make use of volume coupling of high-dimensional data and complex coordinate systems.

## Run PCMS in a container

PCMS provides pre-built containers for x86\_64 and arm64. To use these containers, with docker you cand do the following:

```console
$ docker pull ghcr.io/scorec/pcms:develop
$ docker run --rm -it ghcr.io/scorec/pcms:develop bash
```

| Tag | Contents |
|---|---|
| `develop` | Latest build of `develop`, rebuilt at least weekly |
| `develop-YYYY-MM-DD` | build from specific date |

To run the Python examples, mount a checkout of this repository:

```console
$ git clone https://github.com/SCOREC/pcms.git
$ docker run --rm -v $PWD/pcms:/pcms -w /pcms ghcr.io/scorec/pcms:develop \
    python3 examples/python-api-example/example_conservative_transfer.py
```

If you run into issues due to containers running as root, you may need to use the following environment variables for OpenMPI to work.
As with any container running as root, exercise caution.
`mpirun`:

```console
$ docker run --rm -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ghcr.io/scorec/pcms:develop mpirun -n 2 python3 -c "import pcms"
```

More details about the container build process can be found here
[scorec-spack-packages](https://github.com/SCOREC/scorec-spack-packages#pcms-container-image).

## Build Instructions

There are two supported ways to build PCMS from source:

- [Spack](#build-with-spack) works on any Linux or macOS system.
- On [SCOREC systems](#build-on-scorec-systems), use the module-based scripts.

## Build with Spack

The environment in [`spack/spack.yaml`](spack/spack.yaml) builds PCMS from
your checkout. It fetches the PCMS
recipes from [scorec-spack-packages](https://github.com/SCOREC/scorec-spack-packages),
pins the version of the Spack package repository, and installs prebuilt
dependencies from the SCOREC build cache when they match your system
(Ubuntu 24.04 on x86\_64 or arm64). Elsewhere, Spack builds the dependencies
from source the first time.

1. Install Spack 1.2 or newer:
   ```console
   $ git clone --depth=2 -b releases/v1.2 https://github.com/spack/spack.git
   $ . spack/share/spack/setup-env.sh
   ```
2. Spack needs C, C++, and Fortran compilers. On Ubuntu, install them with
   `sudo apt install gcc g++ gfortran`. Spack finds them on the first build;
   `spack compiler list` shows what it found.
3. Build and install PCMS:
   ```console
   $ git clone https://github.com/SCOREC/pcms.git
   $ cd pcms
   $ spack env activate ./spack
   $ spack install
   ```

This environment has PCMS in *develop* mode after changing the source, you can run
`spack install` again to rebuild it. With the environment active, its view is on your `PATH`, so
`python3 -c "import pcms"` uses the Python bindings.

To build a different configuration, edit the spec in `spack/spack.yaml` and
run `spack concretize -f` before `spack install`. For example:

- `pcms@develop` builds without the optional features (Python, PETSc, NetCDF).
- `+cuda cuda_arch=80` (or `+rocm amdgpu_target=gfx90a`) builds for GPUs.
  Also add `~tests` if the build machine has no GPU, because the tests are
  run during the build.

`spack info pcms` lists available variants.

## Build on SCOREC systems

On the SCOREC RHEL9 machines, PCMS and its dependencies are built with the
system modules and the scripts in
[laces-lab/software-dependencies](https://github.com/laces-lab/software-dependencies).
Follow the instructions there.

## Python API

Point-evaluator construction in Python uses `EvaluationRequest`, for example:

```python
request = pcms.EvaluationRequest.from_coordinates(eval_coords)
evaluator = source_space.create_point_evaluator(request)

optimized_request = pcms.EvaluationRequest.from_function_space(target_space)
evaluator = source_space.create_point_evaluator(optimized_request)
```

## Creating Archive for Release
`git archive --format=tar.gz -o /tmp/pcms.tar.gz --prefix=pcms/ develop`

## Running clang format on SCOREC
```console
module use /opt/scorec/spack/rhel9/v0222_2/lmod/linux-rhel9-x86_64/Core/
module load clang-format
module load llvm
find . -regex '.*\.\(cpp\|cxx\|cc\|c\|hpp\|h\)' -exec clang-format -style=file -i {} \;
```
