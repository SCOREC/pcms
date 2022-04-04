# WDMApp Coupling

Adios2-based coupler for XGC and GENE

## Dependencies

- CMake 3.19+
- MPI
- FFTW 3.3.8+
- redev 1.1.0+ (https://github.com/SCOREC/redev)
- Omega\_h 10+ with MPI enabled (https://github.com/SCOREC/omega_h)
- Catch2 2.\* (for unit tests) (https://github.com/catchorg/Catch2/tree/v2.13.8)

## Build Instructions

SCOREC Rhel7 environment

```
module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
module load \
gcc/10.1.0 \
mpich \
cmake/3.20.0 \
fftw \
gdb
```

Build, install, and test

```
git clone git@github.com:SCOREC/xgc1_data.git #test data
git clone git@github.com:SCOREC/wdmapp_coupling.git

cmake -S wdmapp_coupling -B buildWdmCpl \
-Dredev_ROOT=/path/to/redev/install \
-DOmega_h_ROOT=/path/to/omegah/install \
-DCMAKE_INSTALL_PREFIX=$PWD/buildWdmCpl/install \
-DWDMCPL_TEST_DATA_DIR=$PWD/xgc1_data \
-DCatch2_ROOT=/path/to/catch2/install

cmake --build buildWdmCpl --target install

ctest --test-dir buildWdmCpl --output-on-failure
```

### BUILD TODO
- create a spack environment that's part of this project that can build the whole stack.
  most of the pieces are in place for this, but it will require createing a package for redev
  and of the SCOREC version of Omega\_h
  - scorec version 10.1.0 of Omega\_h is in spack@develop
    https://github.com/spack/spack/blob/8ddaa08ed2aacb4b5e587a33c625492cbdd4886e/var/spack/repos/builtin/packages/omega-h/package.py#L21
