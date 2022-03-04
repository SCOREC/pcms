# WDMApp Coupling

Adios2-based coupler for XGC and GENE

## Dependencies

- CMake 3.19+
- MPI
- FFTW 3.3.8+
- redev 1.1.0+ (https://github.com/SCOREC/redev)
- Omega\_h 10+ with MPI enabled (https://github.com/SCOREC/omega_h)

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
-DWDMCPL_TEST_DATA_DIR=$PWD/xgc1_data

cmake --build buildWdmCpl --target install

ctest --test-dir buildWdmCpl --output-on-failure
```
