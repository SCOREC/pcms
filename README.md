# PCMS: Parallel Coupler For Multimodel Simulations

Adios2-based xgc_coupler for XGC and GENE

## Dependencies

- CMake 3.19+
- MPI
- FFTW 3.3.8+
- redev 3.0.0+ (https://github.com/SCOREC/redev)
- Omega\_h 10.2.0+ with MPI enabled (https://github.com/SCOREC/omega_h)
- Catch2 2.\* (for unit tests) (https://github.com/catchorg/Catch2/tree/v2.13.8)

## Build Instructions

## Build with modules
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
git clone git@github.com:SCOREC/wdmapp_testcases.git #test data
git clone git@github.com:SCOREC/wdmapp_coupling.git

cmake -S wdmapp_coupling -B buildWdmCpl \
-Dredev_ROOT=/path/to/redev/install \
-DOmega_h_ROOT=/path/to/omegah/install \
-DCMAKE_INSTALL_PREFIX=$PWD/buildWdmCpl/install \
-DPCMS_TEST_DATA_DIR=$PWD/wdmapp_testcases \
-DCatch2_ROOT=/path/to/catch2/install

cmake --build buildWdmCpl --target install

ctest --test-dir buildWdmCpl --output-on-failure
```

## Spack based build
1. Install spack
   ```console
   $ mkdir /lore/$USER/spack
   $ cd /lore/$USER/spack
   $ git clone -c feature.manyFiles=true -b releases/v0.20 https://github.com/spack/spack.git
   $ . spack/share/spack/setup-env.sh
   ```
   We can also add the spack setup line into the `~/.bashrc` with `echo ". spack/share/spack/setup-env.sh" >> ~/.bashrc". This will load the spack setup script every time we start our terminal session.

2. Get PCMS spack repo
   The following commands will add the pcms recipe files to spack. They are not currently installed inthe upstream spack repository.
   ```console
   $ git clone https://github.com/jacobmerson/pcms-spack.git
   $ spack repo add pcms-spack/pcms
   ```
   
3. Add Spack binary mirror
   Addding the binary mirrors will avoid some compilation by downloading prebuilt binaries when available.
   ```console
   $ spack mirror add v0.20.1 https://binaries.spack.io/v0.20.1
   $ spack buildcache keys --install --trust
   ```
5. Install PCMS repo
    ```console
    $ mkdir /lore/$USER/pcms-coupler
    $ cd /lore/$USER/pcms-coupler
    $ git clone -b pcms-spack https://github.com/jacobmerson/pcms
    $ cd pcms/spack
    $ spack env create -d env spack.yaml
    $ cd env
    $ spack env activate .
    $ spack install
    ```
    
At this point hopefully, spack will now install all of the relavant dependencies and a baseline build of PCMS. The default environment has PCMS in develop mode. To modify and recompile PCMS you can modify the code and rerun `spack install`.


### BUILD TODO
- create a spack environment that's part of this project that can build the whole stack.
  most of the pieces are in place for this, but it will require createing a package for redev
  and of the SCOREC version of Omega\_h
  - scorec version 10.1.0 of Omega\_h is in spack@develop
    https://github.com/spack/spack/blob/8ddaa08ed2aacb4b5e587a33c625492cbdd4886e/var/spack/repos/builtin/packages/omega-h/package.py#L21

Details instructions for a few systems are available on the wiki.

Another Executing Approach: One would also comment out the BUILDING_TESTING in CMakeFiles.txt included in test folder; 

Assign the full path of testdatas to test_dir in the test_init.cc file; Use "mpirun -np 4 bin/test_init" for the execution.   

## Code Notes

- `Part1` refers to the core (GENE/GEM) application
- `Part3` refers to the edges (XGC) application

