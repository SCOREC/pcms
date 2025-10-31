# PCMS: Parallel Coupler For Multimodel Simulations

Adios2-based xgc_coupler for XGC and GENE

## Dependencies

- CMake 3.19+
- MPI
- Kokkos
- KokkosKernels
- redev 3.0.0+ (https://github.com/SCOREC/redev)
- Omega\_h 10.2.0+ with MPI enabled (https://github.com/SCOREC/omega_h)
- Catch2 2.\* (for unit tests) (https://github.com/catchorg/Catch2/tree/v2.13.8)

## Build Instructions

## Build with modules
SCOREC Rhel7 environment

```
module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno
module load mpich/4.1.1-xpoyz4t
module load cmake/3.26.3-2duxfcd
// no need to load CUDA if not using GPUs
module load cuda/12.1.1-zxa4msk
```

### Build dependencies with CPU

```
git clone --branch 4.6.01 --depth 1 git@github.com:kokkos/kokkos.git
cmake -S kokkos -B build-kokkos \
  -DCMAKE_INSTALL_PREFIX=build-kokkos/install \
  -DCMAKE_CXX_STANDARD=17 \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_OPENMP=OFF \
  -DKokkos_ENABLE_CUDA=OFF \
  -DKokkos_ENABLE_CUDA_LAMBDA=OFF \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=OFF \
  -DBUILD_SHARED_LIBS=OFF
cmake --build build-kokkos --target install

git clone --branch 4.6.01 --depth 1 git@github.com:kokkos/kokkos-kernels.git
cmake -S kokkos-kernels -B build-kokkos-kernels \
  -DCMAKE_INSTALL_PREFIX=build-kokkos-kernels/install \
  -DCMAKE_CXX_STANDARD=17 \
  -DKokkos_ROOT=$PWD/build-kokkos/install/lib64/cmake \
  -DBUILD_SHARED_LIBS=OFF
cmake --build build-kokkos-kernels --target install

git clone git@github.com:SCOREC/omega_h.git
cmake -S omega_h -B build-omega_h \
  -DCMAKE_INSTALL_PREFIX=$PWD/build-omega_h/install \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DBUILD_SHARED_LIBS=OFF \
  -DOmega_h_USE_Kokkos=ON \
  -DOmega_h_USE_CUDA=OFF \
  -DOmega_h_USE_MPI=ON  \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DBUILD_TESTING=OFF  \
  -DKokkos_PREFIX=$PWD/build-kokkos/install/lib64/cmake
cmake --build build-omega_h --target install

git clone git@github.com:ornladios/ADIOS2.git
cmake -S ADIOS2 -B build-ADIOS2 \
  -DCMAKE_INSTALL_PREFIX=build-ADIOS2/install \
  -DADIOS2_USE_CUDA=OFF \
  -DADIOS2_USE_ZFP=off
cmake --build build-ADIOS2 --target install

git clone git@github.com:UO-OACISS/perfstubs.git
cmake -S perfstubs -B build-perfstubs \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_INSTALL_PREFIX=build-perfstubs/install
cmake --build build-perfstubs --target install

git clone git@github.com:catchorg/Catch2
cmake -S Catch2 -B build-Catch2 \
  -DCMAKE_INSTALL_PREFIX=$PWD/build-Catch2/install
cmake --build build-Catch2 --target install

git clone git@github.com:SCOREC/redev.git
cmake -S redev -B build-redev \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_INSTALL_PREFIX=build-redev/install \
  -DADIOS2_DIR=$PWD/build-ADIOS2/install/lib64/cmake/adios2 \
  -Dperfstubs_DIR=$PWD/build-perfstubs/install/lib/cmake
cmake --build build-redev --target install
```
### Build dependencies with CUDA
```
git clone --branch 4.6.01 --depth 1 git@github.com:kokkos/kokkos.git
cmake -S kokkos -B build-kokkos \
  -DCMAKE_INSTALL_PREFIX=build-kokkos/install \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_BUILD_TYPE="Release" \
  -DCMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper \
  -DKokkos_ARCH_AMPERE80=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_OPENMP=off \
  -DKokkos_ENABLE_CUDA=on \
  -DKokkos_ENABLE_CUDA_LAMBDA=on \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=on \
  -DKokkos_ENABLE_DEBUG=off
cmake --build build-kokkos --target install

git clone --branch 4.6.01 --depth 1 git@github.com:kokkos/kokkos-kernels.git
cmake -S kokkos-kernels -B build-kokkos-kernels \
  -DCMAKE_INSTALL_PREFIX=build-kokkos-kernels/install \
  -DCMAKE_CXX_STANDARD=17 \
  -DKokkos_ROOT=$PWD/build-kokkos/install/lib64/cmake \
  -DBUILD_SHARED_LIBS=off
cmake --build build-kokkos-kernels --target install

git clone git@github.com:SCOREC/omega_h.git
cmake -S omega_h -B build-omega_h \
  -DCMAKE_INSTALL_PREFIX=build-omega_h/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=off \
  -DOmega_h_USE_Kokkos=ON \
  -DOmega_h_USE_CUDA=on \
  -DOmega_h_CUDA_ARCH=80 \
  -DOmega_h_USE_MPI=on  \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DBUILD_TESTING=off  \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DKokkos_PREFIX=$PWD/build-kokkos/install/lib64/cmake
cmake --build build-omega_h --target install

git clone git@github.com:ornladios/ADIOS2.git
cmake -S ADIOS2 -B build-ADIOS2 \
  -DCMAKE_INSTALL_PREFIX=build-ADIOS2/install \
  -DADIOS2_USE_CUDA=on \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc \
  -DADIOS2_USE_ZFP=off
cmake --build build-ADIOS2 --target install

git clone git@github.com:UO-OACISS/perfstubs.git
cmake -S perfstubs -B build-perfstubs \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_INSTALL_PREFIX=build-perfstubs/install
cmake --build build-perfstubs --target install

git clone git@github.com:SCOREC/redev.git
cmake -S redev -B build-redev \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_INSTALL_PREFIX=build-redev/install \
  -DADIOS2_DIR=$PWD/build-ADIOS2/install/lib64/cmake/adios2 \
  -Dperfstubs_DIR=$PWD/build-perfstubs/install/lib/cmake
cmake --build build-redev --target install

git clone git@github.com:catchorg/Catch2
cmake -S Catch2 -B build-Catch2 \
  -DCMAKE_INSTALL_PREFIX=$PWD/build-Catch2/install
cmake --build build-Catch2 --target install
```
### Build, install, and test pcms
```
git clone git@github.com:jacobmerson/pcms_testcases.git // test data
git clone git@github.com:SCOREC/pcms.git

cmake -S pcms -B build-pcms \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_BUILD_TYPE="Debug" \
  -DPCMS_TIMEOUT=100 \
  -Dredev_DIR=$PWD/build-redev/install/lib64/cmake/redev \
  -DOmega_h_DIR=$PWD/build-omega_h/install/lib64/cmake/Omega_h/ \
  -Dperfstubs_DIR=$PWD/build-perfstubs/install/lib/cmake \
  -DCatch2_DIR=$PWD/build-Catch2/install/lib64/cmake/Catch2/ \
  -DKokkosKernels_DIR=$PWD/build-kokkos-kernels/install/lib64/cmake/KokkosKernels/ \
  -DPCMS_TEST_DATA_DIR=$PWD/pcms_testcases
cmake --build build-pcms -j 8

ctest --test-dir build-pcms --output-on-failure
```

## Spack based build
1. Install spack
   ```console
   $ mkdir /lore/$USER/spack
   $ cd /lore/$USER/spack
   $ git clone -c feature.manyFiles=true -b releases/v1.0 https://github.com/spack/spack.git
   $ . spack/share/spack/setup-env.sh
   ```
   We can also add the spack setup line into the `~/.bashrc` with `echo ". spack/share/spack/setup-env.sh" >> ~/.bashrc". This will load the spack setup script every time we start our terminal session.

2. Get PCMS spack repo
   The following commands will add the pcms recipe files to spack. They are not currently installed inthe upstream spack repository.
   ```console
   $ git clone https://github.com/jacobmerson/pcms-spack.git
   $ spack repo add pcms-spack/spack_repo/pcms
   ```
   
3. Install PCMS repo
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

## Creating Archive for Release
`git archive --format=tar.gz -o /tmp/pcms.tar.gz --prefix=pcms/ develop`
