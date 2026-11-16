#!/bin/bash -x
(
#cdash output root
d=/users/d_zxg06726/nightlyBuilds/pcms_build
exec > $d/nightly_log.txt 2>&1

source /etc/profile
# source /users/d_zxg06726/.bash_profile

#setup lmod
export PATH=/usr/share/lmod/lmod/libexec:$PATH

#setup spack modules
unset MODULEPATH

module use /opt/scorec/spack/rhel9/v0222_2/lmod/linux-rhel9-x86_64/Core/
module load gcc/13.2.0-4eahhas
module load mpich/4.2.3-62uy3hd
module load cmake/3.30.5-5e54py4
module load cuda/12.6.2-gqq65nw
module load openblas/0.3.28-eubd5ed
export LD_LIBRARY_PATH=/opt/scorec/spack/rhel9/v0222_2/install/linux-rhel9-x86_64_v3/gcc-13.2.0/mpich-4.2.3-62uy3hdwqe75fjfskrrysco5i6gqjblv/lib:$LD_LIBRARY_PATH
export PETSC_OPTIONS="-use_gpu_aware_mpi 0"

cd $d
#remove compilation directories created by previous nightly.cmake runs
[ -d build ] && rm -rf build/

#install kokkos
[ ! -d kokkos ] && git clone --branch 4.6.01 --depth 1 https://github.com/kokkos/kokkos.git
[ -d build-kokkos ] && rm -rf build-kokkos
cmake -S kokkos -B build-kokkos \
  -DCMAKE_INSTALL_PREFIX=build-kokkos/install \
  -DCMAKE_BUILD_TYPE="Release" \
  -DCMAKE_CXX_COMPILER=$d/kokkos/bin/nvcc_wrapper \
  -DKokkos_ARCH_AMPERE80=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_OPENMP=off \
  -DKokkos_ENABLE_CUDA=on \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=on \
  -DKokkos_ENABLE_DEBUG=off
cmake --build build-kokkos -j 8 --target install

#install kokkos kernels
[ ! -d kokkos-kernels ] && git clone --branch 4.6.01 --depth 1 https://github.com/kokkos/kokkos-kernels.git
[ -d build-kokkos-kernels ] && rm -rf build-kokkos-kernels
cmake -S kokkos-kernels -B build-kokkos-kernels \
  -DCMAKE_INSTALL_PREFIX=build-kokkos-kernels/install \
  -DKokkos_ROOT=$d/build-kokkos/install/lib64/cmake \
  -DBUILD_SHARED_LIBS=off
cmake --build build-kokkos-kernels -j 8 --target install

#install ADIOS2
[ ! -d ADIOS2 ] && git clone --branch v2.10.2 https://github.com/ornladios/ADIOS2.git
[ -d build-ADIOS2 ] && rm -rf build-ADIOS2
cmake -S ADIOS2 -B build-ADIOS2 \
  -DCMAKE_INSTALL_PREFIX=$d/build-ADIOS2/install \
  -DADIOS2_USE_CUDA=on \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DADIOS2_USE_ZFP=off
cmake --build build-ADIOS2 -j 8 --target install

#install perfstubs
[ ! -d perfstubs ] && git clone https://github.com/UO-OACISS/perfstubs.git
cd perfstubs && git checkout a1fa3feb1d89214e28047f166500351074b5f0c2 && cd -
[ -d build-perfstubs ] && rm -rf build-perfstubs
cmake -S perfstubs -B build-perfstubs \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_INSTALL_PREFIX=$d/build-perfstubs/install
cmake --build build-perfstubs -j 8 --target install

#install redev
[ ! -d redev ] && git clone https://github.com/SCOREC/redev.git
cd redev && git pull && cd -
[ -d build-redev ] && rm -rf build-redev
cmake -S redev -B build-redev \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_INSTALL_PREFIX=build-redev/install \
  -DADIOS2_DIR=$d/build-ADIOS2/install/lib64/cmake/adios2 \
  -Dperfstubs_DIR=$d/build-perfstubs/install/lib/cmake
cmake --build build-redev -j 8 --target install

#install omega_h
[ ! -d omega_h ] && git clone https://github.com/SCOREC/omega_h.git
cd omega_h && git pull && cd -
[ -d build-omega_h ] && rm -rf build-omega_h
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
  -DKokkos_PREFIX=$d/build-kokkos/install/lib64/cmake
cmake --build build-omega_h -j 8 --target install

#install meshFields
[ ! -d meshFields ] && git clone https://github.com/SCOREC/meshFields
cd meshFields && git pull && cd -
[ -d build-meshFields ] && rm -rf build-meshFields
cmake -S meshFields -B build-meshFields \
  -DCMAKE_INSTALL_PREFIX=build-meshFields/install \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc \
  -DMPIEXEC_EXECUTABLE=mpirun \
  -DKokkos_DIR=$d/build-kokkos/install/lib64/cmake/Kokkos \
  -DOmega_h_DIR=$d/build-omega_h/install/lib64/cmake/Omega_h
cmake --build build-meshFields -j 8 --target install

#install Catch2
[ ! -d Catch2 ] && git clone --branch v3.11.0 https://github.com/catchorg/Catch2
[ -d build-Catch2 ] && rm -rf build-Catch2
cmake -S Catch2 -B build-Catch2 \
  -DCMAKE_INSTALL_PREFIX=$d/build-Catch2/install
cmake --build build-Catch2 -j 8 --target install

[ ! -d petsc ] && git clone --branch v3.24.2 https://gitlab.com/petsc/petsc.git
cd petsc
[ -d cuda-kokkos ] && rm -rf cuda-kokkos
./configure \
  PETSC_ARCH=cuda-kokkos \
  --with-kokkos-dir=$d/build-kokkos/install/ \
  --with-kokkos-kernels-dir=$d/build-kokkos-kernels/install/ \
  --with-cuda=1 \
  --with-shared-libraries=0 \
  --with-openblas-dir="${OPENBLAS_RHEL9_ROOT}"
make all
cd -

#download testcases
[ ! -d pcms_testcases ] && git clone https://github.com/jacobmerson/pcms_testcases.git
cd pcms_testcases && git pull && cd -

touch $d/startedCoreNightly
#run nightly.cmake script
ctest -V --script $d/nightly.cmake
touch $d/doneCoreNightly
)
