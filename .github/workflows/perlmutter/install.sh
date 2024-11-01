#!/bin/bash -e

name=$1
branch=$2

cd $SCRATCH/globus-compute/$name-test

# # kokkos
# rm build-kokkos -rf
# # rm kokkos -rf
# git clone -b 4.2.00 https://github.com/kokkos/kokkos.git
# cmake -S kokkos -B build-kokkos \
#   -DCMAKE_INSTALL_PREFIX=build-kokkos/install \
#   -DCMAKE_BUILD_TYPE="Release" \
#   -DCMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper \
#   -DKokkos_ARCH_AMPERE80=ON \
#   -DKokkos_ENABLE_SERIAL=ON \
#   -DKokkos_ENABLE_OPENMP=off \
#   -DKokkos_ENABLE_CUDA=on \
#   -DKokkos_ENABLE_CUDA_LAMBDA=on \
#   -DKokkos_ENABLE_CUDA_CONSTEXPR=on \
#   -DKokkos_ENABLE_DEBUG=off
# cmake --build build-kokkos -j 24 --target install

# # ADIOS2
# rm build-ADIOS2 -rf
# # rm ADIOS2 -rf
# git clone git@github.com:ornladios/ADIOS2.git
# cmake -S ADIOS2/ -B build-ADIOS2 \
#   -DCMAKE_INSTALL_PREFIX=build-ADIOS2/install \
#   -DADIOS2_USE_CUDA=on \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DCMAKE_C_COMPILER=cc
# cmake --build build-ADIOS2 --target install -j8

# # perfstubs
# rm build-perfstubs -rf
# # rm perfstubs -rf
# git clone git@github.com:UO-OACISS/perfstubs.git
# cmake -S perfstubs -B build-perfstubs \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DCMAKE_C_COMPILER=cc \
#   -DCMAKE_INSTALL_PREFIX=build-perfstubs/install
# cmake --build build-perfstubs -j2 --target install

# # redev
# rm build-redev -rf
# # rm redev -rf
# git clone git@github.com:SCOREC/redev.git
# cmake -S redev -B build-redev \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DCMAKE_C_COMPILER=cc \
#   -DMPIEXEC_EXECUTABLE=`which srun` \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_INSTALL_PREFIX=build-redev/install \
#   -Dperfstubs_DIR=$PWD/build-perfstubs \
#   -DADIOS2_ROOT=build-ADIOS2/install
# cmake --build build-redev -j2 --target install

# # omega_h
# rm build-omega_h -rf
# # rm omega_h -rf
# git clone git@github.com:SCOREC/omega_h.git
# cmake -S omega_h -B build-omega_h \
#   -DCMAKE_INSTALL_PREFIX=build-omega_h/install \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DBUILD_SHARED_LIBS=off \
#   -DOmega_h_USE_Kokkos=ON \
#   -DOmega_h_USE_CUDA=on \
#   -DOmega_h_CUDA_ARCH=80 \
#   -DOmega_h_USE_MPI=on  \
#   -DMPIEXEC_EXECUTABLE=srun \
#   -DBUILD_TESTING=off  \
#   -DCMAKE_C_COMPILER=cc \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DKokkos_PREFIX=build-kokkos/install/lib64/cmake
# cmake --build build-omega_h -j24 --target install

# # Catch2
# rm build-Catch2 -rf
# # rm Catch2 -rf
# git clone https://github.com/catchorg/Catch2
# cmake -S Catch2 -B build-Catch2 \
#   -DCMAKE_INSTALL_PREFIX=$PWD/build-Catch2/install \
#   -DCMAKE_CXX_COMPILER=CC
# cmake --build build-Catch2 -j2 --target install

# pcms
rm pcms -rf
rm build-pcms -rf
rm pcms_testcases -rf

git clone https://github.com/SCOREC/pcms.git
cd pcms && git checkout $branch && cd -
git clone git@github.com:jacobmerson/pcms_testcases.git

cmake -S pcms -B build-pcms \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_BUILD_TYPE=Release \
  -Dredev_DIR=$PWD/build-redev/install/lib64/cmake/redev \
  -Dperfstubs_DIR=$PWD/build-perfstubs \
  -DOmega_h_DIR=$PWD/build-omega_h/install/lib64/cmake/Omega_h/ \
  -DKokkos_DIR=$PWD/build-kokkos/install/lib64/cmake/Kokkos/ \
  -DCatch2_DIR=$PWD/build-Catch2/install/lib64/cmake/Catch2/ \
  -DPCMS_TEST_DATA_DIR=$PWD/pcms_testcases
cmake --build build-pcms -j8