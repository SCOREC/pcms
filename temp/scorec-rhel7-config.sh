#/bin/bash -x

# edit these params
DEVLOC=users
DEVDIR=dev
CC=mpicc
CXX=mpicxx
FTN=mpif90
IMPLICIT_PETSC=OFF
PREFIX=${PWD}/../

cmake ..\
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DADIOS_INCLUDE_DIRS="$PWD/../../cpl/cmake" \
#      ${CONFIG_PARAMS}

