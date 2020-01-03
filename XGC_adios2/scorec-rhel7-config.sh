#/bin/bash -x

# edit these params
DEVLOC=users
DEVDIR=dev
CC=mpicc
CXX=mpicxx
FTN=mpif90
IMPLICIT_PETSC=OFF

# derived vars
USER=`whoami`
DEVROOT=/${DEVLOC}/${USER}/${DEVDIR}
INSTALLROOT=${DEVROOT}/wdm/wdmapp_coupling/install

cmake .. \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_Fortran_COMPILER=$FTN \
      -DPETSC_EXECUTABLE_RUNS=ON \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$INSTALLROOT \
#      -DMPI_INCLUDE_PATH="/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.3-diz4f6ieln25ouifyc7ndtqlfksom6nb/include" \
#      ${CONFIG_PARAMS}

