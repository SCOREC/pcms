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

PREFIX=${INSTALLROOT}/wdm_coupling

cmake .. \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_Fortran_COMPILER=$FTN \
      -DSCOREC_DIR=$PUMI_ROOT/lib/cmake/SCOREC/ \
      -DENABLE_COMPLEX=OFF \
      -DENABLE_TESTING=ON \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DADIOS_INCLUDE_DIRS="$DEVROOT/wdm/wdmapp_coupling/cpl/cmake/adios" \
      -DMPI_INCLUDE_PATH="/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.3-diz4f6ieln25ouifyc7ndtqlfksom6nb/include" \
#      ${CONFIG_PARAMS}

