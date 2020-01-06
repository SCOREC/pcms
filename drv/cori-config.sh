#/bin/bash -x

#edit these parameters

IMPLICIT_PETSC=OFF
PREFIX=${PWD}/../../install/drv

cmake ..\
	-DCMAKE_C_COMPILER=cc	\
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_SCOREC_FIR=$PUMI_ROOT/lib/cmake/SCOREC	\
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_TESTING=ON	\
	-DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DMPI_INCLUDE_PATH="/opt/cray/pe/craype/2.5.18/bin/ftn/include" \

# ${CONFIG_PARAMS}
