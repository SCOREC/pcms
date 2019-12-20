#/bin/bash -x

#edit these parameters

IMPLICIT_PETSC=OFF
PREFIX=${PWD}/../../../coupling_NL/GENE_XGC/circular_NL/XGC

cmake ..\
	-DCMAKE_C_COMPILER=cc	\
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_SCOREC_FIR=$PUMI_ROOT/lib/cmake/SCOREC	\
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_TESTING=ON	\
	-DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DCMAKE_INCLUDE_DIRS="$PWD/../cpl/cmake/adios"\
	-DMPI_INCLUDE_PATH="/opt/cray/pe/craype/2.5.18/bin/ftn/include" \

# ${CONFIG_PARAMS}
