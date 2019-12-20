#include basic makefile definitions
include defs.mk

FC := ftn

#adios - I/O library
ADIOS_INC=$(shell adios_config -f -c)
ADIOS_LIB=$(shell adios_config -f -l)

#pspline lib - spline library
PSPLINE_DIR=/ccs/home/shku/lib/ntcc_pspline/
#/ccs/proj/env003/dazevedo/titan/ntcc_pspline
PSPLINE_INC=-I$(PSPLINE_DIR)/include -I$(PSPLINE_DIR)/mod
PSPLINE_LIB=-L$(PSPLINE_DIR)/lib -lpspline -lezcdf

#performance monitor
PERFMOD_DIR=/ccs/home/shku/lib/camtimers/
#/lustre/atlas2/env003/proj-shared/camtimers
PERFMOD_INC=-I$(PERFMOD_DIR)
PERFMOD_LIB=-L$(PERFMOD_DIR) -ltimers

NETCDF=-I$(NETCDF_DIR)/include \
           -L$(NETCDF_DIR) -lnetcdf -lnetcdff

NETCDF=-I$(NETCDF_DIR)/include \
			            -L$(NETCDF_DIR) -lnetcdf -lnetcdff

VERSION_MOD = -coupling

#reset petsc fc_flags
# FC_FLAGS = -O0 -g -C -traceback -fpic
#FC_FLAGS = -g -C -qopenmp -no-ipo -traceback -fpic
#FC_FLAGS = -zero -g -check all -traceback
#FC_FLAGS = -O3 -qopenmp -no-ipo -traceback
#FC_FLAGS = -O0 -g -no-ipo -traceback
#FC_FLAGS = -O3 -no-ipo -traceback
#FC_FLAGS = -no-ipo -traceback
FC_FLAGS = -fast -D__PGI -fpic # -mp

XGC_FLAGS = -DPSPLINE \
		  -DUSE_BICUB_MOD \
		  -DUSE_ONE_D_I_CUB_MOD \
		  -UCAM_TIMERS \
		  -DADIOS \
		  -DSOLVERLU -llapack -lblas \
		  -DITER_GRID \
		  -DF_USE_MARKER_DEN \
		  -DF_USE_MARKER_DEN2 \
		  -DF0_TOR_LINEAR \
		  -DV_PERP \
		  -DDELTAF_MODE2 \
		  -DINIT_GENE_PERT \
		  -UPOS2_ORG_HACK \
		  -UFORCE_ZEROCHARGE \
		  -UFORCE_ZEROCHARGE_NEW \
		  -UOPTIM_GYRO_AVG_MAT \
		  -UNEWGYROMATRIX \
		  -DPURE_RK4 \
		  -DRK4_ACCURATE_EFIELD \
		  -DSPECIFIC_GENEXGC \
		  -DXGC_COUPLING_CORE_EDGE \
                  -DXGC_COUPLING_CORE_EDGE_FIELD \
		  -DXGC_COUPLING_CORE_EDGE_VARPI2 \
		  -UCCE_TEST_1
		  #-DPOS2_ORG
		  #-UPLANE_MAJOR \


#XGC_FLAGS = -DPSPLINE \
                  -DUSE_BICUB_MOD \
                  -DUSE_ONE_D_I_CUB_MOD \
                  -DCAM_TIMERS \
                  -DADIOS \
                  -DSOLVERLU -llapack -lblas \
                  -DDELTAF_MODE2 \
                  -USPECIFIC_GENEXGC\
                  -UXGC_COUPLING_CORE_EDGE\
                  -UXGC_COUPLING_CORE_EDGE_FIELD\
                  -UXGC_COUPLING_CORE_EDGE_VARPI2\
                  -UCCE_TEST_1
                  #-DPOS2_ORG
                  #-UPLANE_MAJOR \

XGC_INCLUDE = ${PSPLINE_INC} ${PERFMOD_INC} ${PETSC_INCLUDE_OPTS} ${ADIOS_INC} -module . -I./

LIB= ${PSPLINE_LIB} ${ADIOS_LIB} ${PETSC_LIB} ${PERFMOD_LIB}

CMP=${FC}

# Set EXTRA_OBJ to $(IMSL_OBJ) for MY_IMSL wrappers, or $(PORT_OBJ) if not using IMSL
EXTRA_OBJ=$(PORT_OBJ)

# Set MPI_OBJ to $(PAR_OBJ) for parallel code, or $(SER_OBJ) for serial code
MPI_OBJ=$(PAR_OBJ)

# include basic makefile build rules
include rules.mk


