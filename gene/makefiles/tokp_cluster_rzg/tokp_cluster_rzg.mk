###########################################################################
### architecture dependent GENE makefiles, TOK-P cluster version        ###
###########################################################################
#									  #
# the following list of modules has to be loaded first with module load   #
MODULELIST_INTEL=intel/17.0 impi/2017 mkl/2017 fftw petsc-cplx slepc-cplx hdf5-mpi
# 									  #
# for GNU compiler:							  #
MODULELIST_GNU=gcc/6.4 impi/2017 mkl/2017 fftw/gcc petsc-cplx slepc-cplx hdf5-mpi-gcc
#									  #
# for PGI compiler:                                                       #
MODULELIST_PGI=pgi impi mkl
#                                                                         #
MODULELIST_COMMON=perflib idl
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=intel

CHIP = SandyBridge

# set to: mkl, fftw
FFTLIB = fftw

PRECISION = double
OPENMP = no

SLEPC = yes
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no
USE_PERFLIB = none
FUTILS = yes
PRODRUN = yes

#memory per core
MB_PER_CORE=3500

ifeq ($(COMPILER),pgi)
 SLEPC=no
 FUTIILS=no
endif

ifeq ($(COMPILER),gnu)
 COMPILE_MPI_MOD=yes
endif

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS =

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS += -I$(MKLROOT)/include
 LIBS +=  $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   PETSC_ARCH =
   PETSC_DIR = /tokp/work/tbg/soft/petsc/
   SLEPC_DIR = /tokp/work/tbg/soft/slepc/
   ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
    include $(SLEPC_DIR)/conf/slepc_common
   else
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
   endif
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
#  compiling with SLEPc requires libsvml to avoid unresolved dependencies
   LIBS += -L$(IFORT_BASE)/compiler/lib/intel64 -lsvml
   LIBS += $(SLEPC_LIB)
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
else
#Insert only BLAS library
  LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif


###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################

ifeq ($(USE_PERFLIB),perf)
 LIBS += -L$(PERFLIB_HOME)/lib -looperf
 PREPROC += -DWITHPERF=1
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

#LDFLAGS += -Xlinker -M

###########################################################################
### Running                                                             ###
###########################################################################
MPRUN = export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=yes;\
	mpiexec -n $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=

preproc: $(F90PRENAMES) $(PPDIR)/gene.f90
