###########################################################################
### architecture dependent GENE makefiles, draco version                ###
### http://www.mpcdf.mpg.de/services/computing/draco                    ###
###########################################################################
#									  #
# the following list of modules has to be loaded first with module load   #
MODULELIST_INTEL=intel impi mkl fftw-mpi petsc-complex slepc-complex hdf5-mpi
# 									  #
# for GNU compiler:							  #
MODULELIST_GNU=gcc impi mkl fftw-mpi petsc-complex slepc-complex hdf5-mpi
#									  #
# for PGI compiler:                                                       #
MODULELIST_PGI=intel impi pgi mkl
#                                                                         #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=intel
#COMPILER=gnu
#COMPILER=pgi

CHIP = Haswell

# FFTLIB - needed in /src/files.mk
FFTLIB = fftw
#FFTLIB = mkl

PRECISION= double
OPENMP = no

SLEPC = yes

SCALAPACK = yes

DEBUG = no

INCLUDE_SYMBOLS = yes

COMPILER_REPORTS=no

USE_PERFLIB = none

FUTILS = yes

PRODRUN = yes

#memory per core
MB_PER_CORE=2800

#PREPROC += -DWITH_DEBUG_OUTPUT

ifeq ($(COMPILER),pgi)
 SLEPC=no
 FUTILS=no
 FFTLIB=mkl
 MPI_BINDING_HOME=/draco/u/tbg/soft/impi_binding
 PGI_NAME=intel64.pgi.18.1
 MPFC=$(MPI_BINDING_HOME)/f90/intel64/bin/mpi$(PGI_NAME)
 MPCC=$(MPI_BINDING_HOME)/c/intel64/bin/mpi$(PGI_NAME)
 INCPATHS = -I$(MPI_BINDING_HOME)/f90/intel64/include/$(PGI_NAME)
 LIBS += -L$(MPI_BINDING_HOME)/f90/intel64/lib -lmpi$(PGI_NAME)
 #add $(MPI_BINDING_HOME)/f90/intel64/lib and $I_MPI_ROOT/lib64 to LD_LIBRARY_PATH
endif


###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS += -I$(MKLROOT)/include
 LIBS += $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Wl,-rpath,$(FFTW_HOME)/lib
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
    ifeq ($(COMPILER),gnu)
       PETSC_ARCH =
       PETSC_DIR = /draco/u/tbg/soft/petsc-gcc
       SLEPC_DIR = /draco/u/tbg/soft/slepc-gcc
   endif
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK = $(MKL_SCALAPACK_LINKLINE)
else
  LIB_BLAS = $(MKL_BLAS_LINKLINE)
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
# PREPROC += -DWITHPERF=1
# everything is set in the compiler specific makefiles
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

#FOR interactive use:
#MPRUN = srun -n $(N_PES) -t 15 -p express ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=
