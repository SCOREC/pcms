###########################################################################
### architecture dependent GENE makefiles, draco version                ###
### http://www.mpcdf.mpg.de/services/computing/draco                    ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load intel mkl fftw petsc-cplx slepc-cplx			  #
# module load hdf5-mpi	     	                                          #
# 									  #
# for GNU compiler:							  #
# module unload fftw hdf5-mpi						  #
# module load gcc fftw/gcc hdf5-mpi/gcc/1.8.18				  #
#									  #
# Warning: impi/2017.3 (default MPI enviroment (Sep. 2017)) MPI-IO appears#
#          to be buggy; try module load switch impi/2018 if checkpointing #
#          fails                                                          #
#									  #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=intel
#COMPILER=gnu

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

COMPILER_REPORTS = no

USE_PERFLIB = none

FUTILS = yes

PRODRUN = yes

MKLVERSION = 2017
ARCH = intel64

ifeq ($(COMPILER),gnu)
#COMPILE_MPI_MOD = yes
endif

#memory per core
MB_PER_CORE=2800

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS =

#LIBRARIES AND LIBFLAGS
LIBS =

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
    ifeq ($(COMPILER),gnu)
       PETSC_ARCH =
       PETSC_DIR = /draco/u/tbg/soft/petsc-gcc
       SLEPC_DIR = /draco/u/tbg/soft/slepc-gcc
   endif
   ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
    include $(SLEPC_DIR)/conf/slepc_common
   else
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
   endif

   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
else
  LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR) -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
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
	srun -n $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=
