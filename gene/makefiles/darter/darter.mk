###########################################################################
### architecture dependent GENE makefile, Darter XC30 version           ###
### http://www.nics.tennessee.edu/computing-resources/kraken            ###
###########################################################################
#       module load fftw subversion
#       module unload PrgEnv-pgi
#       module load PrgEnv-intel
MODULELIST_INTEL = PrgEnv-intel fftw git xt-libsci
###########################################################################
### BASIC settings                                                      ###
###########################################################################

MPFC = ftn
MPCC = cc

COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)
#CHIP = istanbul

ARCHIVE = ar r

MPRUN = aprun -n $(N_PES) ./$(EXEC)


###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = fftw

PRECISION= double

OPENMP = no

SLEPC= no #yes

SCALAPACK = no #yes

DIAG_MPI_IO = no

USE_PERFLIB = none

FUTILS=no

DEBUG = no

PRODRUN = no

INCLUDE_SYMBOLS = no

COMPILER_REPORTS = no

MB_PER_CORE=900

###########################################################################
#   COMPULSORY LIBRARIES                                                  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS =


#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS =

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif

ifeq ($(FFTLIB),fftw)
 INCPATHS += $(FFTW_INCLUDE_OPTS)
 LIBS += $(FFTW_POST_LINK_OPTS)
#$(FFTW_LIB) -lfftw3f -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)                         #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH =
# PETSC_ARCH = $(COMPILER)-double-complex
 PETSC_DIR = /nics/c/home/pueschel/PETSc-intel
 SLEPC_DIR = /nics/c/home/pueschel/SLEPc-intel
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS += -lsci_$(COMPILER)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(CRAY_HDF5_DIR)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################
#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
