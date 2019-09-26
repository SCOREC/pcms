###########################################################################
### architecture dependent GENE makefiles				###
### Comet@SDSC                                                          ###
### https://www.sdsc.edu/support/user_guides/comet.html                 ###
###########################################################################

###########################################################################
### List of required environment settings/modules:                      ###
### module load fftw hdf5 mkl                                           ###
MODULELIST_INTEL=intel mkl fftw hdf5
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = intel
#uncomment the following lines on CRAY machine
#MPFC = ftn
#MPCC = cc
#COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP = Xeon_E5_IvyBridge

ARCHIVE = ar r

MPRUN = mpirun_rsh -hostfile $(SLURM_NODEFILE) -np $(N_PES) ./$(EXEC)


###########################################################################
### SWITCHES                                                            ###
###########################################################################


# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

#double precision should be default; single precision, however, might
#be sufficient and much faster
PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= no

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

#performance profiling library
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = no

#Provide an upper (RAM) memory limit in MB per core for the code internal
#performance optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=4500
#Note: the system itself often requires a significant fraction
#      which should be taken into account

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS =  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include

#LIBRARIES AND LIBFLAGS
LIBS =

#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 INCPATHS += -I$(MKLROOT)/include
 LIBS +=  $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(FFTWHOME)/include
 LIBS += -L$(FFTWHOME)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if those variables are not
# set by modules (make sure, you use the *complex* versions)
# PETSC_ARCH =
 PETSC_DIR = /home/dtold/lib/petsc
 SLEPC_DIR = /home/dtold/lib/slepc
endif

#Comet: BLAS, LAPACK, etc. are all taken from MKL
ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
#  LIBS += -L$(SCALAPACKHOME)/lib -lscalapack
else
  LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
