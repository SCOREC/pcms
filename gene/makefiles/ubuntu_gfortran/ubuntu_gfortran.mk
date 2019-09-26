###########################################################################
### architecture dependent GENE makefiles                               ###
### Ubuntu 18.04 + gfortran laptop version                              ###
###########################################################################
#									  #
# the following packages have to be installed first:		          #
# sudo apt install gfortran openmpi-bin libopenmpi-dev 	        	  #
# sudo apt install libfftw3-dev liblapack-dev libatlas-base-dev		  #
# sudo apt install scalapack-mpi-test libhdf5-openmpi-dev		  #
#								          #
# for testsuite script							  #
# sudo apt install perl 	                                          #
###########################################################################
### BASIC settings                                                      ###
###########################################################################

COMPILER = gnu
#uncomment following lines for MPICH compiler wrappers
#MPFC=mpif90.mpich
#MPCC=mpicc.mpich
#MPCXX=mpicxx.mpich

CHIP =

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = mpiexec -n $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = fftw

PRECISION= double

OPENMP = no

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed
SLEPC= no

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

USE_PERFLIB = none

FUTILS = yes

PRODRUN = no

#memory per core
MB_PER_CORE=1000

#compile src/mpimod.F90 to include mpi library?
COMPILE_MPI_MOD=yes

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS += -I/usr/local/include -I/usr/include

#LIBRARIES AND LIBFLAGS
LIBS = -L/usr/local/gfortran/lib -L/usr/lib/x86_64-linux-gnu/
LIBS += -L/usr/local/lib -L/usr/lib/ -llapack -lblas

### FFT LIBRARY
ifeq ($(FFTLIB),fftw)
 ifeq ($(PRECISION),double)
   LIBS += -lfftw3
 else
   LIBS += -lfftw3f
 endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if PETSc/SLEPc is installed and those
# variables are not set by modules (make sure to use the *complex* versions!)
# PETSC_ARCH =
# PETSC_DIR =
# SLEPC_DIR =

#by default, PETSC/SLEPC link lines from makefiles/rules.mk are used;
#add lines here if machine-specific modifications are needed
endif

ifeq ($(SCALAPACK),yes)
# uncomment for static libraries (not available on all ubuntu version)
# LIBS += -L/usr/lib -lscalapack-openmpi -lblacs-openmpi -lblacsF77init-openmpi
# shared libraries:
  LIBS += /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.0
#/usr/lib/x86_64-linux-gnu/libblacs-openmpi.so.1 /usr/lib/x86_64-linux-gnu/libblacsF77init-openmpi.so.1
else
#Insert only BLAS library
  LIBS +=
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = /usr/lib/x86_64-linux-gnu/hdf5/openmpi/
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=
