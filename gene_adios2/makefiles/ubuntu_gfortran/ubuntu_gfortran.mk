###########################################################################
### architecture dependent GENE makefiles, Ubuntu+gfortran laptop version #
###########################################################################
#									  #
# the following packages have to be installed first:		          #
# apt-get install gfortran openmpi-bin libopenmpi-dev 	        	  #
# apt-get install libfftw3-dev liblapack-dev libatlas-base-dev		  #
# apt-get install scalapack-mpi-test libhdf5-openmpi-dev		  #
#								          #
# for testsuite script							  #
# apt-get install perl libswitch-perl                                     #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=gnu
CHIP=

FFTLIB = fftw
PRECISION= double
OPENMP = no

MPI_IO = no
SLEPC = no
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = yes
PRODRUN = no

#memory per core
MB_PER_CORE=1000

COMPILE_MPI_MOD=yes

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)
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
 ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
  include $(SLEPC_DIR)/conf/slepc_common
 else
  include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 endif
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS += $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
# uncomment for static libraries (not available on all ubuntu version)
# LIBS += -L/usr/lib -lscalapack-openmpi -lblacs-openmpi -lblacsF77init-openmpi
# shared libraries:
  LIBS += /usr/lib/libscalapack-openmpi.so.1 /usr/lib/libblacs-openmpi.so.1 /usr/lib/libblacsF77init-openmpi.so.1
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

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(USE_PERFLIB),perf)
 LIBS += $(PERFLIB_TS)
 PREPROC += -DWITHPERF=1
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

###########################################################################
### Running                                                             ###
###########################################################################

MPRUN = export OMP_NUMTHREADS=1; export MKL_SERIAL=yes;\
      mpiexec -np $(N_PES) ./$(EXEC)

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=

