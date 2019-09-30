###########################################################################
### architecture dependent GENE makefile, HOPPER (NERSC) version        ###
### see http://www.nersc.gov/users/computational-systems/hopper/        ###
###########################################################################
#									  #
# NOTE: slepc/3.3.0_complex has been removed, and the newer ones          #
#       provided by Cray do not seem to work with the PGI compiler --     #
#       until this is resolved, use the Intel environment instead.	  #
#									  #
# for PrgEnv-intel:							  #
# module unload PrgEnv-pgi                                                #
# module load PrgEnv-intel cray-petsc-complex slepc-complex fftw          #
# module load cray-hdf5-parallel    				          #
#  									  #
# what worked before for PGI:						  #
# module load cray-hdf5-parallel fftw                                     #
# module load petsc/3.3_O_complex slepc/3.3_O_complex                     #
# module load cray-libsci                                                 #
#								          #
#                                                                         #
###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use cray wrappers for compiler names:
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP =

ARCHIVE = ar r

MPRUN = aprun -n $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = fftw

PRECISION = double

DEBUG= no

SLEPC= yes

SCALAPACK = yes

OPENMP = no

USE_PERFLIB = none

DIAG_MPI_IO= no

FUTILS = no

#memory per core (7.38 GB per node)
MB_PER_CORE=1100

INCLUDE_SYMBOLS = no

COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = 

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = 
#FTN links automatically to BLAS

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),fftw)
 ifeq ($(FFTW_INC),)
  $(error run module load fftw first)
 endif
 INCPATHS += -I$(FFTW_INC)
 LIBS += -L$(FFTW_DIR) -lfftw3 -lfftw3f
endif

###########################################################################
#  ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 ifeq ($(COMPILER),intel)
   PETSC_DIR ?= 
   SLEPC_DIR ?= /global/u1/t/tbg/soft/edison/slepc-3.4.3
   ifeq ($(strip $(PETSC_ARCH)),)
     PETSC_ARCH = arch-installed-petsc
   endif
 endif
 ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
  include $(SLEPC_DIR)/conf/slepc_common
 else
  include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 endif
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS += $(PETSC_LIB) $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
  #No flag necessary (automatically set in ftn by cray-libsci)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = 
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif

ifneq ($(USE_PERFLIB),none)
  PREPROC += -DWITHPERF
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST = 
