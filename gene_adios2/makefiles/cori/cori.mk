###########################################################################
### architecture dependent GENE makefile for CORI (NERSC)               ###
###  http://www.nersc.gov/users/computational-systems/cori              ###
###########################################################################

###########################################################################
### List of required environment settings/modules:                      ###
### (copy lines, e.g., in your ~/.cshrc.ext)                            ###
###                                                                     ###
### The Cray compiler produces incorrect results presently, but Intel   ###
### 16 has been fixed. It is thus the recommended option.               ###
###                                                                     ###
### Switch to KNL environment if your code should be run on KNL         ###
### module swap craype-haswell craype-mic-knl	                        ###
###                                                                     ###
### module unload PrgEnv-cray                                           ###
### module load PrgEnv-intel                                            ###
### module load fftw cray-petsc-complex                                 ###
### module load cray-hdf5-parallel                                      ###
###                                                                     ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use cray wrappers for compiler names:
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP = $(CRAY_CPU_TARGET)

ARCHIVE = ar r

NTASKS ?= 32
MPRUN = srun -n $(N_PES)  ./$(EXEC)

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

PRODRUN = yes

ADIOS = yes

COUPLE=no

ifeq ($(CHIP),mic-knl)
MB_PER_CORE=1300
else
MB_PER_CORE=3200
#officially: 4 GB/core
#however, the MPI environment itself typically consumes some fraction
endif

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

#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE =
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 ifeq ($(FFTW_INC),)
  $(error run module load fftw first)
 endif
 INCPATHS += -I$(FFTW_INC)
 LIBS += -L$(FFTW_DIR) -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# set by module load petsc-complex:
 PETSC_DIR ?=
 SLEPC_DIR ?= /global/homes/d/dtold/soft/slepc-3.7.1_cori
 ifeq ($(strip $(PETSC_ARCH)),)
  PETSC_ARCH = arch-installed-petsc
 endif
 ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
  include $(SLEPC_DIR)/conf/slepc_common
 else
  include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 endif
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS +=
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

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
