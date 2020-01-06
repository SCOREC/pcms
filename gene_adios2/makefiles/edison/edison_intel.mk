###########################################################################
### architecture dependent GENE makefile for EDISON (NERSC)             ###
###  http://www.nersc.gov/systems/edison-cray-xc30/                     ###
###########################################################################

###########################################################################
### List of required environment settings/modules:                      ###
### (copy lines, e.g., in your ~/.cshrc.ext)                            ###
###                                                                     ###
### module load fftw cray-libsci                                        ###
### module load gv gnuplot idl                                          ###
### module load cray-hdf5-parallel                                      ###
###                                                                     ###
### For PrgEnv-intel:                                                   ###
### module load cray-petsc-complex/3.7.6.0                              ###
### module load hdf5-parallel/1.10.1                                    ###
###									###
### All petsc modules are linked against hdf5/1.10.0.X, hence the       ###
### different linkage order.                                            ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use cray wrappers for compiler names:
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP = Xeon_E5_IvyBridge

ARCHIVE = ar r

NTASKS ?= 24
MPRUN = srun -n $(SLURM_NTASKS) ./$(EXEC)

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

MB_PER_CORE=2400
#officially: 2.67GB/core
#however, the MPI environment itself typically consumes some fraction

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
ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_PARALLEL_DIR)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(SLEPC),yes)
# set by module load petsc-complex:
 PETSC_DIR ?=
 SLEPC_DIR =/global/homes/m/merlo/soft/software/SLEPc/3.7.4-CrayIntel-2017.01-complex
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

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
