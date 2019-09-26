###########################################################################
### architecture dependent GENE makefile for EDISON (NERSC)             ###
###  http://www.nersc.gov/systems/edison-cray-xc30/                     ###
###########################################################################

###########################################################################
### List of required environment modules to be called by module load    ###
### (can be added once to your ~/.cshrc.ext or ~/.bashrc.ext file)      ###
###                                                                     ###
### For INTEL Compiler / Programming Environment:                       ###

MODULELIST_INTEL=\
modules nsg craype craype-network-aries PrgEnv-intel atp craype-ivybridge \
cray-mpich altd darshan cray-fftw cray-libsci \
cray-petsc-complex/3.7.6.0 hdf5-parallel/1.10.1

### All petsc modules are linked against hdf5/1.10.0.X, hence this      ###
### particular call order.                                              ###
###                                                                     ###
### For CRAY Compiler / Programming Environment:                        ###

MODULELIST_CRAY=\
modules nsg craype craype-network-aries PrgEnv-cray atp craype-ivybridge \
cray-mpich altd darshan cray-fftw cray-libsci \
cray-petsc-complex/3.7.6.0 cray-hdf5-parallel

### For GNU Compiler / Programming Environment:                         ###

MODULELIST_GNU=\
modules nsg craype craype-network-aries PrgEnv-gnu atp craype-ivybridge \
cray-mpich altd darshan cray-fftw cray-libsci \
cray-petsc-complex cray-hdf5-parallel

### Common Modules for post-processing etc.:                            ###
MODULELIST_COMMON=gnuplot idl git
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

FUTILS = no

MB_PER_CORE=2400
#officially: 2.67GB/core
#however, the MPI environment itself typically consumes some fraction

INCLUDE_SYMBOLS = no

COMPILER_REPORTS = no

ifneq ($(filter cray gnu,$(COMPILER)),)
SLEPC = no
endif

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
ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE =
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 ifeq ($(FFTW_INC),)
  $(warning run module load fftw first)
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
 ifeq ($(COMPILER),intel)
  PETSC_DIR ?=
  SLEPC_DIR =/global/homes/m/merlo/soft/software/SLEPc/3.7.4-CrayIntel-2017.01-complex
  ifeq ($(strip $(PETSC_ARCH)),)
   PETSC_ARCH = arch-installed-petsc
  endif
 endif
endif

#SCALAPACK is already loaded in ftn via cray-libsci
ifeq ($(SCALAPACK),yes)
  LIBS +=
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
ifeq ($(COMPILER),cray)
NOOPTLIST = par_geom.o
else
NOOPTLIST =
endif
