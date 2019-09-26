###########################################################################
### architecture dependent GENE makefile for CORI (NERSC)               ###
###  http://www.nersc.gov/users/computational-systems/cori              ###
###########################################################################

###########################################################################
### List of required environment modules to be called by module load    ###
### (can be added once to your ~/.cshrc.ext or ~/.bashrc.ext file)      ###
###                                                                     ###
### Switch to KNL environment first if your code should be run on KNL(!)###
#module swap craype-haswell craype-mic-knl	                        ###
###                                                                     ###
### For INTEL Compiler / Programming Environment                        ###

MODULELIST_INTEL=\
modules nsg craype craype-network-aries PrgEnv-intel \
atp craype-haswell cray-mpich cray-fftw cray-petsc-complex \
cray-hdf5-parallel altd
###                                                                     ###
### For CRAY Compiler / Programming Environment:                        ###

MODULELIST_CRAY=\
modules nsg craype craype-hugepages2M craype-haswell craype-network-aries \
PrgEnv-cray atp  cray-mpich cray-fftw cray-petsc-complex \
cray-hdf5-parallel altd

###                                                                     ###
### For GNU Compiler / Programming Environment:                         ###

MODULELIST_GNU=\
modules nsg craype craype-network-aries PrgEnv-gnu atp craype-haswell \
cray-mpich cray-fftw cray-petsc-complex cray-hdf5-parallel altd

###                                                                     ###
###                                                                     ###
MODULELIST_COMMON=git
###                                                                     ###
### Note: Above module lists are partially taken from the NERSC startup ###
### scripts. Please check for updates of the following files:           ###
### /opt/cray/pe/modules/default/etc/modules.sh                         ###
### /etc/bash.bashrc.local                                              ###
### $HOME/.bash_profile                                                 ###
###                                                                     ###
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
#FFTLIB = mkl

PRECISION = double

DEBUG= no

ifeq ($(filter cray gnu,$(COMPILER)),)
SLEPC = yes
else
SLEPC = no
endif

SCALAPACK = yes

OPENMP = no

USE_PERFLIB = none

FUTILS = no

PRODRUN = yes

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
ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 LIBS += $(MKL_LINKLINE)
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

ifeq ($(SLEPC),yes)
# set by module load petsc-complex:
 PETSC_DIR ?=
 SLEPC_DIR ?= /global/homes/t/tbg/soft/cori/slepc-3.9.2-installed-$(CRAY_CPU_TARGET)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK +=
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH =
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
