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

SLEPC= no

SCALAPACK = yes

OPENMP = no

USE_PERFLIB = none

DIAG_MPI_IO= no

FUTILS = yes

PRODRUN = yes

ADIOS = yes

COUPLE_XGC = yes

WITH_GPTL = yes

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

ifeq ($(ADIOS),yes)
  ADIOS_LIB = ADIOS1_LIB=-L/global/common/sw/cray/cnl7/haswell/adios/1.13.1/intel/19.0.3.199/wuurn4w/lib -ladiosf -L/global/common/sw/cray/cnl7/haswell/zlib/1.2.11/intel/19.0.3.199/bsk25he/lib -L/global/common/sw/cray/cnl7/haswell/bzip2/1.0.6/intel/19.0.3.199/p7h55tt/lib -L/global/common/sw/cray/cnl7/haswell/lz4/1.9.0/intel/19.0.3.199/qf6mmqo/lib -L/global/common/sw/cray/cnl7/haswell/c-blosc/1.16.3/intel/19.0.3.199/tn4uoyl/lib -L/global/common/sw/cray/cnl7/haswell/zfp/0.5.0/intel/19.0.3.199/cle2mg7/lib -L/global/common/sw/cray/cnl7/haswell/sz/1.4.12.3/intel/19.0.3.199/pubgujg/lib -L/global/common/sw/cray/cnl7/haswell/snappy/1.1.7/intel/19.0.3.199/vlbmtjf/lib -L/global/common/sw/cray/cnl7/haswell/zstd/1.4.0/intel/19.0.3.199/k2svlvb/lib -lblosc -Wl,-Bdynamic -lsnappy -lzstd -Wl,-Bstatic -lz -lbz2 -llz4 -Wl,-Bdynamic -lzfp -Wl,-Bstatic -lSZ -lzlib
  ADIOS_INC = $(shell adios_config -c -f)
  ADIOS_DIR = $(shell adios_config -d)
endif

ifeq ($(WITH_GPTL),yes)
  PERFMOD_DIR= /global/homes/m/merlo/soft/GPTL/gptl-v5.5.3-2-gbb58395/intel
  INCPATHS += -I$(PERFMOD_DIR)/include
  LIBS += -L$(PERFMOD_DIR)/lib -lgptl
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
