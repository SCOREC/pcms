###########################################################################
### architecture dependent GENE makefile for PizDaint at CSCS           ###
### /http://user.cscs.ch/hardware/monte_rosa_cray_xe6/index.html        ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
###                                                                     ###
### switch to desired PrgEnv,e.g.:                                      ###
### module unload PrgEnv-cray && module load PrgEnv-intel               ###
###                                                                     ###
### then                                                                ###
###                                                                     ###
### module load fftw                                                    ###
### module load cray-hdf5-parallel/1.8.12                               ###
###                                                                     ###
### N.B. due to a problem with hdf5 we have to downgrade to version 12  ###
###########################################################################

########################################################################
### BASIC settings                                                      ###
###########################################################################

MPFC = ftn
MPCC = cc
COMPILER = intel

CHIP = Broadwell

#set all in submission script
MPRUN = srun -n $(NTASK) --ntasks-per-node $(TASKS_N)  -c $(CPU_T) ./$(EXEC)

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

FUTILS = yes

MB_PER_CORE=env

INCLUDE_SYMBOLS = no

COMPILER_REPORTS = no

PRODRUN=yes
###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.
INCPATHS =

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS =

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
 PETSC_ARCH ?=
 PETSC_DIR ?=
 SLEPC_DIR =/users/merlo/easybuild/daint/$(shell echo $(CHIP) | tr A-Z a-z)/software/SLEPc/3.7.3-CrayIntel-2016.11-complex/
 ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
  include $(SLEPC_DIR)/conf/slepc_common
 else
  include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 endif
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
  LIBS +=
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_DIR)
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
