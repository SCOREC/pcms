###########################################################################
### architecture dependent GENE makefiles for BULL supercomputer helios ###
### https://www.iferc-csc.org/                                          ###
###########################################################################
#                                                                         #
# the following modules have to be loaded first (e.g., in ~/.bash_profile)#
# module load bullxmpi intel                                              #
# module load fftw hdf5_p                                                 #
# module load mxml adios                                                  #
#                                                                         #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER = intel
CHIP = SandyBridge

#ARCHIVE command
ARCHIVE = ar r

MPRUN = export OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
        export MKL_SERIAL=yes;\
        mpiexec -n $(N_PES) $(EXEC)


# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

PRECISION= double

# Choose MPIVENDOR=intel for Intel MPI and bullx for Bull MPI
ifneq ($(I_MPI_ROOT),)
MPIVENDOR = intelmpi
else
MPIVENDOR = bullxmpi
endif

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= yes

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

#performance profiling library
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required for some geometry interfaces
FUTILS = no

# HLST ADIOS CHECKPOINT
ifneq ($(ADIOS_DIR),)
 HAC = yes
else
 HAC = no
endif

#Provide an upper (RAM) memory limit for the code internal performance
#optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=3500
#Note: the system itself often requires a significant fraction
#      which should be taken into account

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

MKLVERSION = 11.2
MKLROOT?= /csc/softs/intel/mkl/

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS
LIBS =

#FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS += $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(FFTW_DIR)/include
 LIBS += -L$(FFTW_DIR)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(findstring intel,$(MPIVENDOR)),intel)
#hdf5 is not installed yet for intelmpi
FUTILS=no
endif

#add FUTILS library before PETSc/SLEPc in order to avoid problems
#with additional hdf5 linkages in these libraries
#which cause empty scalar entries
ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR) -I$(HDF5PATH)/include

endif

#Use self-compiled libraries due to occasional problems in the testsuite
#when using the official ones
ifeq ($(SLEPC),yes)
 PETSC_ARCH =
 PETSC_DIR = /csc/workdir3/tbg/soft/petsc-new/
 SLEPC_DIR = /csc/workdir3/tbg/soft/slepc-new/
 ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
  include $(SLEPC_DIR)/conf/slepc_common
 else
  include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 endif
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
 LIBS += $(MKL_SCALAPACK_LINKLINE)
else
 LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(HAC),yes)
 #  HLST-ADIOS-CHECKPOINT interface
 ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
 ADIOS_LIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
endif


ifeq ($(USE_PERFLIB),perf)
 LIBS += -L/csc/home0/dannert/helios_inst/lib -looperf_r -lpfm -lstdc++
 LD_LIBRARY_PATH += /csc/home0/dannert/helios_inst/lib
endif


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################


