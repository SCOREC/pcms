###########################################################################
### template for architecture dependent GENE makefiles                  ###
### OCCIGEN		    		     				###
### http://www.cines.fr/calcul/materiels/occigen	 	        ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
### module load intel bullxmpi                                          ###
### module load hdf5 fftw3                                              ###
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

COMPILER = intel
CHIP = Haswell

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
#MPRUN = mpiexec -n $(N_PES) ./$(EXEC)
#MPRUN = mpirun -np $(N_PES) ./$(EXEC)
#MPRUN = aprun -n $(N_PES) ./$(EXEC)
MPRUN = srun --mpi=pmi2 -K1 --resv-ports -n $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

#double precision should be default; single precision, however, might
#be sufficient and much faster
PRECISION= double

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

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = no

#Provide an upper (RAM) memory limit in MB per core for the code internal
#performance optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=2000
#Note: the system itself often requires a significant fraction
#      which should be taken into account

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

MKLVERSION ?= 15.0

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS +=

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS +=

#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE =
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(FFTW3_INC_DIR)
 LIBS += -L$(FFTW3_LIB_DIR) -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH =
 PETSC_DIR = /panfs/panasas/cnt0022/lpp6884/SHARED/soft/petsc
 SLEPC_DIR = /panfs/panasas/cnt0022/lpp6884/SHARED/soft/slepc

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
#  LIBS += $(MKL_SCALAPACK_LINKLINE)
  LIBS += $(MKL_SCA_LIBS)
else
#  LIBS += $(MKL_BLAS_LINKLINE)
  LIBS += $(MKL_LIBS)
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
