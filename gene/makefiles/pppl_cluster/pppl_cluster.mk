###########################################################################
### PPPL cluster              		     				###
### researchcomputing.pppl.gov                                          ###
###########################################################################

###########################################################################
#                                                                         #
# Last Tested 06/06/2018 S. Lazerson (lazerson@pppl.gov)                  #
# load the following modules before compilation:	                  #
# GNU Compiler:                                                           #
MODULELIST_GNU=gcc/6.1.0 acml/6.1.0/gfortran64 openmpi/3.0.0 blacs fftw \
hdf5-parallel/1.10.1 scalapack petsc_complex/3.8.3 slepc_complex/3.8.3
#                                                                         #
# Common Modules for post-processing etc.:                                #
MODULELIST_COMMON=python/3.2.2 gnuplot idl git                            #
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################
# Fix for missing vars in modules
FFTWHOME=$(FFTW_HOME)

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = gnu

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP = 

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = srun -t 4:00:00 -p dawson -n $(N_PES) ./$(EXEC)

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
MB_PER_CORE=1750
#Note: the system itself often requires a significant fraction 
#      which should be taken into account

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I$(SRCDIR) 

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LDFLAGS =
LIBS = -L$(ACML_HOME)/lib -lacml

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
 INCPATHS += -I$(FFTWHOME)/include
 LIBS += -L$(FFTWHOME)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
#use default settings from makefiles/rules.mk
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  INCPATHS += -I${BLACS_HOME}/include
  LIBS += -L${SCALAPACK_HOME} -lscalapack
  LIBS    += ${BLACS_HOME}/lib/blacs_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsF77init_MPI-LINUX-0.a
  LIBS    += ${BLACS_HOME}/lib/blacsCinit_MPI-LINUX-0.a
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST = 
