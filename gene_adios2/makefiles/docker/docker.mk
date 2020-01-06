###########################################################################
### Used with Docker    		     									###
###########################################################################


###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = gnu
CHIP =
ARCHIVE = ar r


#MPI command to run the executable $(EXEC) with $N_PES MPI processes -- I don't use this
MPRUN = mpiexec -n $(N_PES) ./$(EXEC)


###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk -- (set to: mkl, fftw or essl)
FFTLIB = fftw

#double precision should be default; single precision, however, might #be sufficient and much faster
PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= yes

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again. (Currently, pure MPI is most likely the best choice)
OPENMP = yes

#performance profiling library
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = yes


# ADIOS
ADIOS = yes

# Coupling with XGC
COUPLE = yes
COUPLE_CODE = $(WDM_CODE)/ecp-wdm/
COUPLING_FILE = "../../src/xgc-coupled.json"

# GPTL
WITH_GPTL = yes
PERFMOD_DIR=/usr/local/gptl-install


#Provide an upper (RAM) memory limit in MB per core for the code internal performance optimization (w.r.t. MPI mapping & alternative routines)
#Note: the system itself often requires a significant fraction which should be taken into account
MB_PER_CORE=1024

#include symbols or compiler reports for debugging (will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################


#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
USR_DIR = /usr
INCPATHS = -I${USR_DIR}/include -I${LOCAL_DIR}/include
LIBS += -L${USR_DIR}/lib  -L${LOCAL_DIR}/lib


#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

#optional parameters to tune mkl link line used for FFT/SCALAPACK
#MKLVERSION = 11.3
#MKLVERSION = 2017
#MKL_DYNAMIC = no

ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS += $(MKL_LINKLINE)
#uses link line from makefile/compilers/mkllinkline.def
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS +=
 #LIBS += -lfftw3 -lfftw3f
 LIBS += -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if those variables are not
# set by modules (make sure, you use the *complex* versions)
 PETSC_ARCH = 
 PETSC_DIR = ${LOCAL_DIR}/petslep-complex
 SLEPC_DIR = ${LOCAL_DIR}/petslep-complex

 include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS += -lscalapack
# uncomment if intel mkl is available; uses link line
# from makefile/compilers/mkllinkline.def:
# LIBS += $(MKL_SCALAPACK_LINKLINE)
else
# uncomment if intel mkl is available; uses link line
# from makefile/compilers/mkllinkline.def:
# LIBS += $(MKL_BLAS_LINKLINE)
endif


ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = ${LOCAL_DIR}
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif


ifeq ($(ADIOS),yes)
  ADIOS_LIB = $(shell adios_config -l -f)
  ADIOS_INC = $(shell adios_config -c -f)
  ADIOS_DIR = $(shell adios_config -d)
endif

ifeq ($(WITH_GPTL),yes)
  INCPATHS += -I$(PERFMOD_DIR)/include
  LIBS += -L$(PERFMOD_DIR)/lib -lgptl
endif


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
