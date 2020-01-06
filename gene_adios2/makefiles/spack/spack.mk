###########################################################################
### BASIC settings                                                      ###
###########################################################################


# Might update these from spack
COMPILER = gnu
MPFC = f90
MPCC = cc
SCALAPACK := other
COUPLE_CODE = $(WDM_CODE)/ecp-wdm/
COUPLING_FILE = ../../src/xgc-coupled.json
PERFMOD_DIR = /usr/local/gptl-install
CHIP =

# Somewhat more relevant
NVCC = nvcc
ARCHIVE = ar r

# Not really planning to change these
FFTLIB = fftw
PRECISION= double
SLEPC= yes
OPENMP = no
FUTILS = yes
ADIOS = yes
COUPLE = yes
WITH_GPTL = yes

# Less important
MB_PER_CORE= 1500
USE_PERFLIB = none
DIAG_MPI_IO= no
DEBUG= no
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no

# I don't actually use this
MPRUN = aprun -n $(N_PES) ./$(EXEC)


##################################################################################################
# Addtional  steps based on the settings above. Trying to not have to edit this
##################################################################################################

#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = 
LIBS = 

ifeq ($(FFTLIB),mkl)
    INCPATHS +=
    LIBS += $(MKL_LINKLINE)
else ifeq ($(FFTLIB),essl)
    INCPATHS +=
    LIBS +=
else ifeq ($(FFTLIB),fftw)
    LIBS += -lfftw3
endif

ifeq ($(SLEPC),yes)
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
    INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
    LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
    LIBS += -lscalapack
else ifeq ($(SCALAPACK),mkl)
    SCALAPACK := yes	
    LIBS += $(MKL_LINKLINE)
else ifeq ($(SCALAPACK),other)
    SCALAPACK := yes	
endif

ifeq ($(FUTILS),yes)
    FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
    HDF5_LIBS = -L$(FUTILSDIR) -lfutils -lhdf5_fortran -lhdf5 -lz
    LIBS += $(HDF5_LIBS)
    INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(WITH_GPTL),yes)
  INCPATHS += -I$(PERFMOD_DIR)/include
  LIBS += -L$(PERFMOD_DIR)/lib -lgptl
endif
