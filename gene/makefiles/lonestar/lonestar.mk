###########################################################################
### architecture dependent GENE makefiles, Lonestar5 @ TACC             ###
### https://portal.tacc.utexas.edu/user-guides/lonestar5                ###
###########################################################################
#									  #
# Stack is too old, so All is installed with custom modules               #
# only gcc is available and functional so do this:			  #
#									  #
# module use /home1/05447/merlo/soft/modules				  #
# module load gcc/7.3.0							  #
# module unload intel							  #
# module load fftw3							  #
#									  #
###########################################################################
### SWITCHES                                                            ###
###########################################################################
COMPILER=gnu

CHIP=Haswell

FFTLIB = fftw

PRECISION= double

OPENMP = no

SLEPC = yes

SCALAPACK = no

DEBUG = no

INCLUDE_SYMBOLS = yes

COMPILER_REPORTS=no

USE_PERFLIB = none

FUTILS = no

PRODRUN = yes

ARCHIVE =ar r

#memory per core
#need to check how much the OS takes. is 64GB for 48
MB_PER_CORE=2600

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

ifeq ($(FFTLIB),fftw)
   FFTW_HOME=$(TACC_FFTW3_DIR)
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 #-lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################
SCALAPACK_HOME=/home1/05447/merlo/soft/scalapack_installer/install/lib
LAPACK_HOME=/home1/05447/merlo/soft/scalapack_installer/build/lapack-3.8.0

LIBS +=  \
        -L$(SCALAPACK_HOME) -lscalapack -lreflapack -lrefblas

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   PETSC_ARCH ?=
   PETSC_DIR = /home1/05447/merlo/soft/PETSC/petsc3.9.2-gcc7.3-complex
   SLEPC_DIR = /home1/05447/merlo/soft/SLEPC/slepc3.9.1-gcc7.3
   ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
    include $(SLEPC_DIR)/conf/slepc_common
   else
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
   endif

   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
#  LIBS += $(MKL_SCALAPACK_LINKLINE)
else
#  LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(TACC_HDF5_DIR)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR) -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif


###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################

ifeq ($(USE_PERFLIB),perf)
# PREPROC += -DWITHPERF=1
# everything is set in the compiler specific makefiles
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

#LDFLAGS += -Xlinker -M

###########################################################################
### Running                                                             ###
###########################################################################
MPRUN = export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=no;\
	ibrun -n $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=