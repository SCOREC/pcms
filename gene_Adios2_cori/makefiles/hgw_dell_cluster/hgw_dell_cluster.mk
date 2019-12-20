###########################################################################
### architecture dependent GENE makefile for HGW Dell Blade Cluster     ###
### accessible from, e.g., work20.ipp-hgw.mpg.de                        ###
###########################################################################

###########################################################################
### IMPORTANT: DO NOT COMPILE ON THE LOGIN NODE                         ###
### LAUNCH AN INTERACTIVE SHELL: mksge -i --nox -n <# of processors>    ###
###                                                                     ###
### and load the following modules:                                     ###
###                                                                     ###
### module load slepc mkl fftw       					###
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = intel

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP =

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
#MPRUN = mpiexec -n $(N_PES) ./$(EXEC)
MPRUN = mpirun -np $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################


# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

PRECISION= double
DEBUG= no
SLEPC= yes
SCALAPACK = yes
OPENMP = no
USE_PERFLIB = none
FUTILS = no
MKLVERSION = 11.0
MB_PER_CORE = 2000
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no
COMPILE_MPI_MOD = yes

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
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
#fftw is linked in petsc
 FFTW_HOME = $(PETSC_DIR)
 INCPATHS += -I$(FFTW_HOME)/include
 LIBS += -L$(FFTW_HOME)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# PETSC_ARCH =
# PETSC_DIR =
# SLEPC_DIR =
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
#Insert only BLAS library
  LIBS += $(MKL_BLAS_LINKLINE)
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
