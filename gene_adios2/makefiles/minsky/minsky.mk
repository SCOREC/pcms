###########################################################################
### template for architecture dependent GENE makefiles                  ###
### <insert machine name>     		     				###
### <insert a web link to the machine description if available>         ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
###                                                                     ###
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = gnu
COMPILE_MPI_MOD=yes
#uncomment the following lines on CRAY machine
#MPFC = ftn
#MPCC = cc
#COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP =

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = mpiexec -n $(N_PES) ./$(EXEC)
#MPRUN = mpirun -np $(N_PES) ./$(EXEC)
#MPRUN = aprun -n $(N_PES) ./$(EXEC)

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
SLEPC= no

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
MB_PER_CORE=
#Note: the system itself often requires a significant fraction
#      which should be taken into account

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no
WITH_CUTILS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

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
# FFTW_HOME=/mpcdf/soft/RHEL7/MKY/fftw/3.3.5/gcc-6.2/smpi-10.1
 FFTW_HOME=/home/tisd/soft/fftw/3.3.5/gcc-6.2
 INCPATHS += -I$(FFTW_HOME)/include
# LIBS += -L$(FFTW_HOME)/lib -lfftw3
 LIBS += $(FFTW_HOME)/lib/libfftw3.a
 #INCPATHS += -I/usr/include
 #LIBS += -L/usr/lib64 -lfftw3
endif

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
#INCPATHS =

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
SCALAPACK_HOME=/home/tisd/build/scalapack-2.0.2
LAPACK_HOME=/home/tisd/build/lapack-3.6.1
#LIBS += -L/usr/lib64 -lpesslsmp -lblacssmp -lessl \
#	-L/opt/ibm/xlsmp/4.1.5/lib -lxlsmp -Wl,-rpath,/opt/ibm/xlsmp/4.1.5/lib \
#	-Wl,-rpath,/opt/ibm/lib \
#	-L$(SCALAPACK_HOME) -lscalapack \
#	-L$(LAPACK_HOME) -llapack \
#	-L/opt/ibm/xlf/15.1.5/lib -lxl

LIBS += -L/usr/lib64 -lessl \
	-L$(SCALAPACK_HOME) -lscalapack \
	-L$(LAPACK_HOME) -llapack \
	-L/opt/ibm/xlf/15.1.5/lib -lxl

#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

#optional parameters to tune mkl link line used for FFT/SCALAPACK
#MKLVERSION = 11.3
#MKLVERSION = 2017
#MKL_DYNAMIC = no

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if those variables are not
# set by modules (make sure, you use the *complex* versions)
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

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS +=
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
