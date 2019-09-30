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
#uncomment the following lines on CRAY machine
#MPFC = ftn
#MPCC = cc
#COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP = Haswell

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
DEBUG= yes

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
FUTILS = yes

PRODRUN = no
ADIOS = yes
COUPLE_XGC = yes

#Provide an upper (RAM) memory limit in MB per core for the code internal
#performance optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE= 3000
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
INCPATHS = -I/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/openblas-0.3.5-7miavkpfij3m35htbo6sabiobj5euoeb/include

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/openblas-0.3.5-7miavkpfij3m35htbo6sabiobj5euoeb/lib -lopenblas

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
 INCPATHS += -I/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/fftw-3.3.8-htqxzzkouh6varwbleedft67f5i6nsxj/include
 LIBS += -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/fftw-3.3.8-htqxzzkouh6varwbleedft67f5i6nsxj/lib/ -lfftw3
endif

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
  LIBS += -L/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/netlib-scalapack-2.0.2-n6xl7lyau3awm2qalbeyvszxm3kijg6k/lib/ -lscalapack -lopenblas
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
 HDF5PATH = /opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/hdf5-1.8.21-nqg53a4knig6ipcmqegnhpy2nnqbgzhy
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz
 #HDF5VAR = /static

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(HDF5PATH)/include/ -I$(FUTILSDIR)
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
