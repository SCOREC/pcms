###########################################################################
### architecture dependent GENE makefile for IRENE/TGCC                 ###
### http://www-hpc.cea.fr/en/complexe/tgcc-doc-util-private.htm         ###
###########################################################################
### NOTE: We currently experience issues with the latest OPENMPI version ##
### installed on IRENE at the beginning of FEBRUARY 2019 which prevent  ###
### parameter scans in GENE.                                            ###
###                                                                     ###
### You may therefore want to use intelmpi instead:                     ###
### module unload python3 mpi/openmpi                                   ###
### module load mpi/intelmpi python3                                    ###
###                                                                     ###
###########################################################################
### List required environment settings/modules here (optional):         ###
###                                                                     ###
### Whitespace separated module list for Intel Compiler                 ###
MODULELIST_INTEL= intel mpi flavor/hdf5/parallel hdf5
###                                                                     ###
### Whitespace separated module list for GNU/GCC Compiler               ###
MODULELIST_GNU= gnu mpi
###                                                                     ###
### Whitespace separated module list for PGI Compiler                   ###
MODULELIST_PGI=
###                                                                     ###
### Whitespace separated module list for CRAY Compiler                  ###
MODULELIST_CRAY=
###                                                                     ###
### COMMON modules, e.g., for diagnostics, documentation etc.           ###
MODULELIST_COMMON = idl python3 texlive
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = intel

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP = SkyLake

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = ccc_mprun -n $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################


# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl, fftw
FFTLIB = mkl

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
USE_PERFLIB = none

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = no
#Note: For unknown reason, compilation currently freezes at futils.f90 on IRENE

#Provide an upper (RAM) memory limit in MB per core for the code internal
#performance optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=3750
#Note: the system itself often requires a significant fraction
#      which should be taken into account, therefore 3750MB instead of 4GB

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no

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

#optional parameter to tune mkl link line used for FFT/SCALAPACK
#MKL_DYNAMIC = no

ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS += $(MKL_LINKLINE)
#uses link line from makefile/compilers/mkllinkline.def
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS +=
 LIBS +=
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# PETSC_DIR = $(PETSC_ROOT)
# SLEPC_DIR = $(SLEPC_ROOT)
 PETSC_ARCH =
 ifndef I_MPI_ROOT
  PETSC_DIR = /ccc/work/cont005/ra4403/ra4403/software/petsc-complex
  SLEPC_DIR = /ccc/work/cont005/ra4403/ra4403/software/slepc-complex
 else
  PETSC_DIR = /ccc/work/cont005/ra4403/ra4403/software/petsc-complex-intelmpi
  SLEPC_DIR = /ccc/work/cont005/ra4403/ra4403/software/slepc-complex-intelmpi
 endif
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
 LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
else
 LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif


ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5_ROOT)
 HDF5_LIBPATH  = -L$(HDF5_LIBDIR)
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
