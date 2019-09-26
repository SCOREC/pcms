###########################################################################
### architecture dependent GENE makefile                                ###
### TOK Cluster           		     				###
### https://wiki.mpcdf.mpg.de/ipphpc/index.php/TOK_batch_cluster_in_Garching
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### List required environment settings/modules here (optional):         ###
###                                                                     ###
### Whitespace separated module list for Intel Compiler                 ###
MODULELIST_INTEL=intel mkl impi fftw-serial hdf5-mpi \
		       petsc-complex slepc-complex
###                                                                     ###
### Whitespace separated module list for GNU/GCC Compiler               ###
MODULELIST_GNU=gcc/8 impi mkl fftw-serial hdf5-mpi \
		       petsc-complex slepc-complex
###                                                                     ###
### Whitespace separated module list for PGI Compiler                   ###
MODULELIST_PGI=
###                                                                     ###
## Common modules                                                       ###
MODULELIST_COMMON=git idl gnuplot anaconda texlive
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler
COMPILER = intel
#COMPILER=gnu

#set chip for proper chip optimization flags (optional)
CHIP = SkyLake

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = srun -n $(N_PES) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

PRODRUN = yes

FFTLIB = mkl
#FFTLIB = fftw

PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= yes

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

OPENMP = no

# * possible choice distributed with GENE are USE_PERFLIB=FR and HT
USE_PERFLIB = none

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = yes

#Provide an upper (RAM) memory limit in MB per core for the code internal
MB_PER_CORE=2000

INCLUDE_SYMBOLS = yes
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
INCPATHS =

#LIBRARIES AND LIBFLAGS
LIBS =

#FFT LIBRARY
#MKL_DYNAMIC = no
ifeq ($(FFTLIB),mkl)
 INCPATHS += -I$(MKLROOT)/include
 LIBS +=  $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Wl,-rpath,$(FFTW_HOME)/lib
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH ?=
 PETSC_DIR ?=/tokp/work/tbg/soft/petsc-complex-$(COMPILER)
 SLEPC_DIR ?=/tokp/work/tbg/soft/slepc-complex-$(COMPILER)
endif

ifeq ($(SCALAPACK),yes)
 LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
else
 LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif


ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =
