###########################################################################
### architecture dependent GENE makefile                                ###
### MareNostrum IV, BCS-CNS                                             ###
### https://www.bsc.es/user-support/mn4.ph                              ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
###                                                                     ###
### module load petsc/3.7.6-complex slepc/3.7.4-complex                 ###
### module load fftw                                                    ###
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

COMPILER = intel
CHIP =  SkyLake

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
#MPRUN = mpiexec -n $(N_PES) ./$(EXEC)
#MPRUN = mpirun -np $(N_PES) ./$(EXEC)
MPRUN = srun ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################


# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl or fftw
FFTLIB = fftw

#double precision should be default; single precision, however, might
#be sufficient and much faster
PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no


SLEPC= yes

SCALAPACK = yes

OPENMP = no

USE_PERFLIB = none

FUTILS = no

MB_PER_CORE=1600

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
#optional parameters to tune mkl link line used for FFT/SCALAPACK
MKLVERSION = 2017
#MKL_DYNAMIC = no

ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS += $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += $(FFTW_INCL)
 LIBS += $(FFTW_LIBS)
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
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

#
