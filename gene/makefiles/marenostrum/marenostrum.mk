###########################################################################
### architecture dependent GENE makefile                                ###
### MareNostrum IV, BCS-CNS                                             ###
### https://www.bsc.es/user-support/mn4.ph                              ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
###                                                                     ###
MODULELIST_COMMON = intel petsc/3.7.6-complex slepc/3.7.4-complex fftw
###                                                                     ###
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

ARCHIVE = ar r

MPRUN = srun ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = fftw

PRECISION= double

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
endif

ifeq ($(SCALAPACK),yes)
 LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
else
 LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 HDF5PATH =
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST =

#
