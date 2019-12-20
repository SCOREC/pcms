############################################################################
### architecture dependent GENE makefiles, marconi version               ###
#https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.1%3A+MARCONI+UserGuide#
############################################################################
#									   #
# the following modules have to be loaded first:		           #
#                                                                          #
#--------------------------------------------------------------------------#
# SELECT ONE PARTITION SPECIFIC ENVIRONMENT FROM THE FOLLOWING CHOICES     #
# - MAKE SURE TO *COMPILE AND RUN* WITH SAME SETTINGS!!!                   #
#                                                                          #
# FOR A1 BROADWELL                                                         #
# module load env-bdw                                                      #
# FOR A2 KNL                                                               #
# module load env-knl profile/knl                                          #
# FOR A3 SKYLAKE                                                           #
# module load env-skl                                                      #
#--------------------------------------------------------------------------#
#                                                                          #
# Necessary libraries:                                                     #
# module load intel intelmpi mkl fftw                                      #
# module load blas lapack scalapack                                        #
# module load szip zlib hdf5                                               #
#                                                                          #
# 									   #
# for pgi or gnu compiler tests:                                           #
# module load pgi                                                          #
# module load gnu                                                          #
#                                                                          #
# for IDL diagnostics:                                                     #
# module load profile/astro idl                                            #
#									   #
############################################################################
### SWITCHES                                                             ###
############################################################################

COMPILER=intel
#COMPILER=pgi

#determine partition automatically by default
ifneq ($(ENV_KNL_HOME),)
 CHIP=KnightsLanding
#effective mem. per core: 86 GB/68 ~ 1.26 GB
 MB_PER_CORE=1200
else
ifneq ($(ENV_SKL_HOME),)
 CHIP=SkyLake
#effective mem. per core: 180 GB/48 = 3.75 GB
 MB_PER_CORE=3700
else
 CHIP=Broadwell
#effective mem. per core: 118 GB/36 = 3.27 GB
 MB_PER_CORE=3200
endif
endif

# FFTLIB - needed in /src/files.mk
#FFTLIB = fftw
FFTLIB = mkl

PRECISION= double
OPENMP = no

SLEPC = yes

SCALAPACK = yes

DEBUG = no

INCLUDE_SYMBOLS = yes

COMPILER_REPORTS=no

USE_PERFLIB = none

FUTILS = yes

PRODRUN = yes

MKLVERSION = 2017
MKL_DYNAMIC = yes

ifeq ($(COMPILER),pgi)
SLEPC=no
FUTILS=no
COMPILE_MPI_MOD = yes
FFTLIB = mkl
endif

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS += -I$(MKLROOT)/include
 LIBS += $(MKL_LINKLINE)
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
   PETSC_ARCH =
   ifeq ($(PRECISION),single)
    PETSC_DIR = /marconi/home/userexternal/tgoerler/soft/petsc-sgl
    SLEPC_DIR = /marconi/home/userexternal/tgoerler/soft/slepc-sgl
   else
    PETSC_DIR = /marconi/home/userexternal/tgoerler/soft/petsc
    SLEPC_DIR = /marconi/home/userexternal/tgoerler/soft/slepc
   endif
   ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
    include $(SLEPC_DIR)/conf/slepc_common
   else
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
   endif

   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
# LIBS += -L$(SCALAPACK_LIB) -lscalapack
# INCPATHS += -I$(SCALAPCK_INC)
else
  LIBS += $(MKL_BLAS_LINKLINE)
endif

#LIBS += -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas
#INCPATHS += -I$(BLAS_INC) -I$(LAPACK_INC)


ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR) -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(ADIOS),yes)
  ADIOS_LIB = $(shell adios_config -l -f)
  ADIOS_INC = $(shell adios_config -c -f)
  ADIOS_DIR = $(shell adios_config -d)
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
	export MKL_SERIAL=yes;\
	mpirun -np $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=

#$(OBJDIR)/collisions.o: $(SRCDIR)/collisions.F90
#	@echo $(MPFC) -O1 $(notdir $<) $(NULLEXT)
#	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -O1 -c -o $@ $<
