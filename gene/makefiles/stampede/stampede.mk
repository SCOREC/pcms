############################################################################
### architecture dependent GENE makefiles, Stampede2 version             ###
### https://portal.tacc.utexas.edu/user-guides/stampede2                 ###
############################################################################
#									   #
# the following modules/environment variables have to be loaded first:	   #
#                                                                          #
# module load fftw                                                         #
# export MACHINE=stampede                                                  #
#									   #
############################################################################
### SWITCHES                                                             ###
############################################################################

COMPILER=intel

CHIP=KnightsLanding
#CHIP=SkyLake

ifeq ($(CHIP),KnightsLanding)
#effective mem. per core: 86 GB/68 ~ 1.26 GB
 MB_PER_CORE=1200
else
 MB_PER_CORE=3700
endif

# FFTLIB - needed in /src/files.mk
#FFTLIB = mkl
FFTLIB = fftw

PRECISION= double

OPENMP = no

SLEPC = yes

SCALAPACK = yes

DEBUG = no

INCLUDE_SYMBOLS = yes

COMPILER_REPORTS=no

USE_PERFLIB = none

# FUTILS and HDF5 are required for some geometry interfaces
FUTILS = no

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
 INCPATHS += -I$(TACC_FFTW3_INC)
 ifeq ($(PRECISION),double)
  LIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpi -lfftw3
 else
  LIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpif -lfftw3f
 endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
#using makefiles/rules.mk block
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
else
  LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
endif

###########################################################################
### Running                                                             ###
###########################################################################
MPRUN = export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=yes;\
	ibrun ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=

#$(OBJDIR)/collisions.o: $(SRCDIR)/collisions.F90
#	@echo $(MPFC) -O1 $(notdir $<) $(NULLEXT)
#	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -O1 -c -o $@ $<
