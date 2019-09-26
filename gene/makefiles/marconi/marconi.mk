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
# THEN load the following libraries via module load                        #
#                                                                          #
# For INTEL Compiler:                                                      #
MODULELIST_INTEL=intel intelmpi mkl fftw blas lapack scalapack szip zlib/1.2.8--gnu--6.1.0 hdf5
#                                                                          #
# For GNU Compiler:                                                        #
MODULELIST_GNU=gnu intel intelmpi mkl fftw blas lapack scalapack szip zlib/1.2.8--gnu--6.1.0 hdf5
#MODULELIST_GNU=gnu/6.1.0 openmpi mkl fftw blas lapack scalapack szip zlib/1.2.8--gnu--6.1.0 hdf5
# 									   #
# For PGI Compiler:				                           #
MODULELIST_PGI=profile/advanced pgi intel intelmpi mkl
#MODULELIST_PGI=profile/advanced pgi openmpi/1-10.3--pgi--16.5 mkl
#                                                                          #
# for IDL diagnostics:                                                     #
MODULELIST_COMMON=profile/astro idl python                                 #
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
ifeq ($(COMPILER),gnu)
 FFTLIB = fftw
else
 FFTLIB = mkl
endif

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

MKL_DYNAMIC = yes

ifeq ($(COMPILER),pgi)
 SLEPC=no
 FUTILS=no
 FFTLIB=mkl
 MPI_BINDING_HOME=/marconi/home/userexternal/tgoerler/soft/impi_binding
 PGI_NAME=impi2017.3.196_pgi16.5
 MPFC=$(MPI_BINDING_HOME)/f90/intel64/bin/mpi$(PGI_NAME)
 MPCC=$(MPI_BINDING_HOME)/c/intel64/bin/mpi$(PGI_NAME)
 INCPATHS = -I$(MPI_BINDING_HOME)/f90/intel64/include/$(PGI_NAME)
 LIBS += -L$(MPI_BINDING_HOME)/f90/intel64/lib -lmpi$(PGI_NAME)
 #add $(MPI_BINDING_HOME)/f90/intel64/lib and $I_MPI_ROOT/lib64 to LD_LIBRARY_PATH
endif

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS +=

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 MKLINCLUDE = $(MKLROOT)/include
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
    PETSC_DIR = /marconi/home/userexternal/tgoerler/soft/petsc-$(CHIP)-sgl
    SLEPC_DIR = /marconi/home/userexternal/tgoerler/soft/slepc-$(CHIP)-sgl
   else
    PETSC_DIR = /marconi/home/userexternal/tgoerler/soft/petsc-$(CHIP)
    SLEPC_DIR = /marconi/home/userexternal/tgoerler/soft/slepc-$(CHIP)
   endif
endif

ifeq ($(SCALAPACK),yes)
  LIB_SCALAPACK += $(MKL_SCALAPACK_LINKLINE)
# LIBS += -L$(SCALAPACK_LIB) -lscalapack
# INCPATHS += -I$(SCALAPCK_INC)
else
  LIB_BLAS += $(MKL_BLAS_LINKLINE)
endif

#LIBS += -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas
#INCPATHS += -I$(BLAS_INC) -I$(LAPACK_INC)


ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz
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
	srun -n $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST=

#$(OBJDIR)/collisions.o: $(SRCDIR)/collisions.F90
#	@echo $(MPFC) -O1 $(notdir $<) $(NULLEXT)
#	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -O1 -c -o $@ $<
