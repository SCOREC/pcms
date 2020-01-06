###########################################################################
### architecture dependent GENE makefile for EPFL BlueGene/Q Lemanicus  ###
###									###
### load the following modules in e.g. .bashrc				###
###									###
### module load hdf5/1.8.8 fftw						###
### module load xlf/14.1  wrapper/xl					###
### module load essl							###
### module load scalapack lapack					###
###									###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
FFTLIB = fftw

#precison double
PRECISION = double

OPENMP = no

#custom installation of PETSC
SLEPC= yes

USE_PERFLIB = none

DEBUG = no

#scalapack is available
SCALAPACK = yes

#Futils for H5 and CHEASE input
FUTILS = yes

#memory per core
MB_PER_CORE=1000

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
BGQ_SYS=/bgsys/drivers/ppcfloor/comm/xl
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR) -I$(BGQ_SYS)/include

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS =-L$(BGQ_SYS)/lib

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=-I/bgsys/ibm_essl/prod/opt/ibmmath/include
 LIBS += -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64 -lesslbg
endif
ifeq ($(FFTLIB),fftw)
 FFTW_DIR =/bgsys/local/fftw/3.3.3/
 INCPATHS += -I$(FFTW_DIR)/include
 LIBS += -L$(FFTW_DIR)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  PETSC_DIR=/home/merlo/soft/petsc-3.5.2
  SLEPC_DIR=/home/merlo/soft/slepc-3.5.3
  PETSC_ARCH=BGQ
  ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
   include $(SLEPC_DIR)/conf/slepc_common
  else
   include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
  endif
  INCPATHS += $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
  LIBS += $(SLEPC_LIB) $(PETSC_LIB)
endif

ifeq ($(SCALAPACK),yes)
  LIBS += -L/bgsys/local/lib64 -lscalapack -llapack
  LIBS += -L/bgsys/local/lib64 -lesslbg
endif

ifeq ($(USE_PERFLIB),perf)
 LIBS += -L/usr/local/lib
endif

ifeq ($(FUTILS),yes)
 FUTILSDIR = /home/merlo/soft/futils/src
 HDF5PATH = /bgsys/local/hdf5/1.8.8

 HDF5LIBPATH  = -L$(HDF5PATH)/lib -L/bgsys/local/lib64
 INCPATHS +=  -I$(FUTILSDIR) -I$(HDF5PATH)/include

 HDF5LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz
 LIBS += $(HDF5LIBPATH) $(HDF5LIBS) -L$(FUTILSDIR)
endif

###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpixlf90_r

CC = mpixlc_r

ARCHIVE = ar r

#FORTAN COMPILER FLAGS
ifeq ($(PRODRUN),yes)
 OPTLEVEL = 4
 IPA= -qipa
else
 OPTLEVEL = 0
 IPA= -qnoipa
endif

ifeq ($(DEBUG),yes)
 OPT= -qcheck -g -qsigtrap=xl__trcedump -qsource
else
 OPT= -qnocheck -qnoddim -qunroll
endif

FFLAGS = $(OPT) $(IPA) -O$(OPTLEVEL) -qmaxmem=-1 -qmoddir=$(OBJDIR)
FFLAGS += -qflag=I:I -qextname=flush -qautodbl=dbl4 -qarch=qp -qtune=qp -qessl -q64

#can get rid of this if xlf2003 is used
FFLAGS +=-qxlf2003=polymorphic


##------- mpirun
run:	$(EXECDIR)/$(EXEC)
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS)\
        runjob --np ${NP} --ranks-per-node 16 : $(EXEC)

##############################################################################