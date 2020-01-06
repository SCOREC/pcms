###########################################################################
### architecture dependent GENE makefile, PPB-110                       ###
###########################################################################
### Load the following module in e.g. ~/.bashrc				###
###  									###
### module load mpich/3.1.3-intel13.1                                   ###
### module laod intel_comp/13.1                                         ###
### module load hdf5/1.8.12_intel13.1                                   ###
### module load slurm                                                   ###
###									###
### We currently rely on custom installation of                         ###
###  - fftw3                                                            ###
###  - petsc (3.6.1)                                                    ###
###  - slepc (3.6.1)                                                    ###
###                                                                     ###
### built in /home/merlo/soft                                           ###
###                                                                     ###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
### No parallel filesystem available, so no parallel IO can be done     ###
### This means switch off any checkpointing and s_checkpointing         ###
### The PVFS filesystem is available but unstable. Feel free to use it  ###
### If you really need checkpoints set many_chpt=T                      ###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###									###
###########################################################################


###########################################################################
### SWITCHES                                                            ###
###########################################################################
COMPILER = intel
CHIP= Haswell

ARCHIVE = ar r

FFTLIB = fftw
PRECISION= double

OPENMP = no

MPRUN = export OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
        export MKL_SERIAL=yes;\
        mpiexec -n $(N_PES) $(EXEC)

SCALAPACK = yes

SLEPC= yes

USE_PERFLIB = none

DEBUG =no

FUTILS=yes

#might be changed according to ypour computer perfomrance
#there are currently 32 GB per node, which is a quad core.
#The following assumes no hyperthreading, if you enable it
#might be worth cjheking its not better to lower the number
PREPROC= -D'MB_PER_CORE=6000'
###########################################################################
#   COMPULSORY LIBRARIES                                                  #
###########################################################################
MKLROOT = $(MKL)/../../

ifeq ($(FFTLIB),fftw)
 INCPATHS += -I/home/merlo/soft/fftw-3.3.4/include
 LIBS += -L/home/merlo/soft/fftw-3.3.4/lib64
 ifeq ($(PRECISION),double)
        LIBS += -lfftw3
 else
        LIBS += -lfftw3f
 endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)                         #
###########################################################################

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH =  $(HDF5)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz
 INCPATHS += -I$(FUTILSDIR) -I$(HDF5PATH)/include

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
endif

ifeq ($(SLEPC),yes)
 PETSC_ARCH=arch-linux2-c-opt
 PETSC_DIR=/home/merlo/soft/petsc-3.7.2/
 SLEPC_DIR=/home/merlo/soft/slepc-3.7.2/

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
