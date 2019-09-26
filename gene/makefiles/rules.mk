##################################################################
# * contains rules for compiling and linking GENE                #
# * sets FFLAGS and precompiler switches based on definitions in #
#   in MACHMKFILE and $(COMPILER).def                            #
##################################################################

include $(MACHMKFILE)
ifneq ($(COMPILER),)
#remove if once all MACHMKFILEs set COMPILER
 include $(BASEDIR)/makefiles/compilers/$(COMPILER).def
 include $(BASEDIR)/makefiles/compilers/mkllinkline.def
endif

##################################################################
# set defaults for system commands
MPCC ?= $(CC)
MPLD ?= $(MPFC)
LD ?= $(FC)
ARCHIVE ?= ar r

ifeq ($(MPRUN),)
 cmd=$(shell basename "`which mpiexec 2>/dev/null`")
 ifneq ($(cmd),mpiexec)
  cmd=$(shell basename "`which mpirun 2>/dev/null`")
  ifneq ($(cmd),mpirun)
   cmd=$(shell basename "`which aprun 2>/dev/null`")
   ifneq ($(cmd),aprun)
     cmd=$(shell basename "`which runjob 2>/dev/null`")
    ifneq ($(cmd),runjob)
     $(warning MPRUN not set in $(MACHINE).mk)
    else
     MPRUN = (runjob -n $(N_PES) --ranks-per-node 16 : ./$(EXEC))
    endif
   else
    MPRUN = (aprun -n $(N_PES) ./$(EXEC))
   endif
  else
   MPRUN = (mpirun -np $(N_PES) ./$(EXEC))
  endif
 else
  MPRUN = (mpiexec -n $(N_PES) ./$(EXEC))
 endif
endif



##################################################################
#backward compatibility
ifeq ($(MODDIR_FFLAGS),)
 MODDIR_FFLAGS = $(MODDIR)
endif

#set common INCPATHS
INCPATHS += -I$(OBJDIR) -I. -I$(SRCDIR)

#set FFLAGS/CFLAGS/PREPROC etc.
PREPROC += $(F2003_MISSING) $(GITDEF)

ifneq ($(MB_PER_CORE),)
 PREPROC += -D'MB_PER_CORE=$(MB_PER_CORE)'
endif

FFLAGS += $(MODDIR_FFLAGS)
FFLAGS += $(SET_FORTRAN_STANDARD)

ifneq ($(strip $(DEBUG)),no)
# DEBUG is yes or noopt
 INCLUDE_SYMBOLS=yes
 ifeq ($(DEBUG),noopt)
  FFLAGS += $(DEBUG_TRACEBACK) -O0
 endif
 ifeq ($(DEBUG),yes)
  FFLAGS += $(DEBUG_FFLAGS)
  ifneq ($(COMPILER),cray)
    FFLAGS += -O0
  endif
 endif
 ifeq ($(DEBUG),check)
  FFLAGS += $(DEBUG_CHECK_FFLAGS) -O0
 endif
 ifeq ($(DEBUG),sanitize)
  FFLAGS += -g $(DEBUG_SANITIZE) -O0
  LDFLAGS+= -g $(DEBUG_SANITIZE) -O0
 endif
else
 ifeq ($(PRODRUN),yes)
  FFLAGS += $(PRODRUN_FFLAGS)
 else
  FFLAGS += $(OPT_FFLAGS)
 endif
endif

ifeq ($(WITH_DEBUG_OUTPUT),yes)
 PREPROC += -DWITH_DEBUG_OUTPUT
endif

ifeq ($(INCLUDE_SYMBOLS),yes)
 FFLAGS += $(SYMBOL_FFLAGS)
 CFLAGS += $(SYMBOL_CFLAGS)
endif
ifeq ($(INCLUDE_SYMBOLS),trace)
 FFLAGS += $(SYMBOL_FFLAGS) $(DEBUG_TRACEBACK)
 CFLAGS += $(SYMBOL_FFLAGS) $(DEBUG_TRACEBACK)
endif

ifeq ($(COMPILER_REPORTS),yes)
 FFLAGS += $(REPORT_FFLAGS)
endif

ifeq ($(FFTLIB),essl)
 PREPROC += -DWITHESSL
endif

ifeq ($(SLEPC),yes)
 #include SLEPC_DIR only if not already done in machine makefile
 ifeq ($(PETSC_FC_INCLUDES),)
   ifeq (,$(wildcard $(SLEPC_DIR)/lib/slepc/conf/slepc_common))
    ifeq (,$(wildcard $(SLEPC_DIR)/conf/slepc_common))
     $(warning CANNOT INCLUDE SLEPC_COMMON)
    else
     include $(SLEPC_DIR)/conf/slepc_common
    endif
   else
    include $(SLEPC_DIR)/lib/slepc/conf/slepc_common
   endif
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
 endif
 PREPROC += -DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
 LIBS += $(LIB_SCALAPACK)
 PREPROC += -DWITHSCAL
else
 LIBS += $(LIB_BLAS)
endif

ifeq ($(FUTILS),yes)
 ifeq ($(HDF_LIBS),)
  FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
  LIBS += $(HDF5_LIBPATH) -L$(FUTILSDIR) -lfutils $(HDF5_LIBS)
  INCPATHS += -I$(FUTILSDIR)
 endif
 PREPROC +=  -DWITHFUTILS
endif

ifeq ($(HAC),yes)
 PREPROC +=  -DWITHHAC
 PREPROC += -DHLST_HAC_USE_REAL=0 -DHLST_HAC_USE_CMSK=1 -DHLST_HAC_USE_ADIM=6
 PREPROC += -DHLST_HAC_EXTRA_METADATA=1
endif

ifeq ($(OUTPUT),yes)
 PREPROC+= -DWITH_OUTPUT
endif

ifeq ($(OPENMP),yes)
 FFLAGS += $(OPENMP_FFLAGS)
 PREPROC += -DWITHOMP
else
 OMP_NUM_THREADS=1
endif

ifeq ($(OMP_NUM_THREADS),)
 OMP_NUM_THREADS=1
endif

ifeq ($(FFTLIB),mkl)
 MKL_FFLAGS = $(FFLAGS)
 #all flags except double precision
 ifeq ($(MKLINCLUDE),)
  ifneq ($(MKLROOT),)
   MKLINCLUDE = $(MKLROOT)/include
  else
   ifneq ($(MKL_HOME),)
    MKLINCLUDE = $(MKL_HOME)/include
   endif
  endif
 endif
endif

ifeq ($(PRECISION),double)
 FFLAGS += $(DOUBLE_PRECISION_FFLAGS)
 PREPROC += -DDOUBLE_PREC $(DOUBLE_PRECISION_PREPREC)
 LDFLAGS += $(DOUBLE_PRECISION_LDFLAGS)
endif


## add -WF, flag to PREPROC for IBM preprocessor
ifeq ($(findstring xlf,$(MPFC)),xlf)
 PREPROC += -DAIX
 CPREPROC += -Wp,-DPOWER
 PREFIX = -WF,
 FPREPROC = $(addprefix $(PREFIX),$(PREPROC))
else
 FPREPROC = $(PREPROC)
 CPREPROC = $(PREPROC)
endif

NOOPTFLAGS ?= $(MODDIR_FFLAGS) -O0 $(DOUBLE_PRECISION_FFLAGS)
LDFLAGS ?= $(FFLAGS)

#FORCHECK related:
empty:=
space:= $(empty) $(empty)
incflag:=$(space)-I
colon:=:
FCK_INCPATHS = $(PPDIR)$(subst $(incflag),$(colon),$(space)$(strip $(INCPATHS)))
FCKCNF = $(FCKDIR)/share/forcheck/f03.cnf


##################################################################

include $(FFILES)

##################################################################

##########################################################################
### VARIABLES
##########################################################################

#Compiler specific MODULELISTs can be specified by the architecture dependent
#makefiles and returned with get_required_modules (see target below)
MODULELIST += $(MODULELIST_$(shell echo $(COMPILER) | tr a-z A-Z))
MODULELIST += $(MODULELIST_COMMON)


##############################################################################
### RULES
##############################################################################

###########################################################################
### Compiling rules (FLAGS should be set before)                        ###
###########################################################################
#note: adding the @ sign in front of $(CC) and $(MPFC) will suppress the
#      full compiler call which increases readability
ifeq ($(GENE_COMPILE_VERBOSE),1)
 ATCMD=
 NULLEXT=>/dev/null
else
 ifeq ($(GENE_COMPILE_VERBOSE),2)
  ATCMD=@echo "$@: " && time
  NULLEXT=>/dev/null
 else
  ATCMD=@
  NULLEXT=
 endif
endif

##################################################################

ifneq ($(MKDIFF),)
 MKDIFFMSG=WARNING: $(MACHMKFILE) differs from makefiles/ version!
else
 MKDIFFMSG=
endif

DIFFBRANCH=$(shell git diff HEAD --numstat $(SRCDIR) 2>/dev/null | wc -l)
ifneq ($(shell git show-ref master 2>/dev/null),)
DIFFMASTER=$(shell git diff master --numstat $(SRCDIR) 2>/dev/null | wc -l)
else
DIFFMASTER="no master found"
endif

##################################################################

# for the mpimod we switch off all flags, only the
# path of where to write the mpi.mod file is given via MODDIR
$(OBJDIR)/mpimod.o: $(SRCDIR)/mpimod.F90
	@echo "--------- Compiling mpi module ------"
	$(ATCMD)$(MPFC) $(INCPATHS) $(KIND_FFLAG) $(MODDIR_FFLAGS) -c -o $@ $<

#same for the MKL interface which is only called if FFTLIB=mkl
#and requires a properly set MKLINCLUDE
$(OBJDIR)/mkl_dfti.o: $(MKLINCLUDE)/mkl_dfti.f90
	@echo $(MPFC) $(MKL_FFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(MKL_FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo $(MPCC) $(CFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPCC) $(CFLAGS) $(INCPATHS) $(CPREPROC) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	@echo $(NVCC) $(CUDAFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(NVCC) $(CUDAFLAGS) $(INCPATHS) $(CPREPROC) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	@echo $(MPFC) $(FFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -c -o $@ $<

$(PPDIR)/%.f90: $(SRCDIR)/%.F90
	@echo $(MPFC) $(FPREPROC) $(ONLY_PREPROCESS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(FFLAGS) $(ONLY_PREPROCESS) $(INCPATHS) $(FPREPROC) $< >$@

#special treatment (no optimizations) for files in NOOPTLIST
NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	@echo $(MPFC) $(NOOPTFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

###########################################################################
ifeq ($(FUTILS),yes)
ifeq ($(FUTILSDIR),$(EXTDIR)/$(MACHINE)/futils/src)
export FUTILSDIR MPFC FC FFLAGS LDFLAGS MPCC CC CFLAGS HDF5PATH ARCHIVE HDF5_LIBS HDF5VAR COMPILER OPT SET_F2003_STANDARD PREPROC_FLAG
$(FUTILSDIR)/libfutils.a:
	@make -C $(EXTDIR) OPT="$(OPT_FFLAGS)" futils
endif
endif

###########################################################################
ifeq ($(HAC),yes)
 LIBS += $(ADIOS_LIB)
 INCPATHS += $(ADIOS_INC)
endif

###########################################################################

.PHONY: get_required_modules


###########################################################################
##include $(MACHMKFILE)
##ifneq ($(COMPILER),)
###remove if once all MACHMKFILEs set COMPILER
## include $(BASEDIR)/makefiles/compilers/$(COMPILER).def
##endif
#------- Create executable ----- #
$(EXECDIR)/$(EXEC): $(OBJLIST) $(COBJLIST) $(OBJDIR)/gene.o
	@echo "Linking "
	$(MPLD) -o $@ $^ $(LDFLAGS) $(LIBS)
	@echo "File modifications (src) w.r.t. local branch HEAD: $(DIFFBRANCH)"
	@echo "File modifications (src) w.r.t. local master HEAD: $(DIFFMASTER)"
	@echo "$(MKDIFFMSG)"

##################################################################
#------- Create GENE library ---- #
$(EXECDIR)/$(LIBN).a: $(LIBELEMENTS)
	$(ARCHIVE) $@ $(LIBELEMENTS)
	ranlib $@
	@echo ""
	@echo "--- Link $(LIBN) to your program with ---"
	@echo "-L$(ABSBASEDIR)/bin -l$(LIBN2) $(LIBS)"
	@echo ""

$(OBJDIR)/test_$(LIBN).o: $(SRCDIR)/test_$(LIBN).F90
	$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -c -o $@ $<

$(EXECDIR)/test_$(LIBN): $(OBJDIR)/test_$(LIBN).o $(EXECDIR)/$(LIBN).a
	$(MPLD) -o $@ $(LIBELEMENTS) $(OBJDIR)/test_$(LIBN).o \
	-L$(ABSBASEDIR)/$(EXECDIR) -l$(LIBN2) $(LDFLAGS) $(LIBS)

##################################################################
ULSTACK = (ulimit -s unlimited 2>&1) > /dev/null

#------- Execute GENE ----------- #
run:	$(EXECDIR)/$(EXEC)
	@echo "Calling GENE with $(N_PES) MPI and $(OMP_NUM_THREADS) OpenMP thread(s)"
	$(ATCMD)$(ULSTACK);\
        cd $(RUNDIR);\
        $(MPRUN)

##################################################################
# Miscellaneous
show_srcnames:
	@echo $(F90FULLSRC)
	$(ATCMD)for x in $(F90FULLSRC); do echo $$x >> flist.txt; done
	@echo "--> check ./flist.txt for full list"

# run the preprocessor
show_pp_srcnames:
	@echo $(F90PRENAMES)
	$(ATCMD)for x in $(F90PRENAMES); do echo $$x >> flist.txt; done
	@echo "--> check ./flist.txt for full list"

show_objlist:
	$(ATCMD)rm -f objlist.txt cobjlist.txt
	$(ATCMD)for x in $(OBJLIST); do echo $$x >> objlist.txt; done
	$(ATCMD)for x in $(COBJLIST); do echo $$x >> cobjlist.txt; done
	@echo "--> check ./objlist.txt and ./cobjlist.txt for full list"

preproc: $(F90PRENAMES) $(PPDIR)/gene.f90

forcheck: preproc
	forchk -I $(FCK_INCPATHS) -l mylistfile -ninf \
	-ff -i4 -dp -allc -declare -nshsrc -ancmpl -anref -shmodvar \
	$(F90PRENAMES) $(PPDIR)/gene.f90 -library $(FCKDIR)/share/forcheck/MPI_2.flb

liblinkline:
	@echo "-L$(ABSBASEDIR)/bin -l$(LIBN2) $(LIBS)"

get_required_modules:
	@echo "module load $(MODULELIST)"
