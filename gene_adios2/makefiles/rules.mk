#################################################################
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
 cmd=$(shell basename $(shell which mpiexec))
 ifneq ($(cmd),mpiexec)
  cmd=$(shell basename $(shell which mpirun))
  ifneq ($(cmd),mpirun)
   cmd=$(shell basename $(shell which aprun))
   ifneq ($(cmd),aprun)
     cmd=$(shell basename $(shell which runjob))
    ifneq ($(cmd),runjob)
     $(error MPRUN not set in $(MACHINE).mk)
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

##################################################################

include $(FFILES)

##################################################################
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
  FFLAGS += $(DEBUG_FFLAGS) -O0
 endif
else
 ifeq ($(PRODRUN),yes)
  FFLAGS += $(PRODRUN_FFLAGS)
 else
  FFLAGS += $(OPT_FFLAGS)
 endif
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
 PREPROC += -DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
 PREPROC += -DWITHSCAL
endif

ifeq ($(FUTILS),yes)
 PREPROC +=  -DWITHFUTILS
endif


ifeq ($(ADIOS),yes)
    PREPROC += -DADIOS
    ifeq ($(ADIOS_KY0),yes)
        PREPROC += -DADIOS_KY0
    endif
endif

ifeq ($(ADIOS2),yes)
 PREPROC += -DADIOS2
 ADIOS2_LIB = $(shell $(ADIOS2_DIR)/adios2-config --fortran-libs)
 ADIOS2_INC = $(shell $(ADIOS2_DIR)/adios2-config --fortran-flags)
 LIBS += $(ADIOS2_LIB)
 INCPATHS += $(ADIOS2_INC)
endif

ifeq ($(COUPLE_XGC),yes)
 COUPLE=yes
 ifneq ($(FUTILS),yes)
  $(error Coupling must be built with FUTILS to read XGC grid)
 endif

 ifneq ($(ADIOS2),yes)
  $(error Coupling must be built with ADIOS. ADIOS must be yes)
 endif
 PREPROC += -DCOUPLE_XGC -DGENE_SIDE -DCOUPLE
endif

ifeq ($(INIT_XGC),yes)
 PREPROC += -DINIT_XGC
endif


ifeq ($(WITH_GPTL),yes)
  PREPROC += -DWITH_GPTL
endif

ifeq ($(OUTPUT),yes)
 PREPROC+= -DOUTPUT
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

#FFLAGS += -ffree-line-length-0


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

DIFFBRANCH=$(shell git diff HEAD --numstat $(SRCDIR) | wc -l 2>/dev/null)
ifneq ($(shell git show-ref master),)
DIFFMASTER=$(shell git diff master --numstat $(SRCDIR) | wc -l 2>/dev/null)
else
DIFFMASTER="no master found"
endif

##################################################################

# for the mpimod we switch off all flags, only the
# path of where to write the mpi.mod file is given via MODDIR
$(OBJDIR)/mpimod.o: $(SRCDIR)/mpimod.F90
	@echo "--------- Compiling mpi module ------"
	$(ATCMD)$(MPFC) $(MODDIR_FFLAGS) -c -o $@ $<

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



##############################################################################
ifeq ($(ADIOS),yes)
 LIBS += $(shell adios_config -l -f)
 INCPATHS += $(shell adios_config -c -f)
 ADIOSBASE = $(shell adios_config -d)
endif

$(OBJDIR)/adios_io.o: $(SRCDIR)/adios_io.F90 $(XMLFILE)
	@echo gpp.py $(notdir $(XMLFILE))
	$(ATCMD)cd $(SRCDIR); python $(ADIOSBASE)/bin/gpp.py $(notdir $(XMLFILE)); cd $(RUNDIR); $(NULLEXT)
	@echo $(MPFC) $(FFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -c -o $@ $<

##############################################################################
# couplng file
ifeq ($(ADIOS),yes)
$(OBJDIR)/coupling_core_gene.o: $(SRCDIR)/coupling_core_gene.F90 $(XMLFILE)
	@echo gpp.py $(notdir $(XMLFILE))
	$(ATCMD)cd $(SRCDIR); python $(ADIOSBASE)/bin/gpp.py $(notdir $(XMLFILE)); cd $(RUNDIR); $(NULLEXT)
	@echo $(MPFC) $(FFLAGS) $(notdir $<) $(NULLEXT)
	$(ATCMD)$(MPFC) $(FFLAGS) $(INCPATHS) $(FPREPROC) -c -o $@ $<
endif

##############################################################################


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

###########################################################################
##include $(MACHMKFILE)
##ifneq ($(COMPILER),)
###remove if once all MACHMKFILEs set COMPILER
## include $(BASEDIR)/makefiles/compilers/$(COMPILER).def
##endif


#.PHONY:

#------- Create executable ----- #
#main: coupler basic
#basic: $(EXECDIR)/$(EXEC)

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


coupler:
	@echo "------  coupler rule - you should not be here, something is wrong --------"
	#rsync -a --exclude="${COUPLE_DIR}/install" $(COUPLE_DIR)/ coupler
	#$(MAKE) -C coupler FFLAGS="" FC="${MPFC}" PARAMS=${COUPLING_FILE} PREFIX=install lib
	#rsync -a coupler/*.mod coupler/libcoupler.a coupler/modules/libcmodules.a $(OBJDIR)
