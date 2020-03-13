###############################################################################
###                                                                         ###
###   Master makefile                                                       ###
###                                                                         ###
###   last changed: 2015-May-04                                             ###
###                                                                         ###
###   calls /makefiles/rules.mk                                             ###
###   which includes $MACHINE.mk, /makefiles/compilers/$COMPILER.def and    ###
###                  src/files.mk                                           ###
###                                                                         ###
###   the $MACHINE 'detector' can be found in /makefiles/machine.def'       ###
###                                                                         ###
###############################################################################

###############################################################################
#check if user wants to compile GENE for the present prob directory only
#by searching for a file "make_local" in that directory
pathsearch = $(firstword $(wildcard $(addsuffix /$(1),$(subst :, ,$(PWD)))))
MK_LOCAL:= $(wildcard make_local)
ifeq ($(MK_LOCAL),make_local)
   #executable should be local
   BASEDIR = ..
   EXECDIR  = .
   ABSBASEDIR = $(PWD)/..
else
   pathsearch = $(firstword $(wildcard $(addsuffix /$(1),$(subst :, ,$(PWD)))))
   CHECKBASE := $(call pathsearch,src)
   ifeq ($(CHECKBASE),)
      #call from prob directory
      ABSBASEDIR = $(PWD)/..
      BASEDIR = ..
   else
      #call from gene base directory
      ABSBASEDIR = $(PWD)
      BASEDIR = .
   endif
   EXECDIR  = $(BASEDIR)/bin
endif

##############################################################################
#try to identify MACHINE if not set
ifeq ($(MACHINE),)
 include $(BASEDIR)/makefiles/machine.def
endif

##############################################################################
#set dependent path variables
MACHMKFILE  = $(EXECDIR)/$(MACHINE).mk
MKFILE = $(BASEDIR)/makefiles/rules.mk

RUNDIR = $(PWD)
SRCDIR  = $(BASEDIR)/src
TOOLDIR = $(BASEDIR)/tools
EXTDIR  = $(ABSBASEDIR)/external
OBJDIR  = $(EXECDIR)/obj_$(MACHINE)
PPDIR   = $(EXECDIR)/preproc
FFILES  = $(SRCDIR)/files.mk
RULES	= $(BASEDIR)/makefiles/rules.mk

GITBRANCH = $(shell git rev-parse HEAD 2>/dev/null)
ifneq ($(shell git show-ref master),)
GITMASTER = $(shell git rev-parse master 2>/dev/null)
endif
GITDEF = -D'GIT_BRANCH="$(GITBRANCH)"' -D'GIT_MASTER="$(GITMASTER)"'

MKDIFF = $(shell diff -q $(BASEDIR)/makefiles/$(MACHINE)/$(MACHINE).mk $(MACHMKFILE) | cut -c1)

EXEC = gene_$(MACHINE)
LIBTARGET = $(EXECDIR)/$(LIBN).a

export RUNDIR SRCDIR TOOLDIR EXTDIR OBJDIR PPDIR FFILES RULES
export BASEDIR ABSBASEDIR EXECDIR MACHINE MACHMKFILE
export MKDIFF EXEC LIBTARGET
export GITBRANCH GITMASTER GITDEF

$(shell $(TOOLDIR)/codemods.sh $(SRCDIR))

##############################################################################
#.PHONY targets are always called - independent of the dependencies being
# up-to-date
.PHONY: clean distclean doc update mach_wrapper install run lib all checkpath
##############################################################################

all: $(EXECDIR)/$(EXEC)

checkpath:
	@test -d $(EXECDIR) || mkdir -p $(EXECDIR)
	@test -f $(MACHMKFILE) || cp $(BASEDIR)/makefiles/$(MACHINE)/*.mk $(EXECDIR)
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
	@test -d $(PPDIR)  || mkdir -p $(PPDIR)
ifneq ($(RUNDIR),$(BASEDIR))
#copy submit scripts for calls from prob directory
	@diff $(BASEDIR)/makefiles/$(MACHINE) $(RUNDIR) | \
	 grep -e 'makefiles/$(MACHINE)' | \
	 grep -v -e 'svn' -e '.mk' | awk -F ': ' '{print $$2}' > list
	@for file in `cat list`; do cp $(BASEDIR)/makefiles/$(MACHINE)/$$file $(RUNDIR) ; done
	@rm list
endif
ifeq ($(MK_LOCAL),make_local)
	@echo compiling a GENE executable in $(RUNDIR)
	@test ! -L $(RUNDIR)/$(EXEC) || (rm -f $(RUNDIR)/$(EXEC))
else
	@echo using GENE executable in bin $(RUNDIR)
ifneq ($(RUNDIR),$(ABSBASEDIR))
	@test ! -d $(RUNDIR)/$(MACHINE) || (echo removing $(RUNDIR)/$(MACHINE) directory && \
	   rm -r $(RUNDIR)/$(MACHINE))
	@test ! -d $(RUNDIR)/obj_$(MACHINE) || (echo removing $(RUNDIR)/obj_$(MACHINE) directory && \
	   rm -r $(RUNDIR)/obj_$(MACHINE))
	@rm -f $(RUNDIR)/$(EXEC)
	@ln -s $(EXECDIR)/$(EXEC) $(RUNDIR)
endif
	@ln -sf $(ABSBASEDIR)/bin/$(EXEC) $(BASEDIR)/testsuite
endif

#------- Create executable ----- #
#$(MAKE) -f $(MKFILE) $@
$(EXECDIR)/$(EXEC): checkpath
	$(MAKE) -f $(MKFILE) $@

#coupler: checkpath
#	$(MAKE) -f $(MKFILE) coupler

#basic: checkpath
#	$(MAKE) -f $(MKFILE) basic


#------- Create GENE library ---- #
lib:	checkpath
	$(MAKE) -f $(MKFILE) $(EXECDIR)/$(LIBN).a

#------- Create personality check program on BGP
BGP_test: checkpath
	@echo Calling make in subdir for building BGP_test
	$(MAKE) -f $(MACHMKFILE) $(EXECDIR)/BGP_test

#------- Run executable  ---- #
run:	checkpath
	$(MAKE) -f $(MKFILE) $@

#-----------------------------------------------------------------------------
clean:
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.a
	-rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.flc $(SRCDIR)/gread*.fh $(SRCDIR)/gwrite*.fh $(SRCDIR)/adios_include.h
	-rm -f $(EXECDIR)/$(EXEC) $(EXECDIR)/*.a
	-rm -f opari.rc opari.tab.c *.opari.inc opari.tab.o
	-rm -f $(SRCDIR)/*.mod.F90
	if [ -x $(EXTDIR)/$(MACHINE)/futils/src ] ; then \
	$(MAKE) -s -C $(EXTDIR)/$(MACHINE)/futils/src clean ; fi

distclean: clean
	-rm -f $(OBJDIR)/*~
	-rm -R -f $(OBJDIR)
	-rm -R -f $(PPDIR)
	-rm -f $(EXECDIR)/*.mk
	if [ -x $(EXTDIR)/$(MACHINE)/futils/src ] ; then \
	$(MAKE) -s -C $(EXTDIR)/$(MACHINE)/futils/src distclean ; fi
	$(MAKE) -s -C doc distclean
	$(MAKE) -s -C $(EXTDIR) distclean

doc:
	$(MAKE) -s -C doc

update:
ifeq ($(VERSION),)
#	@echo please set VERSION: gmake update VERSION=X.Y
else
#	@echo switching to release-$(VERSION);\
#	(svn switch https://solps-mdsplus.aug.ipp.mpg.de/repos/GENE11/tags/release-$(VERSION)/ .) || \
#	(svn switch http://solps-mdsplus.aug.ipp.mpg.de/repos/GENE11/tags/release-$(VERSION)/ .) || \
#	echo ERROR: release-$(VERSION) not found!
endif

mach_wrapper:
	@echo $(MACHINE)

$(OBJDIR)/%.o: checkpath
	$(MAKE) -f $(MKFILE) $@

forcheck preproc: checkpath
	$(MAKE) -f $(MKFILE) $@

show_objlist show_srcnames show_pp_srcnames:
	$(MAKE) -f $(MKFILE) $@

liblinkline:
	@$(MAKE) -s -f $(MKFILE) $@

###############################################################################
