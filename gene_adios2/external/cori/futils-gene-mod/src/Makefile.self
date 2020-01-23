#Makefile for compilation of futils without GENE dependency
#change the following flags and set a prefix if needed. Tests have a
#separate Makefile requiring the MPI wrapper.

#PREFIX=/global/homes/m/merlo/futils-test
HDF5=$(HDF5_PARALLEL_DIR)
F90 = ftn
CC = cc

OPT = -O3
#OPT = -g -traceback
F90FLAGS = $(OPT) -I${HDF5}/lib -I${HDF5}/include
CFLAGS = -O2
LDFLAGS = $(OPT) -L. -L${HDF5}/lib
LIBS =  -lfutils -lhdf5_fortran -lhdf5 -lz

#F2003 HDF5 interface
#IBM compilers
ifeq ($(findstring mpcc,$(CC)),mpcc)
#add -qlanglvl=stdc99 on IBM compilers if not set in CFLAGS
CFLAGS += -qlanglvl=stdc99
endif
ifeq ($(findstring xlf,$(F90)),xlf)
F90FLAGS += -qsuffix=cpp=f90 -qextname=fu_fsize:fu_ftos:fu_stof:fu_stostdout
# the following variable is set by loading the
# HDF5 module.
ifeq ($(HDF5_HAS_FORTRAN2003_INTERFACES),yes)
F90FLAGS += -qlanglvl=2003std -WF,-DHDF5_FORTRAN2003 -qflag=i:w
endif
endif

#mpif90
ifeq ($(findstring mpif90,$(F90)),mpif90)
F90FLAGS += -fpp
ifeq ($(HDF5_HAS_FORTRAN2003_INTERFACES),yes)
F90FLAGS += -stand f2003 -DHDF5_FORTRAN2003
endif
endif

# Vampirtrace compiler wrappers
ifeq ($(findstring vtf90,$(F90)),vtf90)
F90FLAGS += -fpp
ifeq ($(HDF5_HAS_FORTRAN2003_INTERFACES),yes)
F90FLAGS += -stand f2003 -DHDF5_FORTRAN2003
endif
endif

ifeq ($(findstring mpiifort,$(F90)),mpiifort)
F90FLAGS += -fpp
ifeq ($(HDF5_HAS_FORTRAN2003_INTERFACES),yes)
F90FLAGS += -stand f2003 -DHDF5_FORTRAN2003
endif
endif

ifeq ($(findstring ftn,$(F90)),ftn)
F90FLAGS += -fpp
ifeq ($(HDF5_HAS_FORTRAN2003_INTERFACES),yes)
F90FLAGS += -stand f2003 -DHDF5_FORTRAN2003
endif
endif

.SUFFIXES:
.SUFFIXES: .o .c .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

lib: libfutils.a getfile

getfile: getfile.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

libfutils.a: futils.o cutils.o buffer.o
	ar r $@ $?
	ranlib $@

futils.o: append.tpl zappend.tpl \
          putarr.tpl cputarr.tpl putarrnd.tpl cputarrnd.tpl \
          getarr.tpl cgetarr.tpl getarrnd.tpl cgetarrnd.tpl

buffer.o: futils.o

install:
	make distclean
	make -f Makefile.self libfutils.a
	mkdir -p $(PREFIX)/{lib,include}
	cp -p libfutils.a $(PREFIX)/lib
	cp -p *.mod $(PREFIX)/include

test:
	make  "F90_FFLAGS=$(F90FLAGS)" "LD_FLAGS=$(LDFLAGS)" \
	-C ../examples test_s test_p

getfile.o:	libfutils.a

clean:
	-rm -f *.o *~ a.out fort.* *.h5* *.mod *.a getfile

distclean: clean
	rm -f getfile *.a *.mod TAGS
	make distclean -C ../examples