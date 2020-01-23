FUTILS=../src
HDF5=/usr/local/hdf5

F90 = mpif90
OPT = -g -traceback  -CB
F90FLAGS = $(OPT) -I$(FUTILS) -I${HDF5}/lib
CC = icc
CFLAGS = -g
LDFLAGS = $(OPT) -L$(FUTILS) -L${HDF5}/lib
LIBS =  -lfutils -lhdf5_fortran -lhdf5 -lz
CXX = icpc
CXXFLAGS=$(CFLAGS)

.SUFFIXES:
.SUFFIXES: .o .c .f90 .cpp

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all: testr testr++

testr++: main++.o readf.o
	$(F90) -o $@ -nofor_main -cxxlib $(LDFLAGS) main++.o readf.o $(LIBS)

testr:	main.o readf.o
	$(F90) -o $@ -nofor_main $(LDFLAGS) main.o readf.o $(LIBS)

clean:
	rm -f *.o *~ a.out
distclean: clean
	rm -f testr testr++ *.mod *.h5 TAGS
