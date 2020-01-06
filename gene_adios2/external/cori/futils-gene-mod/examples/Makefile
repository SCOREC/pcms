#adapt the following to machine
F90 ?= ftn
CC ?= cc
HDF5 ?= $(HDF5_PARALLEL_DIR)
MPIRUN = srun -n 4 --export='HDF5_USE_FILE_LOCKING=FALSE'

OPT ?= -O3
#OPT = -g -traceback   #  debug flags
F90FLAGS ?= $(OPT) -I../src/ -I${HDF5}/include
CFLAGS ?= -O2
LDFLAGS ?= $(OPT) -L../src -L${HDF5}/lib

LIBS =  -lfutils -lhdf5_fortran -lhdf5 -lz

SERIAL = ex1 ex2 ex3 ex4 ex7 ex8 ex9 ex11 ex12 ex13 ex14 ex15
PARA   = pex1 pex3 pex4 pex5 pex6 pex7 pex8 pex9 ex10
PARA_ALL = $(PARA) pex5r pex10 pex11 pex12 pex13 pex14 pex15

.SUFFIXES:
.SUFFIXES: .o .c .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

examples: $(SERIAL) ex6 ex5 $(PARA_ALL)

lib:
	make -C ../src
	touch lib

ex1:	ex1.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex2:	ex2.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex3:	ex3.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex4:	ex4.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex5:	ex5.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex6:	ex6.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex7:	ex7.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex8:	ex8.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex9:	ex9.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex10: ex10.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex11: ex11.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex12: ex12.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex13: ex13.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex14: ex14.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

ex15: ex15.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex1:	pex1.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex3:	pex3.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex4:	pex4.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex5:	pex5.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex5r:	pex5r.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex6:	pex6.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex7:	pex7.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex8:	pex8.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex9:	pex9.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex10:	pex10.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex11:	pex11.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex12:	pex12.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex14:	pex14.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex13:	pex13.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex15:	pex15.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex16:	pex16.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

pex17:	pex17.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

getfile: getfile.o
	$(F90) $(LDFLAGS) -o $@ $<  $(LIBS)

testmem: testmem.o
	$(F90) $(LDFLAGS) -o $@ $< $(LIBS)

test_s:	$(SERIAL) ex5
	@for x in $(SERIAL); do \
	echo === Running $$x ===;\
	./$$x ;\
        echo ;\
	done
	@echo === Running ex5 ===;\
	./ex5 < ex5.in ;\
        echo ;\

test_p:	$(PARA_ALL)
	@for x in $(PARA); do \
	echo === Running $$x ===;\
	$(MPIRUN) ./$$x ;\
        echo ;\
	done ;\
	echo === Running pex10 ===;\
	$(MPIRUN) -n 12 ./pex10 ;\
        echo ;\
	echo === Running pex11 ===;\
	$(MPIRUN) -n 8 ./pex11 ;\
        echo ;\
	echo === Running pex12 ===;\
	$(MPIRUN) -n 4 ./pex12 2 2;\
        echo ;\
	echo === Running pex13 ===;\
	$(MPIRUN) -n 6 ./pex13 3 2 ;\
        echo ;\
	echo === Running pex14 ===;\
	$(MPIRUN) -n 4 ./pex14 2 2;\
        echo ;\
	echo === Running pex15 ===;\
	$(MPIRUN) -n 6 ./pex15 3 2 ;\
        echo ;\

ex1.o:	lib
ex2.o:	lib
ex3.o:	lib
ex4.o:	lib
ex5.o:	lib
ex6.o:	lib
ex7.o:	lib
ex8.o:	lib
ex9.o:	lib
ex10.o: lib
ex11.o: lib
ex12.o: lib
ex13.o: lib
ex14.o: lib
ex15.o: lib
pex1.o:	lib
pex3.o:	lib
pex4.o:	lib
pex5.o:	lib
pex5r.o: lib
pex6.o:	lib
pex7.o:	lib
pex8.o:	lib
pex9.o:	lib
pex10.o:	lib
pex11.o:	lib
pex12.o:	lib
pex14.o:	lib
pex13.o:	lib
pex15.o:	lib
pex16.o:	lib
pex17.o:	lib
getfile.o:	lib

clean:
#	make -C ../src clean
	rm -f *.o *~ a.out fort.* *.h5*

distclean: clean
#	make -C ../src distclean
	rm -f lib
	rm -f $(SERIAL) ex5 ex6 $(PARA_ALL) *.mod TAGS
