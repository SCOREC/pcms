Serial examples:
===============
ex1:	Extendible 0-dim arrays
ex2:	Append to exsiting datasets created in ex1
ex3:	1d/2d time-dependent profiles using "append"
ex4:	Create particle arrays using "putarr"
ex5:    Example with "putfile"
ex6:    Simple read dataset created by pex1 with "getarr"
ex7:    Write "doubles" and read into "simple precision" arrays
ex8:    Get all its attributes in dataset/group
ex9:    Writes in slices to reduce storage required
ex10:   Test module HASTABLE
ex11:   Write and read complex arrays
ex12:   Test copy_file/move_file (which calls copy_file) from cutils.c

Parallel examples:
=================
pex1:	2d array write with putarr
pex3:	1d/2d time-dependent profiles using "append" (parallel version)
pex4:	Create particle arrays using "putarr" (paralle version)
pex5:   Parallel read particle arrays (created in pex4) using "getarr"
pex5r:	Serial version of pex5
pex6:	Test of optionall "ionode" in parallel "getarr"
pex7:   Test of optionall "ionode" in parallel "putarr"
pex8:   Parallel append 0d vars (history 0d arrays)
pex9:   Parallel write a 2d COMPLEX array (from pex1.f90)
pex10:  Parallel write/read of 2d and 3d arrays partionned on
        2d processor grid
pex11:  Parallel write/read a 3d array partionned on 2d processor grid
        A(n1/P1, n2/P2, n3), with GHOST CELLS on the partitionned dimensions.
pex12:  Create and write array with ghost cells, dynamic version
pex13:  Read array from file created by pex12
pex14:  Write 4d array with ghost cells with 2 dimensions distributed on
        2d processor grid
pex15:  Read 4d array from file created by pex14
