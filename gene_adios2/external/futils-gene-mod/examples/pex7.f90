PROGRAM main
!
!   P-IO by selected nodes
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5'
  INTEGER :: nx=5
  INTEGER :: ierr, fid, me, npes, comm=MPI_COMM_WORLD
  INTEGER :: start, nxp, nyp, i, j, n
  DOUBLE PRECISION, ALLOCATABLE :: array(:), parray(:)
  REAL, ALLOCATABLE :: sarray(:), psarray(:)
  INTEGER, ALLOCATABLE :: iarray(:), piarray(:)
  INTEGER :: ionode
!===========================================================================
!                             1.  Prologue
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Create file collectively
!
  CALL creatf(file, fid, &
       &      desc="A parallel file", &
       &      real_prec='s', &
       &      mpicomm=MPI_COMM_WORLD)
!
!   *All* processors allocate array
!
  ALLOCATE(array(nx), iarray(nx), sarray(nx))
!
!   Only array from ionode will be written
!
  ionode = npes-1
  DO i=1,nx
     array(i) = 10*me + i
  END DO
  iarray = -array
  sarray = 1.0/array
  CALL putarr(fid, "/array", array, ionode=ionode)
  CALL putarr(fid, "/iarray", iarray, ionode=ionode)
  CALL putarr(fid, "/sarray", sarray, ionode=ionode)
!
!   Partitionned array
!
  ALLOCATE(parray(nx), piarray(nx), psarray(nx))
  DO i=1,nx
     parray(i) = me*nx + i
  END DO
  piarray = -parray
  psarray = 1.0/array
  CALL putarr(fid, "/parray", parray, pardim=1)
  CALL putarr(fid, "/piarray", piarray, pardim=1)
  CALL putarr(fid, "/psarray", psarray, pardim=1)
!
!===========================================================================
!                             9.  Epilogue
!
  DEALLOCATE(array)
  DEALLOCATE(parray)
  CALL closef(fid)
  CALL mpi_finalize(ierr)
END PROGRAM main
