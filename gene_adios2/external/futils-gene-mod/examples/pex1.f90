PROGRAM main
!
!   Parallel write a 2d array
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5'
  INTEGER :: nx=5, ny=8
  INTEGER :: ierr, fid, me, npes, comm=MPI_COMM_WORLD
  INTEGER :: start, nxp, nyp, i, j, n
  DOUBLE PRECISION, ALLOCATABLE :: array(:,:)
!===========================================================================
!                             1.  Prologue
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Create file collectively
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL creatf(file, fid, &
       &      desc="A parallel file", &
       &      real_prec='s', &
       &      mpiposix=.FALSE., &
       &      mpicomm=MPI_COMM_WORLD)
!
  CALL putfile(fid, '/README.txt', &
       &   'README.txt', ionode=0)
  CALL putfile(fid, '/Makefile', &
       &   'Makefile', ionode=1)
!===========================================================================
!                             2.  Parallel write file
!
!   Define local array partitionned by columns
  CALL dist1d(0, ny, start, nyp)
  ALLOCATE(array(nx,nyp))
  DO i=1,nx
     DO j=1,nyp
        array(i,j) = 10*i + (start+j)
     END DO
  END DO
  CALL putarr(fid, "/array_col", array, pardim=2)
  PRINT*, 'dataset /array_col created'
  DEALLOCATE(array)
!===========================================================================
!                             9.  Epilogue
!
  CALL closef(fid)
  CALL mpi_finalize(ierr)
END PROGRAM main

SUBROUTINE dist1d(s0, ntot, s, nloc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(in) :: s0, ntot
  INTEGER, INTENT(out) :: s, nloc
  INTEGER :: me, npes, ierr, naver, rem
!
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
  naver = ntot/npes
  rem = MODULO(ntot,npes)
  s = s0 + MIN(rem,me) + me*naver
  nloc = naver
  IF( me.LT.rem ) nloc = nloc+1
!
END SUBROUTINE dist1d
