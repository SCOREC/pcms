PROGRAM main
!
!   Parrallel read dataset created by pex1 with "getarr".
!   Only "ionode" effectively reads it!
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5'
  CHARACTER(len=32) :: name='/array_col'
  INTEGER :: nx=5, ny=8
  INTEGER :: ierr, fid, me, npes, comm=MPI_COMM_WORLD
  INTEGER :: start, nxp, nyp, i, j, k, n, ionode = 1
  DOUBLE PRECISION, ALLOCATABLE :: array(:,:)
  REAL, ALLOCATABLE :: garray(:,:)
!===========================================================================
!
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
  CALL dist1d(0, ny, start, nyp)
!!$  PRINT*, 'PE', me, '  start, nx, nyp', start, nx, nyp
  ALLOCATE(array(nx,nyp), garray(nx,ny))
  array = 0.0
  garray=0.0
!
  CALL openf(file, fid, mpicomm=MPI_COMM_WORLD)
  PRINT*, 'PE', me, file, 'Opened'
  CALL getarr(fid, name, array, pardim=2)
  PRINT*, 'PE', me, ' first getarr ok'
  CALL getarr(fid, name, garray, ionode=ionode)
  PRINT*, 'PE', me, ' second getarr ok'
!
  DO j=0,npes-1
     IF( me.EQ.j) THEN
        WRITE(6,'(a,i3.3)') 'PE', me
        DO i=1,nx
           WRITE(6,'((10f6.0))') (garray(i,k), k=1,ny)
        END DO
     END IF
     CALL mpi_barrier(comm, ierr)
  END DO
!
  DEALLOCATE(array)
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
