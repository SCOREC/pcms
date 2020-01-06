PROGRAM main
!
!   Multi-dim (+UNLIMITED time dim) profiles
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='pprof.h5'
  INTEGER, PARAMETER :: nx=256, ny=512
  INTEGER :: me, npes, offset, nyp
  INTEGER :: fid, rank, dims(7)
  INTEGER :: i, j, istep, ierr, n, nrun=50
  DOUBLE PRECISION, ALLOCATABLE :: xg(:), yg(:), phi2(:,:)
  DOUBLE PRECISION :: dx, dy, pi, time, x0, s
!===========================================================================
!                   1. Prologue
!
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!  1D partition in Y-direction
  CALL dist1d(0, ny, offset, nyp)
  WRITE(*,'(a,i3.3,a,2i6)') 'PE', me, ': offset, nyp', offset, nyp
  ALLOCATE(xg(nx), yg(nyp), phi2(nx,nyp))
!
!  Init HDF5
  CALL creatf(file, fid, 'A simple simulation', mpicomm=MPI_COMM_WORLD)
  CALL creatg(fid, "/profile_2d", "1D profiles")
!
  rank = 0
  CALL creatd(fid, rank, dims, "/profile_2d/time", "Normalized Time")
!
  rank = 2
  dims(1) = nx; dims(2) = nyp
  CALL creatd(fid, rank, dims, "/profile_2d/phi", "Potential", pardim=2)
  IF( me .EQ. 0 ) WRITE(*,'(a)') 'extendible datasets in /profile_1d created'
!===========================================================================
!                   2. Initialization
!
  pi = 4*ATAN(1.0d0)
  dx = 2.*pi/REAL(nx-1)
  dy = pi/real(ny-1)
  DO i=1,nx
     xg(i) = (i-1)*dx
  END DO
  DO j=1,nyp
     yg(j) = (j+offset-1)*dy
  END DO
  CALL putarr(fid, "/profile_2d/xg", xg, "x mesh")
  IF( me .EQ. 0 ) WRITE(*,'(a)') 'x-mesh'
  CALL putarr(fid, "/profile_2d/yg", yg, "y mesh", pardim=1)
  IF( me .EQ. 0 ) WRITE(*,'(a)') 'y-mesh'
!===========================================================================
!                   3. Time loop
!
  x0 = 0.5
  DO istep=1,nrun
!
     time = istep-1
     s = 0.1+0.1*time
     DO i=1,nx
        DO j=1,nyp
           phi2(i,j) = SIN(xg(i)) * COS(yg(j)) * COS(0.04*pi*time)
        END DO
     END DO
!
     CALL append(fid, "/profile_2d/time", time, ionode=0) ! Only rank=0 writes it
!!$     CALL append(fid, "/profile_2d/time", time)
     CALL append(fid, "/profile_2d/phi", phi2, pardim=2)
     IF( me .EQ. 0 ) WRITE(*,'(a,i5,a)') 'Step', istep, ' done'
!
  END DO
!===========================================================================
!                   9. Epilogue
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
