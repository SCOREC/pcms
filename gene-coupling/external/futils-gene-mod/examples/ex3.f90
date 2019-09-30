PROGRAM main
!
!   Multi-dim (+UNLIMITED time dim) profiles
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='prof.h5'
  INTEGER, PARAMETER :: nx=256, ny=256
  INTEGER :: fid, rank, dims(7)
  INTEGER :: i, j, istep, ierr, n, nrun=20
  DOUBLE PRECISION, DIMENSION(nx) :: xgrid, phi
  DOUBLE PRECISION :: xg(nx), yg(ny), phi2(nx,ny)
  DOUBLE PRECISION :: dx, dy, pi, time, x0, s
!===========================================================================
!                   1. Prologue
!
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL creatf(file, fid, 'A simple simulation')
!!$  CALL openf(file, fid)
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' open'
!
  CALL creatg(fid, "/profile_1d", "1D profiles")
  CALL creatg(fid, "/profile_2d", "1D profiles")
  WRITE(*,'(a)') 'group profile_1d created'
!
  rank = 0
  CALL creatd(fid, rank, dims, "/profile_1d/time", "Normalized Time")
  CALL creatd(fid, rank, dims, "/profile_2d/time", "Normalized Time")
  rank = 1
  dims(1:rank) = SHAPE(phi)
  CALL creatd(fid, rank, dims, "/profile_1d/phi", "Potential")
!!$  CALL creatd(fid, rank, dims, "/profile_1d/phi", "Potential", compress=.TRUE.)
  rank = 2
  dims(1:rank) = SHAPE(phi2)
  CALL creatd(fid, rank, dims, "/profile_2d/phi", "2d Potential")
!!$  CALL creatd(fid, rank, dims, "/profile_2d/phi", "2d Potential", compress=.TRUE.)
!
  WRITE(*,'(a)') 'extendible datasets in /profile_1d created'
!===========================================================================
!                   2. Initialization
!
! 1D
  dx = 1.0d0/REAL(nx-1)
  DO i=1,nx
     xgrid(i) = (i-1)*dx
  END DO
  CALL putarr(fid, "/profile_1d/xgrid", xgrid, "X mesh")
!
! 2D
  pi = 4*ATAN(1.0d0)
  dx = 2.*pi/REAL(nx-1)
  dy = pi/real(ny-1)
  DO i=1,nx
     xg(i) = (i-1)*dx
  END DO
  DO j=1,ny
     yg(j) = (j-1)*dy
  END DO
  CALL putarr(fid, "/profile_2d/xg", xg, "x mesh")
  CALL putarr(fid, "/profile_2d/yg", yg, "y mesh")
!
  WRITE(*,'(a)') 'fixed sized datasets in /profile_1d created'
!===========================================================================
!                   3. Time loop
!
  x0 = 0.5
  DO istep=1,nrun
!
     time = istep-1
     s = 0.1+0.1*time
     DO i=1,nx
        phi(i) = EXP( -((xgrid(i)-x0)/s)**2)
     END DO
!
     DO i=1,nx
        DO j=1,ny
           phi2(i,j) = SIN(xg(i)) * COS(yg(j)) * COS(0.04*pi*time)
        END DO
     END DO
!
     CALL append(fid, "/profile_1d/time", time)
     CALL append(fid, "/profile_1d/phi", phi)
     CALL append(fid, "/profile_2d/time", time)
     CALL append(fid, "/profile_2d/phi", phi2)
!
  END DO
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
END PROGRAM main
