PROGRAM main
!
!   Create 3D (+time as 3rd extendible dim) datasets with "append"
!   Read slices (in time) of 3d array
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32)  :: file='prof_slice.h5'
  INTEGER, PARAMETER :: nx=256, ny=256, nrun=20, slice=2
  INTEGER            :: fid, rank, dims(7), offsets(3)
  INTEGER            :: i, j, istep, is, nerrs
  DOUBLE PRECISION   :: xg(nx), yg(ny), phi3(nx,ny,nrun), phi_read(nx,ny,slice)
  DOUBLE PRECISION   :: dx, dy, pi, time, x0, s
!===========================================================================
!
!  Create file using DB reals
  CALL creatf(file, fid, real_prec='d')
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' open'
!
!  Create group
  CALL creatg(fid, "/profile_2d", "2D profiles")
  WRITE(*,'(a)') "group profile_12d created"
!
!  Create extendible dataset
  CALL creatd(fid, 2, (/nx, ny/), "/profile_2d/phi")
!
!  Mesh
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
!  Time loop
  x0 = 0.5
  DO istep=1,nrun
     time = istep-1
     DO i=1,nx
        DO j=1,ny
           phi3(i,j,istep) = SIN(xg(i)) * COS(yg(j)) * COS(0.04*pi*time)
        END DO
     END DO
     CALL append(fid, "/profile_2d/phi", phi3(:,:,istep))
  END DO
  CALL closef(fid)
!
!   Reopen for read
  CALL openf(file, fid)
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' reopen'
!
!   Read 2d profiles and check
  offsets = 0
  DO istep=1,nrun,slice
     offsets(3) = istep-1
     CALL getarr(fid, "/profile_2d/phi", phi_read, offsets=offsets)
     nerrs=0
     DO i=1,nx
        DO j=1,ny
           DO is=1,slice
              IF(phi_read(i,j,is) .NE. phi3(i,j,is+istep-1) ) nerrs=nerrs+1
           END DO
        END DO
     END DO
     WRITE(*,'(a,i3,i6)') 'istep, nerrs', istep, nerrs
  END DO
!
!  Clean up
  CALL closef(fid)
END PROGRAM main
