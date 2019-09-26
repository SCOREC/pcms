PROGRAM main
!
!    Simple read dataset created by pex1 with "getarr"
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='para.h5'
  CHARACTER(len=32) :: name='/array_col'
  INTEGER, PARAMETER :: nx=10, ny=10
  INTEGER :: fid, i
  DOUBLE PRECISION, DIMENSION(nx,ny) :: array
!===========================================================================
!
  CALL openf(file, fid)
  CALL getarr(fid, name, array)
!
  DO i=1,nx
     WRITE(*,'(10f5.0)') array(i,:)
  END DO
!
  CALL closef(fid)
END PROGRAM main
