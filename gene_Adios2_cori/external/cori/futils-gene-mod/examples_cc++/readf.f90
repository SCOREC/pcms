  SUBROUTINE readf(file, name, a, nx, ny)
!
!    Simple read dataset created by pex1 with "getarr"
!
    USE futils
    IMPLICIT NONE
    CHARACTER(len=*) :: file, name
    INTEGER, INTENT(in) :: nx, ny
    INTEGER :: fid, i
    REAL(8), DIMENSION(nx,ny) :: a
!
    CALL openf(file, fid)
    a = 0.0
    CALL getarr(fid, name, a)
!
    CALL closef(fid)
  END SUBROUTINE readf
