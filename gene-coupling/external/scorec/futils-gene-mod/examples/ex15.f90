PROGRAM main
!
!   Test "getatt" with optional argumenet "err"
!
  USE hdf5
  USE futils
  IMPLICIT NONE
  CHARACTER(len=128) :: file='sim.h5', label
  INTEGER :: fid
  INTEGER :: n, ierr
  CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: attnames
  CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: atttypes
  INTEGER(SIZE_T), DIMENSION(:), ALLOCATABLE :: attsizes
  INTEGER :: i, natts
  INTEGER :: ival
  REAL :: sval
  DOUBLE PRECISION :: rval
  CHARACTER(len=32) :: cval
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
  CALL openf(file, fid)
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' open'
!
  CALL getatt(fid, '/', 'prec', ival, ierr)
  IF( ierr .EQ. -2 ) THEN
     CALL getatt(fid, '/', 'prec', cval, ierr)
     PRINT*, 'prec is a string = ', TRIM(cval)
  ELSE
     PRINT*, 'prec is a integer =', ival
  END IF
  CALL closef(fid)
!
END PROGRAM main
