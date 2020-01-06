PROGRAM main
!
!   Test getdims
!
  USE futils
  IMPLICIT NONE
  DOUBLE PRECISION  :: arr1(10), arr2(10,20)
  DOUBLE COMPLEX    :: zarr3(10,20,30)
  CHARACTER(len=32) :: file='getdims.h5'
  INTEGER :: fid
  INTEGER :: rank, dims(7)   !  Max rank is 7 in Fortran
!
  CALL creatf(file, fid)
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' created'
!
  arr1=1.0
  arr2=2.0
  zarr3=(3.0, -3.0)
!
  CALL putarr(fid, '/arr1', arr1)
  CALL putarr(fid, '/arr2', arr2)
  CALL putarr(fid, '/zarr3', zarr3)
!
  CALL getdims(fid, '/arr1', rank, dims)
  WRITE(*,'(a,7i6)') 'Dims of arr1 are ', dims(1:rank)
  CALL getdims(fid, '/arr2', rank, dims)
  WRITE(*,'(a,7i6)') 'Dims of arr2 are ', dims(1:rank)
  CALL getdims(fid, '/zarr3', rank, dims)
  WRITE(*,'(a,7i6)') 'Dims of zarr3 are', dims(1:rank)
!
  CALL closef(fid)
END PROGRAM main
