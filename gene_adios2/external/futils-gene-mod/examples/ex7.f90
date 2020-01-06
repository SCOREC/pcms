PROGRAM main
!
!    Write "doubles" and read into "simple precision" arrays
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='file.h5'
  CHARACTER(len=256) :: str
  INTEGER :: fid, i, ierr
  INTEGER, PARAMETER :: n=20
  REAL*8 :: darr(n)
  REAL*8 :: sarr(n)
!
  CALL creatf(file, fid, real_prec='d')
  CALL RANDOM_NUMBER(darr)
  WRITE(*,'(a/(10(1pe12.3)))') 'DARR =', darr
  CALL putarr(fid, 'doubles', darr)
  CALL attach(fid, 'doubles', 'Comment', 'Written by ex7!')
!!$  CALL closef(fid)
  CALL closeall(ierr)
  PRINT*, 'closeall ierr :', ierr
!
!============
  CALL openf(file, fid)
  CALL getarr(fid, 'doubles', sarr)
  CALL getatt(fid, 'doubles', 'Comment', str)
  WRITE(*,'(a,a)') 'Comment = ', TRIM(str)
  WRITE(*,'(a/(10(1pe12.3)))') 'SARR =', sarr
  WRITE(*,'(a/(10(1pe12.3)))') 'DIFF =', darr-sarr
  CALL closef(fid)
  CALL closeall(ierr)
  PRINT*, 'closeall ierr :', ierr
END PROGRAM main
