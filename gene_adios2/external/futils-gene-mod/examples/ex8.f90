PROGRAM main
!
!   Create a new dataset with many attributes
!   Then scan it to get these attributes
!
  USE hdf5
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='sim.h5', label, str='F'
  CHARACTER(len=16), DIMENSION(:), ALLOCATABLE :: attnames
  CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: atttypes
  INTEGER(SIZE_T), DIMENSION(:), ALLOCATABLE :: attsizes
  INTEGER :: fid
  INTEGER, DIMENSION(1000) :: iarr
  INTEGER :: nsd=33, ndim=3, i, natts
  REAL :: pi
  INTEGER :: ival
  REAL :: sval
  DOUBLE PRECISION :: rval
  CHARACTER(len=32) :: cval
!
!   Create a new file
!
  CALL creatf(file, fid, 'A file with many attributes in a dataset',&
       &       real_prec='s')
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' created'
!
!   Put array into a new dataset
!
  label = "ADIR"
  DO i=1,nsd
     iarr(i) = i
  END DO
  CALL putarr(fid, label, iarr, "Memcom ADIR")
!
!   Attach some attributes to /ADIR
!
  pi = 4.0*ATAN(1.0d0)
  CALL attach(fid, label, "NSD", nsd)
  CALL attach(fid, label, "NDIM", ndim)
  CALL attach(fid, label, "PI", pi)
  CALL attach(fid, label, "PREC", str)
!
!   Close and reopen the file
!
  CALL closef(fid)
  CALL openf(file, fid)
!
!   Get all the attributes of ADIR
!
  natts = numatts(fid, label)
  WRITE(*,'(/a,i3)') 'Number of attributes in '//TRIM(label)//' =', natts
!
  ALLOCATE(attnames(natts), atttypes(natts), attsizes(natts))
!
  CALL allatts(fid, label, attnames, atttypes, attsizes)
  DO i=1,natts
     SELECT CASE (atttypes(i))
     CASE ('I')
        CALL getatt(fid, label, attnames(i), ival)
        WRITE(*,'(a8, i8)') TRIM(attnames(i))//' = ', ival
     CASE ('C')
        CALL getatt(fid, label, attnames(i), cval)
        WRITE(*,'(a8, a)') TRIM(attnames(i))//' = ', &
             &   cval(1:attsizes(i))
     CASE ('S')
        CALL getatt(fid, label, attnames(i), sval)
        WRITE(*,'(a8, 1pe12.3)') TRIM(attnames(i))//' = ', sval
     CASE ('R')
        CALL getatt(fid, label, attnames(i), rval)
        WRITE(*,'(a8, 1pe12.3)') TRIM(attnames(i))//' = ', rval
     END SELECT
  END DO
!
!   Epilogue
!
  DEALLOCATE(attnames)
  DEALLOCATE(atttypes)
  DEALLOCATE(attsizes)
  CALL closef(fid)
END PROGRAM main
