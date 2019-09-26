PROGRAM main
!
!   Create an HDF5 file and test copy_file/move_file
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file0='file.h5'
  CHARACTER(len=32) :: file1='file.h5_cp'
  CHARACTER(len=32) :: file2='file.h5_mv'
  INTEGER :: fid
  LOGICAL :: ex
  REAL, DIMENSION(1000,1000) :: array, a
!
!   Create the hdf5 file
  CALL RANDOM_NUMBER(array)
  CALL creatf(file0,fid)
  CALL putarr(fid,'/RANDOM',array)
  CALL closef(fid)
!
!   move/copy the file
  CALL copy_file(file0, LEN_TRIM(file0), file1, LEN_TRIM(file1))
  CALL move_file(file0, LEN_TRIM(file0), file2, LEN_TRIM(file2))
  INQUIRE(file=TRIM(file0), exist=ex)
  IF( ex ) THEN
     PRINT*, TRIM(file0), " is still here!"
  ELSE
     PRINT*, TRIM(file0), " is gone!"
  END IF
!
!   Check contents of file1
  CALL openf(file1,fid)
  CALL getarr(fid,'/RANDOM',a)
  CALL closef(fid)
  PRINT*, 'Min/max of diff', MINVAL(array-a), MAXVAL(array-a)
END PROGRAM main
