PROGRAM main
!
!   Put files in HDF5 file
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=256) :: file='file.h5', path, name, base, dir, line
  INTEGER :: n, fid, ierr
!===========================================================================
!                   1. Prologue
!
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL creatf(file, fid, 'A simple simulation')
  CALL creatg(fid, "/inputs", "Inputs specific to the run")
!===========================================================================
!                    2. Put some files
!
  path = 'README.txt'
  CALL split(path, dir, base)
  name = '/inputs/' // base(1:LEN_TRIM(path))
  CALL putfile(fid, name, path, compress=.TRUE.)
  PRINT*, 'dataset ', TRIM(name), ' written'
!
  path = 'Makefile'
  CALL split(path, dir, base)
  name = '/inputs/' // base(1:LEN_TRIM(path))
  CALL putfile(fid, name, path, compress=.TRUE.)
  PRINT*, 'dataset ', TRIM(name), ' written'
!===========================================================================
!                    3. Save standard input
!
  DO
     READ(*,'(a)', END=110) line    ! Capture stdin
     WRITE(90, '(a)') TRIM(line)
  END DO
110 CONTINUE
  REWIND(90)
!
!   Now read stdin using "READ(90" instead of "READ(*,"
!
  INQUIRE(unit=90, name=path)
  CALL putfile(fid, '/inputs/STDIN', TRIM(path)) ! Save stdin in /inputs/STDIN
  PRINT*, 'dataset /inputs/STDIN written'
  CLOSE(90,status='delete')                      ! Clean up
!===========================================================================
!                   9. Epilogue
!
9900 CONTINUE
  CALL closef(fid)
END PROGRAM main
