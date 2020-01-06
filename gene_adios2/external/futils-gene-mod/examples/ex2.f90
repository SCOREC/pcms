PROGRAM main
!
!  Append datasets to an existing file
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='sim.h5'
  INTEGER :: fid
  INTEGER :: ierr, n, nprev, istep, ibuf, nrun=200
  INTEGER, PARAMETER :: bufsize=20
  LOGICAL :: ok
  DOUBLE PRECISION :: buf(bufsize, 0:2) ! To store hist. arrays for scalars
  DOUBLE PRECISION :: time, ekin, epot
  !
!===========================================================================
!                   1. Prologue
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL openf(file, fid)
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' open'
  CALL getatt(fid, '/', 'Success', ok)
  PRINT*, 'Success', ok
!===========================================================================
!                   2. Time loop
!
  CALL getsize(fid, '/scalars/time', nprev)
  WRITE(*,'(a,i6)') 'Number of previous steps =', nprev
  ibuf=0
  DO istep=1,nrun
     time = nprev+istep
     ekin = COS(0.2*time)*EXP(0.01*time)
     epot = SIN(0.2*time)*(1.0-EXP(0.01*time))
!
     ibuf = ibuf+1
     buf(ibuf,0) = time
     buf(ibuf,1) = ekin
     buf(ibuf,2) = epot
     IF( ibuf.EQ.bufsize .OR. istep.EQ.nrun) THEN ! Dump the buffers to file
        CALL append(fid, "/scalars/time", buf(1:ibuf,0))
        CALL append(fid, "/scalars/ekin", buf(1:ibuf,1))
        CALL append(fid, "/scalars/epot", buf(1:ibuf,2))
        ibuf = 0
     END IF
  END DO
  CALL getsize(fid, "/scalars/time", n)
  WRITE(*,'(a,i6)') 'Number of steps =', n
  CALL attach(fid, '/', 'Success', .NOT.ok)
!
  CALL putarr(fid,'/buf1', buf)
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
!
END PROGRAM main
