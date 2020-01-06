PROGRAM main
!
!   Store scalar (0-dim) history arrays with unlimited dim for
!   time
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='sim.h5'
  INTEGER :: fid, n, istep, ibuf, ierr, nrun=120
  INTEGER :: rank, dims(7)
  INTEGER, PARAMETER :: bufsize=20
  DOUBLE PRECISION :: buf(bufsize, 0:2) ! To store hist. arrays for scalars
  DOUBLE PRECISION :: time, ekin, epot
  CHARACTER(len=16) :: libver
  INTEGER :: l
!===========================================================================
!                   1. Prologue
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
!!$  CALL creatf(file, fid, 'A simple simulation')
  CALL creatf(file, fid, 'A simple simulation', real_prec='d')
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' created'
!
  CALL geth5ver(libver, l)
  WRITE(*,'(a,a)') 'HDF5 library version used: ', libver(1:l)
!
  CALL creatg(fid, "/scalars", "Time evolution for Scalars")
  WRITE(*,'(a)') 'group scalars created'
!
  rank = 0
  CALL creatd(fid, rank, dims, "/scalars/time", "Normalized Time")
  CALL creatd(fid, rank, dims, "/scalars/ekin", "Kinetic Energy")
  CALL creatd(fid, rank, dims, "/scalars/epot", "Potential Energy")
  WRITE(*,'(a)') 'extendible datasets in /scalars created'
!===========================================================================
!                   2. Time loop
!
  ibuf=0
  DO istep=1,nrun
     time = istep
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
  WRITE(*,'(a, i6)') 'Number of steps =', n
  CALL attach(fid, '/', 'Success', .TRUE.)
!
  PRINT*, "Is '/' a group?", isgroup(fid, '/')
  PRINT*, "Is '/scalars' a group?", isgroup(fid, '/scalars')
  PRINT*, "Is '/scalars/time' a group?", isgroup(fid, '/scalars/time')
!
  CALL putarr(fid,'/buf0', buf)
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
END PROGRAM main
