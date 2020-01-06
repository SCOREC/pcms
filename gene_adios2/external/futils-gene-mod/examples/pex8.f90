PROGRAM main
!
!  Parallel append 0d vars (history 0d arrays)
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=256) :: file='sim.h5'
  INTEGER :: ierr, fid, me, npes
  INTEGER :: rank, dims(7), istep, ibuf, n, nrun=10000
  INTEGER, PARAMETER :: bufsize=1000
  DOUBLE PRECISION :: buf(bufsize, 0:2) ! To store hist. arrays for scalars
  DOUBLE PRECISION :: time, ekin, epot, t0, t1, mb, mbs
!===========================================================================
!                             1.  Prologue
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Create file collectively
!
  CALL creatf(file, fid, 'A simple simulation', &
       &      real_prec='s', &
       &      mpiposix=.FALSE., &
       &      mpicomm=MPI_COMM_WORLD)
!
!   Create datasets collectively
!
  CALL creatg(fid, "/scalars", "Time evolution for Scalars")
  WRITE(*,'(a)') 'group scalars created'
!
  rank = 0
  CALL creatd(fid, rank, dims, "/scalars/time", "Normalized Time")
  CALL creatd(fid, rank, dims, "/scalars/ekin", "Kinetic Energy")
  CALL creatd(fid, rank, dims, "/scalars/epot", "Potential Energy")
  WRITE(*,'(a)') 'extendible datasets in /scalars created'
!
!===========================================================================
!                   2. Time loop
!
  nrun = MAX(nrun, bufsize)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  t0 = mpi_wtime()
  ibuf=0
  DO istep=1,nrun
     time = istep
     ekin = COS(0.2*time)*EXP(0.001*time)
     epot = SIN(0.2*time)*(1.0-EXP(0.001*time))
!
     IF( bufsize .LE. 1 ) THEN
        IF( me.EQ.0 ) PRINT*, 'Dump the buffers at istep = ', istep
        CALL append(fid, "/scalars/time", time)
        CALL append(fid, "/scalars/ekin", ekin)
        CALL append(fid, "/scalars/epot", epot)
!!$        CALL append(fid, "/scalars/time", time, ionode=0)
!!$        CALL append(fid, "/scalars/ekin", ekin, ionode=0)
!!$        CALL append(fid, "/scalars/epot", epot, ionode=0)
     ELSE
        ibuf = ibuf+1
        buf(ibuf,0) = time
        buf(ibuf,1) = ekin
        buf(ibuf,2) = epot
        IF( ibuf.EQ.bufsize .OR. istep.EQ.nrun) THEN ! Dump the buffers to file
           IF( me.EQ.0 ) PRINT*, 'Dump the buffers at istep = ', istep
           CALL append(fid, "/scalars/time", buf(1:ibuf,0))
           CALL append(fid, "/scalars/ekin", buf(1:ibuf,1))
           CALL append(fid, "/scalars/epot", buf(1:ibuf,2))
!!$           CALL append(fid, "/scalars/time", buf(1:ibuf,0), ionode=0)
!!$           CALL append(fid, "/scalars/ekin", buf(1:ibuf,1), ionode=0)
!!$           CALL append(fid, "/scalars/epot", buf(1:ibuf,2), ionode=0)
           ibuf = 0
        END IF
     END IF
  END DO
  CALL getsize(fid, "/scalars/time", n)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  t1 = mpi_wtime() - t0
  mb = 4.*3.*n/1024.0/1024.0
  mbs = mb/t1
  IF( me .EQ. 0 ) THEN
     WRITE(*,'(a, i6)')   'Number of steps  =', n
     WRITE(*,'(a, i6)')   'Buffer size      =', bufsize
     WRITE(*,'(a, f8.3)') 'Elapsed time (s) =', t1
     WRITE(*,'(a, f8.3)') 'Size of file (MB)=', mb
     WRITE(*,'(a, f8.3)') 'BW (MB/s)        =', mbs
  END IF
!
!===========================================================================
!                             9.  Epilogue
!
  CALL closef(fid)
  CALL mpi_finalize(ierr)
END PROGRAM main
