PROGRAM main
!
!   Particle array using fixed dims (parallel version)
!
  USE futils
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  CHARACTER(len=256) :: file='part.h5'
  CHARACTER(len=32) :: name
  INTEGER, PARAMETER :: npart=1024*1024, nattrs=8
!!$  INTEGER, PARAMETER :: npart=1024, nattrs=2
  INTEGER :: fid, istep, nrun=1
  INTEGER :: me, npes, n, ierr, nptot, nploc, s0
  DOUBLE PRECISION :: r0=100.0, a=20.0, time, pi, dphi, s(npart), theta(npart)
  DOUBLE PRECISION :: part(nattrs, npart)   ! (r, z, phi) coordinates
  DOUBLE PRECISION :: t0, t1, twrite, mb, mbs_write
!===========================================================================
!                   1. Prologue
!
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL creatf(file, fid, 'A simple simulation', &
       &          mpicomm=MPI_COMM_WORLD, real_prec='d', mpiposix=.FALSE.)
  CALL creatg(fid, "/part", "Particles Coordinates")        ! Group
!===========================================================================
!                   2. Time loop
!
  nptot = npart*npes
  CALL dist1d(1, nptot, s0, nploc)
  CALL initp(nattrs, npart, part, s0)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  t0 = mpi_wtime()
!
  DO istep=1,nrun
     time = istep
     WRITE(name,'(a,i3.3)') "/part/", istep
     CALL putarr(fid, name, part, pardim=2) ! dataset
     CALL attach(fid, name, "time", time)   ! Attr on dataset
     CALL attach(fid, name, "step", istep)
!!$     WRITE(*,'(a,i4,a)') 'Step', istep, ' done'
  END DO
  nptot = npes*npart
  CALL attach(fid, "/part", "nptot", nptot)
  CALL attach(fid, "/part", "nattrs", nattrs)
  CALL attach(fid, "/part", "nsteps", nrun)
!
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  twrite = mpi_wtime()-t0
  mb = 8.0*REAL(SIZE(part))/1024.0/1024*npes*nrun
  mbs_write = mb/twrite
  IF( me.EQ. 0) THEN
     WRITE(*,'(a,5(f8.3,a))') 'Write performance:', &
          &   twrite, ' s,',mb, ' MB,', mbs_write, ' MB/s'
  END IF
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
  CALL mpi_finalize(ierr)
END PROGRAM main

SUBROUTINE initp(natts, n, p, s)
  IMPLICIT NONE
  INTEGER :: natts, n,s
  DOUBLE PRECISION :: x, p(natts,n)
  INTEGER :: i, j
  x=s-1.0d0
  DO j=1,n
     x=x+1.0d0
     DO i=1,natts
        p(i,j) = 1.d0/x
     END DO
  END DO
END SUBROUTINE initp

SUBROUTINE dist1d(s0, ntot, s, nloc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(in) :: s0, ntot
  INTEGER, INTENT(out) :: s, nloc
  INTEGER :: me, npes, ierr, naver, rem
!
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
  PRINT*, 'me, ntot, npes ', me, ntot, npes
  naver = ntot/npes
  rem = MODULO(ntot,npes)
  s = s0 + MIN(rem,me) + me*naver
  nloc = naver
  IF( me.LT.rem ) nloc = nloc+1
!
END SUBROUTINE dist1d
