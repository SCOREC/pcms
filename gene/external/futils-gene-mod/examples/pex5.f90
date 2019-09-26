PROGRAM main
!
! Parallel read of particle array
!
  USE futils
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  CHARACTER(len=256) :: file='part.h5'
  CHARACTER(len=32) :: name
  INTEGER :: nattrs, nptot, npart
  INTEGER :: fidr
  INTEGER :: me, npes, n, ierr, nrun, s, istep, nerrs, i
  DOUBLE PRECISION, ALLOCATABLE  :: part(:,:)   ! (r, z, phi) coordinates
  DOUBLE PRECISION :: t0, t1, tread, mb, mbs_read
!===========================================================================
!!
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Parallel read of "part.h5"
  CALL openf(file, fidr, mpicomm=MPI_COMM_WORLD)
!
  CALL getatt(fidr, '/part', 'nptot', nptot)
  CALL getatt(fidr, '/part', 'nattrs', nattrs)
  CALL getatt(fidr, '/part', 'nsteps', nrun)
!
  CALL dist1d(1, nptot, s, npart)
  ALLOCATE(part(nattrs,npart))
!
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  t0 = mpi_wtime()
!
  DO istep=1,nrun
     WRITE(name,'(a,i3.3)') "/part/", istep
     CALL getarr(fidr, name, part, pardim=2)
  END DO
!
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  tread = mpi_wtime()-t0
  mb = 8.0*REAL(SIZE(part))/1024.0/1024*npes*nrun
  mbs_read = mb/tread
  IF( me.EQ. 0) THEN
     WRITE(*,'(a,5(f8.3,a))') 'Read performance:', &
          &   tread, ' s,',mb, ' MB,', mbs_read, ' MB/s'
  END IF
!
  CALL check(nattrs, npart, part, nerrs, s)
  n = nerrs
  CALL mpi_reduce(n, nerrs, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  IF ( me .EQ. 0 ) THEN
     WRITE(*, '(a,i12)' ) 'Number of errors return from CHECK:', nerrs
  END IF
!
  CALL closef(fidr)
  CALL mpi_finalize(ierr)
END PROGRAM main

SUBROUTINE dist1d(s0, ntot, s, nloc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(in) :: s0, ntot
  INTEGER, INTENT(out) :: s, nloc
  INTEGER :: me, npes, ierr, naver, rem
!
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
  naver = ntot/npes
  rem = MODULO(ntot,npes)
  s = s0 + MIN(rem,me) + me*naver
  nloc = naver
  IF( me.LT.rem ) nloc = nloc+1
!
END SUBROUTINE dist1d

SUBROUTINE check(natts, n, p, nerrs, s)
  IMPLICIT NONE
  INTEGER :: natts, n, nerrs, s
  DOUBLE PRECISION :: x, p(natts,n)
  INTEGER :: i, j
  nerrs = 0
  x=s-1.0d0
  DO j=1,n
     x=x+1.0d0
     DO i=1,natts
        IF( p(i,j) .NE. 1.0d0/x ) nerrs=nerrs+1
     END DO
  END DO
END SUBROUTINE check
