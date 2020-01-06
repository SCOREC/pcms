PROGRAM main
!
! Parallel read of particle array
!
  USE futils
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  CHARACTER(len=32) :: filein='part.h5'
  CHARACTER(len=32) :: name
  INTEGER :: nattrs, nptot, npart
  INTEGER :: fidr
  INTEGER :: me, npes, n, ierr, nrun, s, istep, nerrs
  DOUBLE PRECISION, ALLOCATABLE  :: part(:,:)   ! (r, z, phi) coordinates
  DOUBLE PRECISION :: t0, t1, tread, mb, mbs_read
!===========================================================================
!!
  CALL openf(filein, fidr)
  CALL getatt(fidr, '/part', 'nptot', nptot)
  CALL getatt(fidr, '/part', 'nattrs', nattrs)
  CALL getatt(fidr, '/part', 'nsteps', nrun)
!
  npart = nptot
  ALLOCATE(part(nattrs,npart))
  WRITE(*,'(a,i3,2i8)') 'nattrs, npart, nrun =', nattrs, npart, nrun
!
  DO istep=1,nrun
     WRITE(name,'(a,i3.3)') "/part/", istep
     CALL getarr(fidr, name, part)
     CALL check(nattrs, npart, part, nerrs, 1)
     WRITE(*, '(a,i12)' ) 'Number of errors return from CHECK:', nerrs
     WRITE(*,'(a,i4,a)') 'Step', istep, ' done'
  END DO
!
  CALL closef(fidr)
END PROGRAM main

SUBROUTINE check(natts, n, p, nerrs, s)
  IMPLICIT NONE
  INTEGER :: natts, n, nerrs, s
  DOUBLE PRECISION :: x, p(natts,n)
  INTEGER :: i, j
  nerrs = 0
  x=0.0
  DO j=1,n
     x=x+1.0d0
     DO i=1,natts
        IF( p(i,j) .NE. 1.0d0/x ) nerrs=nerrs+1
     END DO
  END DO
END SUBROUTINE check
