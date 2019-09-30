PROGRAM main
!
!   Parallel write/read of 2d and 3d arrays partionned on
!   2d processor grid
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5'
  INTEGER :: ierr, fid, me, npes, comm=MPI_COMM_WORLD
!
  INTEGER, PARAMETER :: ndims=2
  INTEGER, DIMENSION(ndims) :: dims, coords
  LOGICAL :: periods(ndims), reorder
  INTEGER :: cart, cartcol, cartrow
!
  INTEGER, DIMENSION(ndims) :: offsets, np
  INTEGER :: n1=7, n2=9, n3=8
  INTEGER :: i, j, k, iglob, jglob, kglob
  DOUBLE PRECISION, ALLOCATABLE :: array2(:,:), array3(:,:,:), brray3(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: carray2(:,:), carray3(:,:,:), cbrray3(:,:,:)
  DOUBLE PRECISION :: a
  DOUBLE COMPLEX :: za
  INTEGER :: nerrors
!===========================================================================
!                             1.  Prologue
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Create cartesian topololy
!
  dims(1) = 3; dims(2:ndims) = 0;
  periods = (/.FALSE., .FALSE./)
  reorder = .FALSE.
  CALL mpi_dims_create(npes, ndims, dims, ierr)
  CALL mpi_cart_create(comm, ndims, dims, periods, reorder, cart, ierr)
  CALL mpi_cart_coords(cart, me, ndims, coords, ierr)
  CALL mpi_cart_sub(cart, (/.TRUE., .FALSE. /), cartcol, ierr)
  CALL mpi_cart_sub(cart, (/.FALSE., .TRUE. /), cartrow, ierr)
  IF( me .EQ. 0 ) WRITE(*,'(a,i2,i2)') 'Processor grid', dims(1), dims(2)
!
!   Create file collectively
!
  CALL creatf(file, fid, &
       &      desc="A parallel file", &
       &      real_prec='s', &
       &      mpiposix=.FALSE., &
       &      mpicomm=cart)
!===========================================================================
!                             2.  Parallel write file
!
!   Define local 2d array
!
  CALL dist1d(cartcol, 0, n1, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n2, offsets(2), np(2))
!
  WRITE(*,'(a,i3.3,a,10i5)') 'PE', me, ': coords, offsets, np', &
       &                      coords, offsets, np
!
  ALLOCATE(array2(np(1),np(2)))
  ALLOCATE(carray2(np(1),np(2)))
  DO i=1,np(1)
     iglob = offsets(1)+i
     DO j=1,np(2)
        jglob = offsets(2)+j
        array2(i,j) = 10*iglob + jglob
     END DO
  END DO
  carray2 = CMPLX(array2, -array2, KIND(array2))
!
!   Define local 3d array A(n1, n2/P1, n3/P2)
!
  CALL dist1d(cartcol, 0, n2, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n3, offsets(2), np(2))
!
  WRITE(*,'(a,i3.3,a,10i5)') 'PE', me, ': coords, offsets, np', &
       &                      coords, offsets, np
!
  ALLOCATE(array3(n1,np(1),np(2)))
  ALLOCATE(carray3(n1,np(1),np(2)))
  DO i=1,n1
     DO j=1,np(1)
        jglob = offsets(1)+j
        DO k=1,np(2)
           kglob = offsets(2)+k
           array3(i,j,k) = 100*i + 10*jglob + kglob
        END DO
     END DO
  END DO
  carray3 = CMPLX(-array3/10.0d0, array3, KIND(array3))
!
!   Define local 3d array B(n1/P1, n2, n3/P2)
!
  CALL dist1d(cartcol, 0, n1, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n3, offsets(2), np(2))
!
  WRITE(*,'(a,i3.3,a,10i5)') 'PE', me, ': coords, offsets, np', &
       &                      coords, offsets, np
!
  ALLOCATE(brray3(np(1),n2,np(2)))
  ALLOCATE(cbrray3(np(1),n2,np(2)))
  DO i=1,np(1)
     iglob = offsets(1)+i
     DO j=1,n2
        DO k=1,np(2)
           kglob = offsets(2)+k
           brray3(i,j,k) = 100*iglob + 10*j + kglob
        END DO
     END DO
  END DO
  cbrray3 = CMPLX(-brray3/10.0, brray3)
!
!   Write to file
!
  CALL putarrnd(fid, "/parray2", array2, (/1,2/))
  CALL putarrnd(fid, "/pcarray2", carray2, (/1,2/))
!
  CALL putarrnd(fid, "/parray3", array3, (/2,3/))
  CALL putarrnd(fid, "/pcarray3", carray3, (/2,3/))
!
  CALL putarrnd(fid, "/parray3b", brray3, (/1,3/))
  CALL putarrnd(fid, "/pcarray3b", cbrray3, (/1,3/))
  WRITE(*,'(a,i3.3,a)') 'PE', me, ': All write ok!'
!===========================================================================
!                             3.  Parallel read file
!
!   Close and reopen file
  CALL closef(fid)
  CALL openf(file, fid, mpicomm=cart)
!
!   Read file
  array2 = 0.0d0
  array3 = 0.0d0
  carray2 = (0.0d0, 0.0d0)
  carray3 = (0.0d0, 0.0d0)
  CALL getarrnd(fid, "/parray2", array2, (/1,2/))
  CALL getarrnd(fid, "/pcarray2", carray2, (/1,2/))
  CALL getarrnd(fid, "/parray3", array3, (/2,3/))
  CALL getarrnd(fid, "/pcarray3", carray3, (/2,3/))
!
!  Check read arrays
  CALL dist1d(cartcol, 0, n1, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n2, offsets(2), np(2))
  nerrors = 0
  DO i=1,np(1)
     iglob = offsets(1)+i
     DO j=1,np(2)
        jglob = offsets(2)+j
        a = 10*iglob + jglob
        IF( a .NE. array2(i,j) ) nerrors = nerrors+1
     END DO
  END DO
  WRITE(*,'(a25,i5)') 'nerrors reading /parray2', nerrors
!
  nerrors = 0
  DO i=1,np(1)
     iglob = offsets(1)+i
     DO j=1,np(2)
        jglob = offsets(2)+j
        za = CMPLX(array2(i,j), -array2(i,j))
        IF( za .NE. carray2(i,j) ) nerrors = nerrors+1
     END DO
  END DO
  WRITE(*,'(a25,i5)') 'nerrors reading /pcarray2', nerrors
!
  CALL dist1d(cartcol, 0, n2, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n3, offsets(2), np(2))
  nerrors = 0
  DO i=1,n1
     DO j=1,np(1)
        jglob = offsets(1)+j
        DO k=1,np(2)
           kglob = offsets(2)+k
           a = 100*i + 10*jglob + kglob
           IF( a .NE. array3(i,j,k) )  nerrors = nerrors+1
        END DO
     END DO
  END DO
  WRITE(*,'(a25,i5)') 'nerrors reading /parray3', nerrors
!
  nerrors = 0
  DO i=1,n1
     DO j=1,np(1)
        jglob = offsets(1)+j
        DO k=1,np(2)
           kglob = offsets(2)+k
           za = CMPLX(-array3(i,j,k)/10.0, array3(i,j,k))
           IF( za .NE. carray3(i,j,k) )  nerrors = nerrors+1
        END DO
     END DO
  END DO
  WRITE(*,'(a25,i5)') 'nerrors reading /pcarray3', nerrors
!===========================================================================
!                             9.  Epilogue
!
  DEALLOCATE(array2, array3, brray3)
  DEALLOCATE(carray2, carray3, cbrray3)
!
  CALL closef(fid)
  CALL mpi_finalize(ierr)
END PROGRAM main

SUBROUTINE dist1d(comm, s0, ntot, s, nloc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(in) :: s0, ntot
  INTEGER, INTENT(out) :: s, nloc
  INTEGER :: comm, me, npes, ierr, naver, rem
!
  CALL MPI_COMM_SIZE(comm, npes, ierr)
  CALL MPI_COMM_RANK(comm, me, ierr)
  naver = ntot/npes
  rem = MODULO(ntot,npes)
  s = s0 + MIN(rem,me) + me*naver
  nloc = naver
  IF( me.LT.rem ) nloc = nloc+1
!
END SUBROUTINE dist1d
