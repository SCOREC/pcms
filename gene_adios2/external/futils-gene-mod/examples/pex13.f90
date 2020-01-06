PROGRAM main
!
!  Read array from file created by pex12
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5', str
  INTEGER :: ierr, fid, me, npes, lstr
  INTEGER, PARAMETER :: ndims=2, comm=MPI_COMM_WORLD
  INTEGER, DIMENSION(ndims) :: dims, coords
  LOGICAL :: periods(ndims), reorder
  INTEGER :: cart, cartcol, cartrow
!
  INTEGER, DIMENSION(ndims) :: offsets, np
  INTEGER :: n1=7, n2=9, n3=3, lun
  INTEGER :: i, j, iglob, jglob, k, nerrors
  DOUBLE PRECISION, ALLOCATABLE :: array2(:,:,:)
  DOUBLE PRECISION :: a
!===========================================================================
!                             1.  Prologue
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Read 2d processor grid P1xP2 from command line
  IF( command_argument_count() .EQ. 2 ) THEN
     CALL get_command_argument(1, str, lstr, ierr)
     READ(str(1:lstr),'(i3)') dims(1)
     CALL get_command_argument(2, str, lstr, ierr)
     READ(str(1:lstr),'(i3)') dims(2)
  ELSE
     IF( me.EQ.0 ) PRINT*, 'Requires 2 arguments: P1, P2!'
     CALL mpi_abort(MPI_COMM_WORLD, -1, ierr)
  END IF
!
!   Create cartesian topololy
  IF( me .EQ. 0 ) WRITE(*,'(a,i3,i3)') '2d rocessor grid', dims
  periods = (/.FALSE., .FALSE./)
  reorder = .FALSE.
  CALL mpi_cart_create(comm, ndims, dims, periods, reorder, cart, ierr)
  CALL mpi_cart_coords(cart, me, ndims, coords, ierr)
  CALL mpi_cart_sub(cart, (/.TRUE., .FALSE. /), cartcol, ierr)
  CALL mpi_cart_sub(cart, (/.FALSE., .TRUE. /), cartrow, ierr)
!===========================================================================
!                             2.   Array parallel partition
!
!   Local array offsets and sizes
  CALL dist1d(cartcol, 0, n1, offsets(1), np(1))
  CALL dist1d(cartrow, 0, n2, offsets(2), np(2))
  WRITE(*,'(a,i3.3,a,10i5)') 'PE', me, ': coords, offsets, np', &
       &                      coords, offsets, np
!
!   Allocate local arrays including ghost cells
  ALLOCATE(array2(0:np(1)+1, 0:np(2)+1, n3))
!
!===========================================================================
!                             3.   Parallel read
!
  CALL openf(file, fid,  mpicomm=cart)
  array2=-1   ! Initialize with -1
  CALL getarrnd(fid, "/parray2", array2, (/1,2/), garea=(/1,1/))
  lun = 90+me
  WRITE(lun,*) 'Matrice on proc,', me
  DO k=1,n3
     WRITE(lun,'(a,i2)') 'k =', k
     DO i=0,np(1)+1
        WRITE(lun,'(20f6.0)') array2(i,:, k)
     END DO
  END DO
!
!   Check read array
!
  nerrors = 0
  DO i=1,np(1)
     iglob = offsets(1)+i
     DO j=1,np(2)
        jglob = offsets(2)+j
        DO k=1,n3
           a = 100*iglob + 10*jglob + k
           IF( a .NE. array2(i,j,k) ) nerrors = nerrors+1
        END DO
     END DO
  END DO
  WRITE(*,'(a,6i5)') 'nerrors reading /parray2', nerrors
  CALL closef(fid)
!===========================================================================
!                             9.  Epilogue
  DEALLOCATE(array2)
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
