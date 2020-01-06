PROGRAM main
!
!   Parallel write/read a 3d array partionned on 2d processor grid
!   A(n1/P1, n2/P2, n3), with GHOST CELLS on the partitionned dimensions.
!
  USE futils
  IMPLICIT NONE
  INCLUDE "mpif.h"
  CHARACTER(len=32) :: file='para.h5'
  INTEGER :: ierr, fid, me, npes
  INTEGER, PARAMETER :: ndims=2
  INTEGER, PARAMETER ::  n1p=3, n2p=2, n3=2  ! Dimension of local array
  REAL, DIMENSION(0:n1p+1,0:n2p+1,n3) :: arrayg      ! Local array including ghost area
  COMPLEX, DIMENSION(0:n1p+1,0:n2p+1,n3) :: carrayg
  INTEGER, DIMENSION(ndims) :: dims, coords
  LOGICAL :: periods(ndims), reorder
  INTEGER :: cart, i, j, k, iglob, jglob, nerrors(2)
  REAL :: a
  COMPLEX :: ca
!===========================================================================
!                             1.  Prologue
!
!   Init MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
!
!   Create cartesian topololy
!
  dims    = (/2, 4/)
  periods = (/.FALSE., .TRUE./)
  reorder = .FALSE.
  IF( PRODUCT(dims) .NE. npes ) THEN
     IF( me .EQ. 0 ) THEN
        PRINT*,  PRODUCT(dims), " processors required!"
        CALL mpi_abort(MPI_COMM_WORLD, -1, ierr)
     END IF
  END IF
  CALL mpi_cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, cart, ierr)
  CALL mpi_cart_coords(cart, me, ndims, coords, ierr)
!
!   Define local array
!
  arrayg = me
  DO i=1,n1p
     iglob = coords(1)*n1p + i
     DO j=1,n2p
        jglob = coords(2)*n2p + j
        DO k=1,n3
           arrayg(i,j,k) = 100*iglob + 10*jglob + k
        END DO
     END DO
  END DO
  carrayg = CMPLX(arrayg, REAL(me+1))
!===========================================================================
!                             2.  Parallel write file
!
!   Create file collectively, passing the comm. with cartesian topopily
  CALL creatf(file, fid, mpicomm=cart)
!
!   Write to file collectively using "nd" version of "putarr".
  CALL putarrnd(fid, '/parray0',  arrayg,  (/1,2/), &
       &        desc='With ghost area')
  CALL putarrnd(fid, '/parray',  arrayg,  (/1,2/), garea=(/1,1/), &
       &        desc='Without ghost area')
  CALL putarrnd(fid, '/pcarray', carrayg, (/1,2/), garea=(/1,1/), &
       &        desc='Without ghost area')
!===========================================================================
!                             3.  Parallel read file
!
!   Close and reopen file
  CALL closef(fid)
  CALL openf(file, fid, mpicomm=cart)
!
!   Read file
  arrayg = 0
  carrayg = 0
  CALL getarrnd(fid, '/parray',  arrayg,  (/1,2/), garea=(/1,1/))
  CALL getarrnd(fid, '/pcarray', carrayg, (/1,2/), garea=(/1,1/))
!
!  Check read arrays
  nerrors = 0
  DO i=1,n1p
     iglob = coords(1)*n1p + i
     DO j=1,n2p
        jglob = coords(2)*n2p + j
        DO k=1,n3
           a = 100*iglob + 10*jglob + k
           ca = CMPLX(a, REAL(me+1))
           IF( a .NE. arrayg(i,j,k))   nerrors(1) = nerrors(1)+1
           IF( ca .NE. carrayg(i,j,k)) nerrors(2) = nerrors(2)+1
        END DO
     END DO
  END DO
  WRITE(*,'(a,6i5)') 'nerrors reading /parray and /pcarray', nerrors
!
!===========================================================================
!                             9.  Epilogue
!
!   Clean up and quit
  CALL closef(fid)
  CALL mpi_finalize(ierr)
!
END PROGRAM main
