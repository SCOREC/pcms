PROGRAM main
!
!   Write and read  complex arrays
!
  USE futils
  IMPLICIT NONE
  INTEGER, PARAMETER :: d1=10, d2=4, d3=3
  CHARACTER(len=32) :: file='complex.h5'
  INTEGER :: fid, i, j, k, ierr
  COMPLEX :: carr(d1), carr2(d1,d2), carr3(d1,d2,d3), cval
  DOUBLE COMPLEX :: zarr(d1), zarr2(d1,d2), zarr3(d1,d2,d3), zval
  DOUBLE PRECISION :: time
  INTEGER :: rank, dims(7), nprev
  DOUBLE COMPLEX :: cphi(d1), cscal
!===========================================================================
!                   1. Prologue
!
  CALL creatf(file, fid, 'Complex putarr and getarr', real_prec='d')
!!$  CALL creatf(file, fid, 'Complex putarr and getarr')
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' created'
!
  CALL creatg(fid, "/1d", "1d arrays")        ! Group
  CALL creatg(fid, "/2d", "2d arrays")
  CALL creatg(fid, "/3d", "3d arrays")
!
!===========================================================================
!                   2. Create datasets
!
!  Fixed dims with PUTARR
!
  DO i=1,d1
     carr(i) = CMPLX(i,-i)
     zarr(i) = 10*CMPLX(i,-i)
  END DO
  CALL putarr(fid, '/1d/carr', carr)
  CALL putarr(fid, '/1d/zarr', zarr)
!
  DO i=1,d1
     DO j=1,d2
        carr2(i,j) = CMPLX(i,-j)
        zarr2(i,j) = CMPLX(i,-j)/10.0
     END DO
  END DO
  CALL putarr(fid, '/2d/carr', carr2)
  CALL putarr(fid, '/2d/zarr', zarr2)
!
  DO i=1,d1
     DO j=1,d2
        DO k=1,d3
           carr3(i,j,k) = CMPLX(i,-j)*10*k
           zarr3(i,j,k) = CMPLX(-i,j)/(10*k)
        END DO
     END DO
  END DO
  CALL putarr(fid, '/3d/carr', carr3)
  CALL putarr(fid, '/3d/zarr', zarr3)
!
!   Extendible dims with created/append
!
  rank = 0
  CALL creatd(fid, rank, dims, '/1d/time', 'Time')
  CALL creatd(fid, rank, dims, '/1d/cscal', '0d array', iscomplex=.TRUE.)
  rank=1; dims(1:rank) = SHAPE(cphi)
  CALL creatd(fid, rank, dims, '/1d/cphi', 'Complex field', iscomplex=.TRUE.)
  DO i=1,5
     time = i-1
     cscal = CMPLX(10*i, -i)
     cphi(1:d1) = CMPLX(i, 2*i)
     CALL append(fid, '/1d/time', time)
     CALL append(fid, '/1d/cscal', cscal)
     CALL append(fid, '/1d/cphi', cphi)
  END DO
!
!   Use append to write zarr1
!
  CALL creatd(fid, 0, (/d1/), '/1d/zarr_append', iscomplex=.TRUE.)
  DO i=1,d1
     CALL append(fid, '/1d/zarr_append', zarr(i))
  END DO
!
!   Use append to write zarr2
!
  CALL creatd(fid, 1, (/d1/), '/2d/zarr_append', iscomplex=.TRUE.)
  DO j=1,d2
     CALL append(fid, '/2d/zarr_append', zarr2(:,j))
  END DO
!
!   Use append to write zarr3
!
  CALL creatd(fid, 2, (/d1,d2/), '/3d/zarr_append', iscomplex=.TRUE.)
  DO k=1,d3
     CALL append(fid, '/3d/zarr_append', zarr3(:,:,k))
  END DO
!
  CALL closef(fid)
!===========================================================================
!                   3. Read and check datasets
!
  carr = (0., 0.)
  zarr = (0., 0.)
  CALL openf(file, fid)
  CALL getarr(fid, '/1d/carr', carr)
  CALL getarr(fid, '/1d/zarr', zarr)
  ierr =0
  DO i=1,d1
     cval = CMPLX(i,-i); IF( carr(i) .NE. cval ) ierr = ierr+1
     zval =10*CMPLX(i,-i); IF( zarr(i) .NE. zval ) ierr = ierr+1
  END DO
  PRINT*, 'Check carr/zarr, number of errors', ierr
!
  carr2 = (0., 0.)
  zarr2 = (0., 0.)
  CALL getarr(fid, '/2d/carr', carr2)
  CALL getarr(fid, '/2d/zarr', zarr2)
  ierr =0
  DO i=1,d1
     DO j=1,d2
        cval = CMPLX(i,-j); IF( carr2(i,j) .NE. cval ) ierr = ierr+1
        zval = CMPLX(i,-j)/10.0; IF( zarr2(i,j) .NE. zval ) ierr = ierr+1
     END DO
  END DO
  PRINT*, 'Check carr2/zarr2, number of errors', ierr
!
  carr3 = (0., 0.)
  zarr3 = (0., 0.)
  CALL getarr(fid, '/3d/carr', carr3)
  CALL getarr(fid, '/3d/zarr', zarr3)
  ierr =0
  DO i=1,d1
     DO j=1,d2
        DO k=1,d3
           cval = CMPLX(i,-j)*10*k; IF( carr3(i,j,k) .NE. cval ) ierr = ierr+1
           zval = CMPLX(-i,j)/(10*k); IF( zarr3(i,j,k) .NE. zval ) ierr = ierr+1
        END DO
     END DO
  END DO
  PRINT*, 'Check carr3/zarr3, number of errors', ierr
!
!   Append to the existent extendible datasets
  CALL getsize(fid, '/1d/time', nprev)
  WRITE(*,'(a,i6)') 'Number of previous steps ', nprev
  DO i=nprev+1,nprev+5
     time = i-1
     cscal = CMPLX(10*i, -i)
     cphi(1:d1) = CMPLX(i, 2*i)
     CALL append(fid, '/1d/time', time)
     CALL append(fid, '/1d/cscal', cscal)
     CALL append(fid, '/1d/cphi', cphi)
  END DO
  CALL getsize(fid, '/1d/time', nprev)
  WRITE(*,'(a,i6)') 'Number of total steps', nprev
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
END PROGRAM main
