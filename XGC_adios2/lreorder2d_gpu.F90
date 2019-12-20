      subroutine lreorder2d_gpu(m,n,A,lda, iperm, streamid_in)
      use cudafor
      implicit none

      integer,intent(in) :: m,n,lda
      integer*8, target :: A(lda,*)
      integer :: iperm(*)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      integer*8 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: j, ierr

      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int(m/nthread)+1)
      do j=1,n
        call lreorder_gpu<<<nblock,nthread,0,streamid>>>(m,              &
     &                    A(1,j), tmp, iperm)

        if (present(streamid_in)) then
          ierr = cudaMemcpyAsync(A(1,j), tmp, m,                         &
     &               cudaMemcpyDeviceToDevice,streamid )
          call assert(ierr.eq.0,'lreorder2d_gpu:cudaMemcpyAsync',ierr)
        else
          ierr = cudaMemcpy(A(1,j), tmp, m )
          call assert(ierr.eq.0,'lreorder2d_gpu:cudaMemcpy',ierr)
        endif
      enddo

      return
      end subroutine lreorder2d_gpu

       
