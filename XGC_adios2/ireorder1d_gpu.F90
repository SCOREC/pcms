      subroutine ireorder1d_gpu(m,A, iperm, streamid_in )
      use cudafor
      integer :: m
      integer :: A(m)
      integer :: iperm(m)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      integer :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer:: nblock = 32
      integer :: ierr


      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int( m/nthread ) + 1)

      call ireorder_gpu<<<nblock,nthread,0,streamid>>>(m,                &
     &                       A, tmp, iperm)

      if (present(streamid_in)) then
        ierr = cudaMemcpyAsync(A, tmp, m,                                &
     &               cudaMemcpyDeviceToDevice, streamid )
        call assert(ierr.eq.0,'ireorder1d_gpu:cudaMemcpyAsync',ierr)
      else
        ierr = cudaMemcpy(A, tmp, m )
        call assert(ierr.eq.0,'ireorder1d_gpu:cudaMemcpy',ierr)
      endif

      return
      end subroutine ireorder1d_gpu

       
