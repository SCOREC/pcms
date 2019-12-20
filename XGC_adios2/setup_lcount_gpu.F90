#ifdef USE_GPU
      attributes(global)                                                 &
#endif
      subroutine setup_lcount_gpu(nx,ny,ncopy,lcount,xstart)
      use cudafor
      implicit none
      integer, value :: nx,ny,ncopy
      integer(kind=4) :: lcount(nx*ny, ncopy)
      integer :: xstart(nx*ny+1)

#ifdef USE_GPU
      attributes(device) :: lcount, xstart
#endif

      integer :: ith,iblock, gith
      integer :: nthread, nblock, gnthread 
      integer :: i,j, ip, isize

      ith = threadIdx%x + (threadIdx%y-1)*blockDim%x
      iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x
      nthread = blockDim%x * blockDim%y

      nblock = gridDim%x * gridDim%y
      gnthread = nthread * nblock
      gith = ith + (iblock-1)*nthread

      do i=gith,(nx*ny),gnthread
        ip = xstart(i)
        do j=1,ncopy
          isize = lcount(i,j)
          lcount(i,j) = ip
          ip = ip + isize
        enddo
      enddo

      return

      end subroutine setup_lcount_gpu
