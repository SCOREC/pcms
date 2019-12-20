       attributes(global) &
     & subroutine sreorder_gpu(m, xin, xout, iperm)
       use cudafor
       integer, value :: m
       real*4 :: xin(m),xout(m)
       integer ::iperm(m)
       attributes(device) :: xin, xout, iperm

       integer :: i, ith, iblock, nthread, gnthreads, gith

       ith = threadIdx%x + (threadIdx%y-1)*blockDim%x
       iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x
       nthread = blockDim%x * blockDim%y

       gnthreads = nthread * (gridDim%x * gridDim%y)
       gith = ith + (iblock-1)*nthread

       do i=gith,m,gnthreads
         xout(i) = xin( iperm(i) )
       enddo

       return
       end subroutine sreorder_gpu


          
         
