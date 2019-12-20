#ifdef USE_GPU
       attributes(global)                                                &
#endif
       subroutine isum_col_gpu( m,n, lda, lcount, icount )
       use cudafor
       implicit none 
       integer, value :: m,n, lda
       integer :: lcount(lda,n)
       integer :: icount(m)

#ifdef USE_GPU
       attributes(device) :: icount, lcount
#endif

       integer :: i,j, icount_local
       integer :: gith,ith, iblock, nthread, gnthread

       ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
       iblock = blockIdx%x + (blockIdx%y - 1)*gridDim%x
       nthread = blockDim%x * blockDim%y
       gnthread = nthread * ( gridDim%x * gridDim%y )
       gith = ith + (iblock-1)*nthread

       do i=gith,m,gnthread
           icount_local  = 0
           do j=1,n
              icount_local = icount_local + lcount(i, j)
           enddo
           icount(i) = icount_local
       enddo
       return
       end subroutine isum_col_gpu
