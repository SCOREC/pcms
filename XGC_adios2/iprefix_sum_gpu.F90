
#ifdef USE_GPU
      attributes(global)                                                 &
#endif
      subroutine iprefix_sum_gpu_kernel(ix,n)
      use cudafor
      integer, value :: n
      integer :: ix(*)

#ifdef USE_GPU
      attributes(device) :: ix
#endif

!     --------------------------------------------
!     compute prefix sum   with starting value ival
!     [ival, ix(1)+ival, ix(1)+ix(2) + ival, + ... + sum( ix(:) ) + ival]
!     --------------------------------------------
      integer :: j, ix_j, cumulative_sum


!     ---------------------------------
!     use inefficient algorithm for now
!     ---------------------------------
      if ((threadIdx%x == 1).and.(blockIdx%x == 1) .and.                 &
          (threadIdx%y == 1).and.(blockIdx%y == 1) ) then

        cumulative_sum = 0
        do j=1,n
          ix_j = ix(j)
          cumulative_sum = cumulative_sum + ix_j
          ix(j) = cumulative_sum
        enddo
      endif

      return
      end subroutine iprefix_sum_gpu_kernel


#ifdef USE_GPU
      attributes(global)                                                &
#endif
      subroutine ishift_and_add_gpu_kernel(x_in, x_out,n,ival)
      use cudafor
      implicit none
      integer, value :: n, ival
      integer :: x_in(n), x_out(n+1)

#ifdef USE_GPU
      attributes(device) :: x_in, x_out
#endif

      integer :: i, ith, nthreads, iblock, nblocks
      integer :: gith, gnthreads

      ith = threadIdx%x + (threadIdx%y-1)*blockDim%x
      nthreads = blockDim%x * blockDim%y

      iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x
      nblocks = gridDim%x * gridDim%y

      gnthreads = nblocks * nthreads
      gith = ith + (iblock-1)*nthreads

      if (gith == 1) then
        x_out(1) = ival
      endif

      do i=gith,n,gnthreads
        x_out(i+1) = x_in(i) + ival
      enddo

      return
      end subroutine ishift_and_add_gpu_kernel

      subroutine iprefix_sum_gpu( n, ix, ival )
      use cudafor
#ifdef USE_THRUST
      use thrust_mod
#endif
      implicit none
      integer, value :: n
      integer :: ix(n+1)
      integer, value :: ival

      integer :: ix_tmp(n+1)
#ifdef USE_GPU
      attributes(device) :: ix, ix_tmp
#endif

      ix_tmp = ix

!      print*,'before prefix sum, ix_tmp ', ix_tmp

#ifdef USE_GPU

#ifdef USE_THRUST
      call thrustscan(ix_tmp,n)
#else
!     --------------------
!     very slow and simple
!     --------------------
      call iprefix_sum_gpu_kernel<<<1,1>>>(ix_tmp, n)
#endif


#else
      call iprefix_sum_gpu_kernel( ix_tmp, n )
#endif

!      print*,'after prefix sum, ix_tmp ', ix_tmp

!     ------------------------------------------------
!     shift and add
!     converrt [x1, x2, ... xn] to [ival, x1+ival, ... xn+ival]
!     ------------------------------------------------

#ifdef USE_GPU
      call ishift_and_add_gpu_kernel<<<32,192>>>(ix_tmp,ix, n, ival)
#else
      call ishift_and_add_gpu_kernel(ix_tmp,ix,n,ival)
#endif
      return
      end subroutine iprefix_sum_gpu
