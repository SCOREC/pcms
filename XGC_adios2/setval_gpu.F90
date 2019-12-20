       attributes(global) &
       subroutine setval_kernel_gpu( n, v, dvalue)
       use precision_mod_gpu
       use cudafor
       implicit none

       integer, value :: n
       real(kind=work_p), value :: dvalue
       real(kind=work_p) :: v(*)


       integer :: i, ith, nthreads

       ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) +  &
               ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
                     (blockDim%x * blockDim%y)

       nthreads = (blockDim%x * blockDim%y)*(gridDim%x * gridDim%y)
       do i=ith,n,nthreads
         v(i) = dvalue
       enddo
       return
       end subroutine setval_kernel_gpu


       attributes(host) &
       subroutine setval_gpu( n, v, dvalue, streamid_in )
       use precision_mod_gpu
       use cudafor
       implicit none

       integer, intent(in) :: n
       real(kind=work_p), intent(in) :: dvalue
       real(kind=work_p), intent(inout) :: v(*)
       attributes(device) :: v


       integer, intent(in), optional :: streamid_in
       integer :: streamid

       type(dim3) :: tgrid, tblock

       integer, parameter :: nthread = 256

       tblock%x = nthread
       tblock%y = 1
       tblock%z = 1

       tgrid%x = min(32*1024, int(n/nthread)+1)
       tgrid%y = 1
       tgrid%z = 1

       streamid = 0
       if (present(streamid_in)) streamid = streamid_in

       call setval_kernel_gpu<<<tgrid,tblock,0,streamid>>>(n,v,dvalue)

       return
       end subroutine setval_gpu

