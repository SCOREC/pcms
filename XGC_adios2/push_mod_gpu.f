! collection of codes related to GPU
! try bottom up approach
!
module precision_mod_gpu

! single precision
integer, parameter :: single_p = kind(1.0)

! double precision
integer, parameter :: double_p = kind(1.0d0)

! working precision, set wp = dp or wp = sp

integer, parameter :: work_p = double_p
integer, parameter :: eq_mpsi = 151
integer, parameter :: eq_mr = 151
integer, parameter :: eq_mz = 151
end module precision_mod_gpu



module util_mod_gpu
contains
      integer function get_gpu_streamid()
      use cudafor

      integer, allocatable, dimension(:) :: streamid_array
      integer :: ierr
      save

! ----------------------
! use the default stream
! ----------------------
      get_gpu_streamid = 0
      return
      end function get_gpu_streamid

       

       attributes(global) &
       subroutine setval_kernel_gpu( n, v, dvalue)
       use precision_mod_gpu
       use cudafor
       implicit none

       integer, value :: n
       real(kind=work_p), value :: dvalue
       real(kind=work_p) :: v(*)


       integer :: i, ith, nthreads

       ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) + &
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

        attributes(global) &
        subroutine sum2_kernel_gpu(m,n,A,Asum)
        use cudafor
        use precision_mod_gpu
        implicit none

        integer, value :: m,n
        real(kind=work_p), intent(in) :: A(m,n)
        real(kind=work_p), intent(inout) :: Asum(m)
        attributes(device) :: A, Asum

        integer :: i,j, ith,nthreads
        real(kind=work_p) :: dsum

        ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) + &
                ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
                      (blockDim%x * blockDim%y)

        i = ith
        if ((1 <= i) .and. (i <= m)) then
          dsum = 0
          do j=1,n
            dsum = dsum + A(i,j)
          enddo
          Asum(i) = dsum
        endif
        

        return
        end subroutine sum2_kernel_gpu




        attributes(host) &
        subroutine sum2_gpu(m,n, A_gpu, Asum, streamid_in )
        use cudafor
        use precision_mod_gpu
        implicit none

        integer, intent(in) :: m,n
        real(kind=work_p), intent(in) :: A_gpu(m,n)
        real(kind=work_p), intent(inout) :: Asum(m)

        integer, intent(in), optional :: streamid_in
        integer :: streamid

        real(kind=work_p) :: Asum_gpu(m)
        attributes(device) :: Asum_gpu
        attributes(device) :: A_gpu

        real(kind=work_p), parameter :: dzero = 0.0d0
        type(dim3) :: tgrid, tblock

        tblock%x = 256
        tblock%y = 1
        tblock%z = 1

        tgrid%x = int( m/tblock%x ) + 1
        tgrid%y = 1
        tgrid%z = 1

        if (tgrid%x > 64*1024) then
          tgrid%y = int(sqrt( dble(tgrid%x) )) + 1
          tgrid%x = tgrid%y
        endif


! print*,'m = ', m, 'tgrid%x ',tgrid%x, 'tgrid%y ',tgrid%y


        streamid = 0
        if (present(streamid_in)) streamid = streamid_in
        call sum2_kernel_gpu<<<tgrid,tblock,0,streamid>>>(m,n, &
     & A_gpu,Asum_gpu)
        Asum = Asum_gpu
        return
        end subroutine sum2_gpu
        attributes(global) &
        subroutine sum4_kernel_gpu(l1,u1,l2,u2,l3,u3,l4,u4,A,Asum)
        use precision_mod_gpu
        use cudafor
        implicit none

        integer, value :: l1,u1,l2,u2,l3,u3,l4,u4
        real(kind=work_p) :: A(l1:u1,l2:u2,l3:u3,l4:u4)
        real(kind=work_p) :: Asum(l1:u1,l2:u2,l3:u3)
        attributes(device) :: A, Asum

        integer :: i1,i2,i3,i4
        logical :: has_work
        real(kind=work_p) :: dsum



        i1 = (threadIdx%x -1) + l1
        i2 = (blockIdx%x -1) + l2
        i3 = (blockIdx%y -1) + l3

        has_work = (l1 <= i1).and.(i1 <= u1) .and. &
                   (l2 <= i2).and.(i2 <= u2) .and. &
                   (l3 <= i3).and.(i3 <= u3)
        if (has_work) then
            dsum = 0.0d0
            do i4=l4,u4
              dsum = dsum + A(i1,i2,i3,i4)
            enddo
            !print*,'i1,i2,i3, dsum ',i1,i2,i3,dsum
            Asum(i1,i2,i3) = dsum
        endif

        return
        end subroutine sum4_kernel_gpu




        attributes(host) &
        subroutine sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4,A_gpu, &
     & Asum, streamid_in )
        use precision_mod_gpu
        use cudafor
        implicit none

        integer, intent(in) :: l1,u1,l2,u2,l3,u3,l4,u4
        real(kind=work_p) :: A_gpu(l1:u1,l2:u2,l3:u3,l4:u4)
        real(kind=work_p) :: Asum(l1:u1,l2:u2,l3:u3)

        integer, intent(in), optional :: streamid_in
        integer :: streamid

        real(kind=work_p), allocatable, dimension(:,:,:) :: Asum_gpu
        integer :: n1,n2,n3

        attributes(device) :: Asum_gpu
        attributes(device) :: A_gpu


        type(dim3) :: tgrid, tblock


        n1 = u1-l1+1
        n2 = u2-l2+1
        n3 = u3-l3+1

        tblock%x = n1
        tblock%y = 1
        tblock%z = 1

        tgrid%x = n2
        tgrid%y = n3
        tgrid%z = 1


        allocate( Asum_gpu(l1:u1,l2:u2,l3:u3) )
        Asum_gpu = 0

! print*,'before sum4_kernel_gpu'
        streamid = 0
        if (present(streamid_in)) streamid = streamid_in
        call sum4_kernel_gpu<<<tgrid,tblock,0,streamid>>>(l1,u1,l2,u2,l3,u3,l4,u4, &
                                       A_gpu,Asum_gpu)
! print*,'after sum4_kernel_gpu '
        Asum = Asum_gpu

! print*,'after Asum  = Asum_gpu '

        deallocate( Asum_gpu )
! print*,'after deallocate(Asum_gpu) '
        return
        end subroutine sum4_gpu
end module util_mod_gpu

! --------------------------------
! code related to particle sorting 
! --------------------------------
      module thrust_mod

      interface thrustsort
        subroutine sort_int(input,N) &
     & bind(C,name="sort_int_wrapper")
        use iso_c_binding
        integer(c_int),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine sort_float(input,N) &
     & bind(C,name="sort_float_wrapper")
        use iso_c_binding
        real(c_float),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine sort_double( input,N) &
     & bind(C,name="sort_double_wrapper")
        use iso_c_binding
        real(c_double),device:: input(*)
        integer(c_int),value:: N
        end subroutine

      end interface


      interface thrustscan
        subroutine scan_int(input,N) &
     & bind(C,name="scan_int_wrapper")
        use iso_c_binding
        integer(c_int),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine scan_float(input,N) &
     & bind(C,name="scan_float_wrapper")
        use iso_c_binding
        real(c_float),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine scan_double( input,N) &
     & bind(C,name="scan_double_wrapper")
        use iso_c_binding
        real(c_double),device:: input(*)
        integer(c_int),value:: N
        end subroutine
      end interface thrustscan

      end module thrust_mod


       module reorder_gpu_mod
       use cudafor
       implicit none


       interface reorder_gpu
         module procedure ireorder_gpu, sreorder_gpu, &
     & dreorder_gpu, lreorder_gpu
       end interface 

       interface reorder1d_gpu
         module procedure ireorder1d_gpu, sreorder1d_gpu, &
     & dreorder1d_gpu, lreorder1d_gpu
       end interface 

       interface reorder2d_gpu
         module procedure ireorder2d_gpu, sreorder2d_gpu, &
     & dreorder2d_gpu, lreorder2d_gpu
       end interface 

       contains

       attributes(global) &
     & subroutine lreorder_gpu(m, xin, xout, iperm)
       use cudafor
       integer, value :: m
       integer*8 :: xin(m),xout(m)
       integer :: iperm(m)
       attributes(device) :: xin,xout,iperm

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
       end subroutine lreorder_gpu
       attributes(global) &
     & subroutine ireorder_gpu(m, xin, xout, iperm)
       use cudafor
       integer, value :: m
       integer :: xin(m),xout(m)
       integer :: iperm(m)
       attributes(device) :: xin,xout,iperm

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
       end subroutine ireorder_gpu


          
         
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


          
         
       attributes(global) &
     & subroutine dreorder_gpu(m, xin, xout, iperm)
       use cudafor
       integer, value :: m
       real*8 :: xin(m),xout(m)
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
       end subroutine dreorder_gpu


          
         

      subroutine lreorder1d_gpu(m,A, iperm, streamid_in )
      use cudafor
      integer :: m
      integer*8 :: A(m)
      integer :: iperm(m)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      integer*8 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: ierr


      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int(m/nthread)+1)
      call lreorder_gpu<<<nblock,nthread,0,streamid>>>(m,A,tmp,iperm)

      if (present(streamid_in)) then
        ierr = cudaMemcpyAsync(A, tmp, m, &
     & cudaMemcpyDeviceToDevice,streamid)
        call assert(ierr.eq.0,'lreorder1d_gpu:cudaMemcpyAsync',ierr)
      else
        ierr = cudaMemcpy(A, tmp, m )
        call assert(ierr.eq.0,'lreorder1d_gpu:cudaMemcpy',ierr)
      endif

      return
      end subroutine lreorder1d_gpu

       
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

      call ireorder_gpu<<<nblock,nthread,0,streamid>>>(m, &
     & A, tmp, iperm)

      if (present(streamid_in)) then
        ierr = cudaMemcpyAsync(A, tmp, m, &
     & cudaMemcpyDeviceToDevice, streamid )
        call assert(ierr.eq.0,'ireorder1d_gpu:cudaMemcpyAsync',ierr)
      else
        ierr = cudaMemcpy(A, tmp, m )
        call assert(ierr.eq.0,'ireorder1d_gpu:cudaMemcpy',ierr)
      endif

      return
      end subroutine ireorder1d_gpu

       
      subroutine sreorder1d_gpu(m,A, iperm, streamid_in )
      use cudafor
      integer :: m
      real*4 :: A(m)
      integer :: iperm(m)

      integer, intent(in),optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      real*4 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: ierr

      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int(m/nthread)+1)
      call sreorder_gpu<<<nblock,nthread,0,streamid>>>(m,A,tmp,iperm)

      if (present(streamid_in)) then
        ierr = cudaMemcpyAsync(A, tmp, m, &
     & cudaMemcpyDeviceToDevice,streamid )
        call assert(ierr.eq.0,'sreorder1d_gpu:cudaMemcpyAsync',ierr)
      else
        ierr = cudaMemcpy(A, tmp, m )
        call assert(ierr.eq.0,'sreorder1d_gpu:cudaMemcpy',ierr)
      endif

      return
      end subroutine sreorder1d_gpu

       
      subroutine dreorder1d_gpu(m,A, iperm, streamid_in )
      use cudafor
      integer :: m
      real*8 :: A(m)
      integer :: iperm(m)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      real*8 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: ierr


      streamid = 0
      if (present(streamid_in)) streamid = streamid_in


      nblock = min(32*1024,int( m/nthread) + 1 )
      call dreorder_gpu<<<nblock,nthread,0,streamid>>>(m,A,tmp,iperm)


      if (present(streamid_in)) then
        ierr = cudaMemcpyAsync(A, tmp, m, &
     & cudaMemcpyDeviceToDevice, streamid )
        call assert(ierr.eq.0,'dreorder1d_gpu:cudaMemcpyAsync',ierr)
      else
        ierr = cudaMemcpy(A, tmp, m )
        call assert(ierr.eq.0,'dreorder1d_gpu:cudaMemcpy',ierr)
      endif

      return
      end subroutine dreorder1d_gpu

       

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
        call lreorder_gpu<<<nblock,nthread,0,streamid>>>(m, &
     & A(1,j), tmp, iperm)

        if (present(streamid_in)) then
          ierr = cudaMemcpyAsync(A(1,j), tmp, m, &
     & cudaMemcpyDeviceToDevice,streamid )
          call assert(ierr.eq.0,'lreorder2d_gpu:cudaMemcpyAsync',ierr)
        else
          ierr = cudaMemcpy(A(1,j), tmp, m )
          call assert(ierr.eq.0,'lreorder2d_gpu:cudaMemcpy',ierr)
        endif
      enddo

      return
      end subroutine lreorder2d_gpu

       
      subroutine ireorder2d_gpu(m,n,A,lda, iperm, streamid_in )
      use cudafor
      implicit none

      integer,intent(in) :: m,n,lda
      integer, target :: A(lda,*)
      integer :: iperm(*)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      integer :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer:: nblock 
      integer :: j, ierr

      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min( 32*1024, int(m/nthread)+1)
      do j=1,n
        call ireorder_gpu<<<nblock,nthread,0,streamid>>>(m, &
     & A(1,j), tmp, iperm)

        if (present(streamid_in)) then
          ierr = cudaMemcpyAsync(A(1,j), tmp, m, &
     & cudaMemcpyDeviceToDevice,streamid )
          call assert(ierr.eq.0,'ireorder2d_gpu:cudaMemcpyAsync',ierr)
        else
          ierr = cudaMemcpy(A(1,j), tmp, m )
          call assert(ierr.eq.0,'ireorder2d_gpu:cudaMemcpy',ierr)
        endif
      enddo

      return
      end subroutine ireorder2d_gpu

       
      subroutine sreorder2d_gpu(m,n,A,lda, iperm, streamid_in )
      use cudafor
      implicit none

      integer,intent(in) :: m,n,lda
      real*4, target :: A(lda,*)
      integer :: iperm(*)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      real*4 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: j, ierr

      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int(m/nthread)+1)
      do j=1,n
        call sreorder_gpu<<<nblock,nthread,0,streamid>>>(m, &
     & A(1,j), tmp, iperm)

        if (present(streamid_in)) then
           ierr = cudaMemcpyAsync(A(1,j), tmp, m, &
     & cudaMemcpyDeviceToDevice,streamid )
           call assert(ierr.eq.0,'sreorder2d_gpu:cudaMemcpyAsync',ierr)
        else
           ierr = cudaMemcpy(A(1,j), tmp, m )
           call assert(ierr.eq.0,'sreorder2d_gpu:cudaMemcpy',ierr)
        endif
      enddo

      return
      end subroutine sreorder2d_gpu

       
      subroutine dreorder2d_gpu(m,n,A,lda, iperm, streamid_in )
      use cudafor
      implicit none

      integer,intent(in) :: m,n,lda
      real*8, target :: A(lda,*)
      integer :: iperm(*)

      integer, intent(in), optional :: streamid_in
      integer :: streamid

      attributes(device) :: A, iperm

      real*8 :: tmp(m)
      attributes(device) :: tmp

      integer,parameter :: nthread = 256
      integer :: nblock 
      integer :: j, ierr


      streamid = 0
      if (present(streamid_in)) streamid = streamid_in

      nblock = min(32*1024,int( m/nthread)+1 )
      do j=1,n
        call dreorder_gpu<<<nblock,nthread,0,streamid>>>(m, &
     & A(1,j), tmp, iperm)

        if (present(streamid_in)) then
          ierr = cudaMemcpyAsync(A(1,j), tmp, m, &
     & cudaMemcpyDeviceToDevice, streamid)
          call assert(ierr.eq.0,'dreorder2d_gpu:cudaMemcpyAsync',ierr)

        else
          ierr = cudaMemcpy(A(1,j), tmp, m)
          call assert(ierr.eq.0,'dreorder2d_gpu:cudaMemcpy',ierr)
        endif
      enddo

      return
      end subroutine dreorder2d_gpu

       

       end module reorder_gpu_mod
      module gen_perm_gpu_mod
      use cudafor
      use precision_mod_gpu
      implicit none

      contains

       attributes(global) &
       subroutine isetval_kernel_gpu( n, v, dvalue)
       use precision_mod_gpu
       use cudafor
       implicit none

       integer, value :: n
       integer(kind=4), value :: dvalue
       integer(kind=4) :: v(*)


       integer :: i, ith, nthreads

       ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) + &
               ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
                     (blockDim%x * blockDim%y)

       nthreads = (blockDim%x * blockDim%y)*(gridDim%x * gridDim%y)
       do i=ith,n,nthreads
         v(i) = dvalue
       enddo
       return
       end subroutine isetval_kernel_gpu


       attributes(host) &
       subroutine isetval_gpu( n, v, dvalue, streamid_in )
       use precision_mod_gpu
       use cudafor
       implicit none

       integer, intent(in) :: n
       integer(kind=4), intent(in) :: dvalue
       integer(kind=4), intent(inout) :: v(*)

       integer, intent(in), optional :: streamid_in
       integer :: streamid

       attributes(device) :: v

       integer, parameter :: nthread = 256

       type(dim3) :: tgrid, tblock

       tblock%x = nthread
       tblock%y = 1
       tblock%z = 1

       tgrid%x = min(32*1024,int(n/nthread)+1)
       tgrid%y = 1
       tgrid%z = 1

       streamid = 0
       if (present(streamid_in)) streamid = streamid_in

       call isetval_kernel_gpu<<<tgrid,tblock,0,streamid>>>(n,v,dvalue)

       return
       end subroutine isetval_gpu

       attributes(global) &
       subroutine isum_col_gpu( m,n, lda, lcount, icount )
       use cudafor
       implicit none 
       integer, value :: m,n, lda
       integer :: lcount(lda,n)
       integer :: icount(m)

       attributes(device) :: icount, lcount

       integer :: i,j, icount_local
       integer :: gith,ith, iblock, nthread, gnthread

       ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
       iblock = blockIdx%x + (blockIdx%y - 1)*gridDim%x
       nthread = blockDim%x * blockDim%y
       gnthread = nthread * ( gridDim%x * gridDim%y )
       gith = ith + (iblock-1)*nthread

       do i=gith,m,gnthread
           icount_local = 0
           do j=1,n
              icount_local = icount_local + lcount(i, j)
           enddo
           icount(i) = icount_local
       enddo
       return
       end subroutine isum_col_gpu

      attributes(global) &
      subroutine iprefix_sum_gpu_kernel(ix,n)
      use cudafor
      integer, value :: n
      integer :: ix(*)

      attributes(device) :: ix

! --------------------------------------------
! compute prefix sum with starting value ival
! [ival, ix(1)+ival, ix(1)+ix(2) + ival, + ... + sum( ix(:) ) + ival]
! --------------------------------------------
      integer :: j, ix_j, cumulative_sum


! ---------------------------------
! use inefficient algorithm for now
! ---------------------------------
      if ((threadIdx%x == 1).and.(blockIdx%x == 1) .and. &
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


      attributes(global) &
      subroutine ishift_and_add_gpu_kernel(x_in, x_out,n,ival)
      use cudafor
      implicit none
      integer, value :: n, ival
      integer :: x_in(n), x_out(n+1)

      attributes(device) :: x_in, x_out

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
      implicit none
      integer, value :: n
      integer :: ix(n+1)
      integer, value :: ival

      integer :: ix_tmp(n+1)
      attributes(device) :: ix, ix_tmp

      ix_tmp = ix

! print*,'before prefix sum, ix_tmp ', ix_tmp


! --------------------
! very slow and simple
! --------------------
      call iprefix_sum_gpu_kernel<<<1,1>>>(ix_tmp, n)



! print*,'after prefix sum, ix_tmp ', ix_tmp

! ------------------------------------------------
! shift and add
! converrt [x1, x2, ... xn] to [ival, x1+ival, ... xn+ival]
! ------------------------------------------------

      call ishift_and_add_gpu_kernel<<<32,192>>>(ix_tmp,ix, n, ival)
      return
      end subroutine iprefix_sum_gpu
      attributes(global) &
      subroutine setup_lcount_gpu(nx,ny,ncopy,lcount,xstart)
      use cudafor
      implicit none
      integer, value :: nx,ny,ncopy
      integer(kind=4) :: lcount(nx*ny, ncopy)
      integer :: xstart(nx*ny+1)

      attributes(device) :: lcount, xstart

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

       attributes(global) &
       subroutine gen_perm_gpu_pass1( nx,ny,xmin,ymin, &
     & inv_dx,inv_dy, n,gid,xy,xydim, ncopy,lcount)
       use cudafor
       use precision_mod_gpu
       integer, value :: nx,ny,ncopy
       real(kind=work_p), value :: xmin, ymin, inv_dx, inv_dy

       integer, value :: n, xydim
       integer*8 :: gid(n)
       real(kind=work_p) :: xy(xydim,2)
       integer(kind=4) :: lcount( nx*ny, ncopy)

       attributes(device) :: gid,xy,lcount

       integer :: ith, nthread, iblock, nblock
       integer :: gith, gnthread
       integer :: i, k, ix,iy , idummy,icopy


! ------------------------------------
! assume lcount(nx*ny,:) to be zero
! ------------------------------------

       ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
       nthread = blockDim%x * blockDim%y
       iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x

       nblock = gridDim%x * gridDim%y
       gnthread = nblock * nthread
       gith = ith + (iblock-1)*nthread

       icopy = 1 + mod(iblock-1,ncopy)

! -------------------------
! perform geometric hashing
! -------------------------

       do i=gith,n,gnthread
        if (gid(i) <= 0) then
         ix = nx
         iy = ny
        else
         ix = 1+int( (xy(i,1) - xmin) * inv_dx )
         iy = 1+int( (xy(i,2) - ymin) * inv_dy )

         if ((iy < 0).or.(iy > ny)) then
! print*,'i,iy,xy(i,2) ',i,iy,xy(i,2)
           iy = max(1,min(ny,iy))
         endif

         if ((ix < 0).or.(ix > ny)) then
! print*,'i,ix,xy(i,1) ',i,ix,xy(i,1)
           ix = max(1,min(ny,ix))
         endif

         ix = max(1,min(nx,ix))
         iy = max(1,min(ny,iy))

        endif


         k = ix + (iy-1)*nx
         idummy = atomicadd( lcount(k,icopy), 1 )
       enddo

       return
       end subroutine gen_perm_gpu_pass1
       attributes(global) &
       subroutine gen_perm_gpu_pass2(nx,ny,xmin,ymin, &
     & inv_dx,inv_dy, n,gid,xy,xydim,ncopy,lcount,iperm,xstart)
       use cudafor
       use precision_mod_gpu
       implicit none

       integer, value :: nx,ny,n,ncopy
       real(kind=work_p),value :: xmin, ymin, inv_dx, inv_dy
       integer*8 :: gid(n)
       integer, value :: xydim
       real(kind=work_p) :: xy(xydim,2)
       integer(kind=4) :: lcount(nx*ny, ncopy)
       integer :: iperm(n)
       integer :: xstart(nx*ny+1)

       attributes(device) :: gid, xy, lcount, xstart, iperm

       integer :: ith,nthread,iblock, icopy
       integer :: nblock,gnthread,gith
       integer :: i,ix,iy,k,ip

       ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
       nthread = blockDim%x * blockDim%y
       iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x

       nblock = gridDim%x * gridDim%y
       gnthread = nblock * nthread
       gith = ith + (iblock-1)*nthread

       icopy = 1 + mod(iblock-1,ncopy)

! -------------------------
! perform geometric hashing
! -------------------------
       do i=gith,n,gnthread
        if (gid(i) <= 0) then
         ix = nx
         iy = ny
        else
         ix = 1+int( (xy(i,1) - xmin) * inv_dx )
         iy = 1+int( (xy(i,2) - ymin) * inv_dy )


         if ((iy < 0).or.(iy > ny)) then
! print*,'i,iy,xy(i,2) ',i,iy,xy(i,2)
           iy = max(1,min(ny,iy))
         endif

         if ((ix < 0).or.(ix > ny)) then
! print*,'i,ix,xy(i,1) ',i,ix,xy(i,1)
           ix = max(1,min(ny,ix))
         endif

         ix = max(1, min(nx, ix) )
         iy = max(1, min(ny, iy) )
        endif

         k = ix + (iy-1)*nx

         ip = atomicadd( lcount(k,icopy), 1 )
         iperm(ip) = i
       enddo

       return
       end subroutine gen_perm_gpu_pass2
       subroutine gen_perm_gpu( nx, ny, xmin,ymin, inv_dx,inv_dy, &
     & n, gid, xy, xydim, iperm, xstart )
       use cudafor
       use sml_module
       use precision_mod_gpu
       implicit none

       integer, intent(in) :: nx, ny
       real(kind=work_p), intent(in) :: xmin,ymin, inv_dx,inv_dy
       integer, intent(in) :: n, xydim
       integer*8, intent(in) :: gid(n)
       real(kind=work_p) :: xy(xydim, 2)
       integer, intent(inout) :: iperm(n)
       integer, intent(inout) :: xstart(nx*ny+1)


       integer, parameter :: nblock = 32
       integer, parameter :: ncopy = nblock
       integer, parameter :: nthread = 32*6

       integer(kind=4), allocatable, dimension(:,:) :: lcount

       attributes(device) :: gid, xy, iperm, xstart
       attributes(device) :: lcount

       integer, parameter :: block_limit = 2**15-1

       integer :: mm,nn,lda
       integer :: istat, istart_val
       
       integer, parameter :: idebug = 0
       integer(kind=4), allocatable, dimension(:,:) :: h_lcount
       integer :: i,j,k,iblock
 

! ---------------------------------------------
! 1st pass to count number of particles in cell
! ---------------------------------------------
       allocate( lcount( nx*ny,ncopy ) , stat=istat)
       if (istat.ne.0) then
         write(*,*) sml_mype, 'gen_perm_gpu: allocate lcount return istat=',istat
         write(*,*) 'nx, ny, nblock ', nx, ny, nblock
         stop '** error ** '
       endif

       call isetval_gpu( size(lcount), lcount, 0)


       call gen_perm_gpu_pass1<<<nblock,nthread>>>( &
     & nx,ny, xmin,ymin, &
     & inv_dx,inv_dy, n,gid,xy,xydim, ncopy, lcount)

       if (idebug >= 2) then
          allocate( h_lcount(nx*ny,ncopy) )
          h_lcount = lcount

! print*,'after gen_perm_gpu_pass1 '
          do iblock=1,ncopy
            do j=1,ny
            do i=1,nx
              k = i + (j-1)*nx
! print 9010,i,j,iblock,h_lcount(k,iblock)
 9010 format(1x,'h_lcount(',i8,',',i8,',',i8,') = ',i9 )
            enddo
            enddo
          enddo
! print*,'sum(h_lcount) ', sum(h_lcount)
          deallocate( h_lcount )
         endif

! ------------------------
! setup pointers in xstart
! ------------------------

! ---------------------------------------
! compute xstart(k) = sum( lcount(k,:) )
! ---------------------------------------

       mm = size(lcount,1)
       nn = size(lcount,2)
       lda = size(lcount,1)

       call isum_col_gpu<<<nblock,nthread>>>(mm,nn,lda, lcount, xstart)
! --------------------------------------------
! compute xstart(:) as prefix sum on icount(:)
! --------------------------------------------
       istart_val = 1
       mm = size(xstart)-1
       call iprefix_sum_gpu( mm, xstart, istart_val)
     

! ------------------------------------------
! setup lcount(:,:) as pointers into iperm
! ------------------------------------------
       call setup_lcount_gpu<<<nblock,nthread>>>(nx,ny, &
     & ncopy,lcount,xstart)

! --------------------------
! 2nd pass to populate iperm
! --------------------------

       call gen_perm_gpu_pass2<<<nblock,nthread>>>( nx,ny,xmin,ymin, &
     & inv_dx,inv_dy,n,gid,xy,xydim,ncopy,lcount,iperm,xstart )


       deallocate( lcount, stat=istat )
       if (istat.ne.0) then
          write(*,*) sml_mype, 'gen_perm_gpu: dealloc(lcount),istat=',istat
          stop '** error ** '
       endif
       return
       end subroutine gen_perm_gpu

       end module gen_perm_gpu_mod



module dimensions_mod_gpu
! integer, parameter :: ptl_maxnum_dim = 22*100*100 + 5001
! integer, parameter :: ptl_nummax_dim = ptl_maxnum_dim
! integer, parameter :: maxnum_dim = ptl_maxnum_dim
! integer, parameter :: nummax_dim = maxnum_dim
! integer, parameter :: lost_maxnum_dim = 100 + (ptl_maxnum_dim / 10)
! integer, parameter :: lost_nummax_dim = lost_maxnum_dim

! integer, parameter :: nnode_dim = 20694
! integer, parameter :: ntriangle_dim = 41087

! integer, parameter :: guess_table_n1_dim = 512
! integer, parameter :: guess_table_n2_dim = 512
! integer, parameter :: guess_list_size = 7005393

 ! -----------------------------------
 ! nseg_dim used in boundary_class_gpu
 ! -----------------------------------
 integer, parameter :: iseg_dim = 100

 ! ------------------------------
 ! nrho_dim used in psn_class_gpu
 ! ------------------------------
! integer, parameter :: nrho_dim = 8
! integer, parameter :: nphi_dim = 32

 ! -------------------------------------
 ! nthreads_dim used in diag_module_gpu
 ! -------------------------------------
 integer, parameter :: grid_dim = 384
 integer, parameter :: block_dim = 64

 integer, parameter :: nthreads_dim = grid_dim * block_dim
end module dimensions_mod_gpu

module sml_module_gpu
use precision_mod_gpu
use sml_module, only : &
   sml_e_charge_host => sml_e_charge, &
   sml_epsilon0_host => sml_epsilon0, &
   sml_prot_mass_host => sml_prot_mass, &
   sml_elec_mass_host => sml_elec_mass, &
   sml_n_vf_diag_host => sml_n_vf_diag, &
   sml_nlarmor_host => sml_nlarmor, &
   sml_nrk_host => sml_nrk, &
   sml_boundary_diagonal_host => sml_boundary_diagonal, &
   sml_j2ev_host => sml_j2ev, &
   sml_ev2j_host => sml_ev2j

implicit none
  !! delta-f weight evolution switch. 
  !! false for whole (grad_b, curl_b,and ExB), true for exb only
  logical :: sml_dwdt_exb_only 

  !! delta-f weight evolution switch. 
  !! false for (1-w) factor true for (1) factor. default is .false.
  logical :: sml_dwdt_fix_bg 

  logical :: sml_turb_efield ! set zero turbulece efield
  logical :: sml_00_efield ! set zero 00 efield

  !! delta-f switch. 0 for off, 1 for on. delta-f simulation 
  !! is not verified yet, especially for output file.
  logical :: sml_deltaf
  logical :: sml_electron_on 
  logical :: sml_extra_dwdt ! extra dwdt 
  integer :: sml_deltaf_f0_mode

  !! Bounce routine switch 0 for off, 
  !! 1 for inner boundary, 
  !! 2 for both boundaries
  integer :: sml_bounce 
  integer :: sml_bounce_zero_weight

  !! Neutral routine switch .f. for off
  logical :: sml_neutral

  !! simulation boundary in R-Z space
  real (kind=work_p) :: sml_2pi
  real (kind=work_p) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z 



  !! for sml_deltaf_f0_mode==-1 1/Ln , 1/Lt - artificial f0
  real (kind=work_p) :: sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e 



  real(kind=work_p), parameter :: sml_e_charge = sml_e_charge_host
  real(kind=work_p), parameter :: sml_epsilon0 = sml_epsilon0_host
  real(kind=work_p), parameter :: sml_prot_mass = sml_prot_mass_host
  real(kind=work_p), parameter :: sml_elec_mass = sml_elec_mass_host

  real(kind=work_p), parameter :: sml_j2ev = sml_j2ev_host
  real(kind=work_p), parameter :: sml_ev2j = sml_ev2j_host

  real(kind=work_p), parameter :: sml_boundary_diagonal = sml_boundary_diagonal_host
  integer, parameter :: sml_n_vf_diag = sml_n_vf_diag_host
  integer, parameter :: sml_nlarmor = sml_nlarmor_host
  integer, parameter :: sml_nrk = sml_nrk_host

  !! Inner and outer boundary for initial loading
  real(kind=work_p) :: sml_inpsi, sml_outpsi 

  integer :: sml_comm !! communicator
  integer :: sml_comm_null !! MPI_COMM_NULL for adios posix
  integer :: sml_totalpe !! total number of processors
  integer :: sml_mype !! process index
  integer :: sml_nthreads
  integer :: sml_gstep
  integer :: sml_ipc, ipc_gpu
  integer :: sml_sheath_mode
 
  real (kind=work_p) :: sml_time !! simulation time 
  !! -1 for inverted B 1 for normal B <-- given by sml_invert_B
  real (kind=work_p) :: sml_bp_sign 
  !! -bp_sign for co-current bt, 
  !! bp_sign for couter-current bt <- given by sml_co_curr_bt
  real (kind=work_p) :: sml_bt_sign 

  logical :: sml_ignore_drift_near_wall
  real (kind=work_p):: sml_ignore_drift_r0
  real (kind=work_p):: sml_ignore_drift_z0
  real (kind=work_p):: sml_ignore_drift_slope1
  real (kind=work_p):: sml_ignore_drift_slope2


  attributes(device) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z 
  attributes(device) :: sml_bounce, sml_bounce_zero_weight
  attributes(device) :: sml_2pi, ipc_gpu, epc_gpu
  attributes(device) :: sml_turb_efield, sml_00_efield
  attributes(device) :: sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e
  attributes(device) :: sml_deltaf_f0_mode, sml_dwdt_fix_bg
  attributes(device) :: sml_dwdt_exb_only, sml_deltaf, sml_electron_on, sml_extra_dwdt, sml_neutral
  attributes(device) :: sml_inpsi, sml_outpsi
  attributes(device) :: sml_comm, sml_comm_null 
  attributes(device) :: sml_totalpe,sml_mype, sml_nthreads, sml_gstep, sml_sheath_mode, sml_ipc
  attributes(device) :: sml_time, sml_bp_sign, sml_bt_sign
  attributes(device) :: sml_ignore_drift_near_wall, sml_ignore_drift_r0, sml_ignore_drift_z0, &
                        sml_ignore_drift_slope1, sml_ignore_drift_slope2
  contains

  attributes(host) &
  subroutine update_device_sml()
  use sml_module, only : &
     sml_bd_min_r_host => sml_bd_min_r, &
     sml_bd_max_r_host => sml_bd_max_r, &
     sml_bd_min_z_host => sml_bd_min_z, &
     sml_bd_max_z_host => sml_bd_max_z, &
     sml_bounce_host => sml_bounce, &
     sml_bounce_zero_weight_host => sml_bounce_zero_weight, &
     sml_2pi_host => sml_2pi, &
     sml_00_efield_host => sml_00_efield, &
     sml_turb_efield_host => sml_turb_efield, &
     sml_f0_1_Ln_host => sml_f0_1_Ln, &
     sml_f0_1_Lt_host => sml_f0_1_Lt, &
     sml_f0_1_Lt_e_host => sml_f0_1_Lt_e, &
     sml_deltaf_f0_mode_host => sml_deltaf_f0_mode, &
     sml_deltaf_host => sml_deltaf, &
     sml_electron_on_host => sml_electron_on, &
     sml_extra_dwdt_host => sml_extra_dwdt, &
     sml_neutral_host => sml_neutral, &
     sml_dwdt_exb_only_host => sml_dwdt_exb_only, &
     sml_dwdt_fix_bg_host => sml_dwdt_fix_bg, &
     sml_comm_host => sml_comm, &
     sml_comm_null_host => sml_comm_null, &
     sml_totalpe_host => sml_totalpe, &
     sml_mype_host => sml_mype, &
     sml_nthreads_host => sml_nthreads, &
     sml_gstep_host => sml_gstep, &
     sml_ipc_host => sml_ipc, &
     sml_sheath_mode_host => sml_sheath_mode, &
     sml_time_host => sml_time, &
     sml_inpsi_host => sml_inpsi, &
     sml_outpsi_host => sml_outpsi, &
     sml_bp_sign_host => sml_bp_sign, &
     sml_bt_sign_host => sml_bt_sign, &
     sml_ignore_drift_near_wall_host => sml_ignore_drift_near_wall, &
     sml_ignore_drift_r0_host => sml_ignore_drift_r0, &
     sml_ignore_drift_z0_host => sml_ignore_drift_z0, &
     sml_ignore_drift_slope1_host => sml_ignore_drift_slope1, &
     sml_ignore_drift_slope2_host => sml_ignore_drift_slope2

  implicit none 

  sml_bd_min_r = sml_bd_min_r_host
  sml_bd_max_r = sml_bd_max_r_host
  sml_bd_min_z = sml_bd_min_z_host
  sml_bd_max_z = sml_bd_max_z_host

  sml_bounce = sml_bounce_host
  sml_bounce_zero_weight = sml_bounce_zero_weight_host

  sml_2pi = sml_2pi_host

  sml_00_efield = sml_00_efield_host
  sml_turb_efield = sml_turb_efield_host

  sml_f0_1_Ln = sml_f0_1_Ln_host
  sml_f0_1_Lt = sml_f0_1_Lt_host
  sml_f0_1_Lt_e = sml_f0_1_Lt_e_host

  sml_deltaf_f0_mode = sml_deltaf_f0_mode_host
  sml_deltaf = sml_deltaf_host
  sml_electron_on = sml_electron_on_host
  sml_extra_dwdt = sml_extra_dwdt_host
  sml_neutral = sml_neutral_host
  sml_dwdt_exb_only = sml_dwdt_exb_only_host
  sml_dwdt_fix_bg = sml_dwdt_fix_bg_host

  sml_comm = sml_comm_host
  sml_comm_null = sml_comm_null_host
  sml_totalpe = sml_totalpe_host
  sml_mype = sml_mype_host
  sml_nthreads = sml_nthreads_host
  sml_gstep = sml_gstep_host
  sml_ipc = sml_ipc_host
  ipc_gpu=sml_ipc_host
  sml_sheath_mode = sml_sheath_mode_host
  sml_inpsi = sml_inpsi_host
  sml_outpsi = sml_outpsi_host

  sml_time = sml_time_host
  sml_bp_sign = sml_bp_sign_host
  sml_bt_sign = sml_bt_sign_host

  sml_ignore_drift_near_wall = sml_ignore_drift_near_wall_host
  sml_ignore_drift_r0 = sml_ignore_drift_r0_host
  sml_ignore_drift_z0 = sml_ignore_drift_z0_host
  sml_ignore_drift_slope1 = sml_ignore_drift_slope1_host
  sml_ignore_drift_slope2 = sml_ignore_drift_slope2_host
  if(sml_extra_dwdt_host) print *, "Warning: GPU implementation is not correct for this case"
  return
  end subroutine update_device_sml



  attributes(host) &
  subroutine update_host_sml()
  use sml_module, only : &
     sml_bd_min_r_host => sml_bd_min_r, &
     sml_bd_max_r_host => sml_bd_max_r, &
     sml_bd_min_z_host => sml_bd_min_z, &
     sml_bd_max_z_host => sml_bd_max_z, &
     sml_bounce_host => sml_bounce, &
     sml_bounce_zero_weight_host => sml_bounce_zero_weight, &
     sml_2pi_host => sml_2pi, &
     sml_00_efield_host => sml_00_efield, &
     sml_turb_efield_host => sml_turb_efield, &
     sml_f0_1_Ln_host => sml_f0_1_Ln, &
     sml_f0_1_Lt_host => sml_f0_1_Lt, &
     sml_f0_1_Lt_e_host => sml_f0_1_Lt_e, &
     sml_deltaf_f0_mode_host => sml_deltaf_f0_mode, &
     sml_deltaf_host => sml_deltaf, &
     sml_electron_on_host => sml_electron_on, &
     sml_dwdt_exb_only_host => sml_dwdt_exb_only, &
     sml_dwdt_fix_bg_host => sml_dwdt_fix_bg, &
     sml_comm_host => sml_comm, &
     sml_comm_null_host => sml_comm_null, &
     sml_totalpe_host => sml_totalpe, &
     sml_mype_host => sml_mype, &
     sml_time_host => sml_time, &
     sml_inpsi_host => sml_inpsi, &
     sml_outpsi_host => sml_outpsi, &
     sml_bp_sign_host => sml_bp_sign, &
     sml_bt_sign_host => sml_bt_sign, &
     sml_ipc_host => sml_ipc, &
     sml_ignore_drift_near_wall_host => sml_ignore_drift_near_wall, &
     sml_ignore_drift_r0_host => sml_ignore_drift_r0, &
     sml_ignore_drift_z0_host => sml_ignore_drift_z0, &
     sml_ignore_drift_slope1_host => sml_ignore_drift_slope1, &
     sml_ignore_drift_slope2_host => sml_ignore_drift_slope2


  implicit none

  sml_bd_min_r_host = sml_bd_min_r
  sml_bd_max_r_host = sml_bd_max_r
  sml_bd_min_z_host = sml_bd_min_z
  sml_bd_max_z_host = sml_bd_max_z

  sml_bounce_host = sml_bounce
  sml_bounce_zero_weight_host = sml_bounce_zero_weight

  sml_2pi_host = sml_2pi

  sml_00_efield_host = sml_00_efield
  sml_turb_efield_host = sml_turb_efield

  sml_f0_1_Ln_host = sml_f0_1_Ln
  sml_f0_1_Lt_host = sml_f0_1_Lt
  sml_f0_1_Lt_e_host = sml_f0_1_Lt_e

  sml_deltaf_f0_mode_host = sml_deltaf_f0_mode
  sml_deltaf_host = sml_deltaf
  sml_electron_on_host = sml_electron_on
  sml_dwdt_exb_only_host = sml_dwdt_exb_only
  sml_dwdt_fix_bg_host = sml_dwdt_fix_bg

  sml_comm_host = sml_comm
  sml_comm_null_host = sml_comm_null
  sml_totalpe_host = sml_totalpe
  sml_mype_host = sml_mype

  sml_inpsi_host = sml_inpsi
  sml_outpsi_host = sml_outpsi

  sml_time_host = sml_time
  sml_bp_sign_host = sml_bp_sign
  sml_bt_sign_host = sml_bt_sign

  sml_ignore_drift_near_wall_host = sml_ignore_drift_near_wall
  sml_ignore_drift_r0_host = sml_ignore_drift_r0
  sml_ignore_drift_z0_host = sml_ignore_drift_z0
  sml_ignore_drift_slope1_host = sml_ignore_drift_slope1
  sml_ignore_drift_slope2_host = sml_ignore_drift_slope2

  return
  end subroutine update_host_sml

end module sml_module_gpu
module neu_module_gpu
use neu_module, only: neu_weight_sum_lost
use precision_mod_gpu
use cudafor
implicit none

  integer, parameter :: neu_weight_sum_lost_gpu_dim = 256
  real (kind=work_p) :: neu_weight_sum_lost_gpu(neu_weight_sum_lost_gpu_dim)
  attributes(device) :: neu_weight_sum_lost_gpu

contains

  attributes(host) &
  subroutine update_device_neu()
  implicit none

  integer :: ierr,icount


  neu_weight_sum_lost_gpu = 0

  return
  end subroutine update_device_neu 

  attributes(host) &
  subroutine update_host_neu()
  implicit none

  real(kind=work_p) :: tmp_neu_weight_sum_lost(neu_weight_sum_lost_gpu_dim)
  integer :: lb


  if (allocated(neu_weight_sum_lost)) then
   if (size(neu_weight_sum_lost) >= 1) then
    tmp_neu_weight_sum_lost = neu_weight_sum_lost_gpu
    lb = lbound(neu_weight_sum_lost,1)
    neu_weight_sum_lost(lb) = neu_weight_sum_lost(lb) + &
                            sum(tmp_neu_weight_sum_lost)
   endif
  endif

  return
  end subroutine update_host_neu
end module neu_module_gpu


module eq_module_gpu
use precision_mod_gpu
real (kind=work_p) :: eq_min_r,eq_max_r,eq_min_z,eq_max_z 
real (kind=work_p) :: eq_x_psi, eq_x_z, eq_axis_r, eq_axis_z, eq_x_slope, eq_x_r
real (kind=work_p) :: eq_x2_r, eq_x2_z, eq_x2_slope
attributes(device) :: eq_min_r,eq_max_r,eq_min_z,eq_max_z 
attributes(device) :: eq_x_psi, eq_x_z, eq_axis_r, eq_axis_z, eq_x_slope, eq_x_r
attributes(device) :: eq_x2_r, eq_x2_z, eq_x2_slope


  ! data structure for equilbirum profile
  type eq_ftn_type
     integer :: shape
     real (kind=work_p):: inx(3), iny(3)
     real (kind=work_p):: sv(6)
     ! for arbitrary profile function - use pspline
     ! character (len=100) :: filename
     ! type (EZspline1_r8) :: spl
     real (kind=work_p) :: min, max
  end type eq_ftn_type

type(eq_ftn_type) :: eq_tempi, eq_tempe, eq_den, eq_flowi, eq_flowe

attributes(device) :: eq_tempi, eq_tempe, eq_den, eq_flowi, eq_flowe

contains



  attributes(host) &
  subroutine copy_eq_ftn_type(eq_tempi_gpu, eq_tempi_host)
  ! --------------------------------------------------
  ! copy structure type(eq_ftn_type) from host to gpu
  ! --------------------------------------------------
  use eq_module, only: eq_ftn_type_host => eq_ftn_type
  implicit none

  type(eq_ftn_type), intent(inout) :: eq_tempi_gpu
  type(eq_ftn_type_host),intent(in) :: eq_tempi_host

  attributes(device) :: eq_tempi_gpu

  type(eq_ftn_type) :: p_eq_tempi


  p_eq_tempi%inx = eq_tempi_host%inx
  p_eq_tempi%iny = eq_tempi_host%iny
  p_eq_tempi%sv = eq_tempi_host%sv
  p_eq_tempi%shape = eq_tempi_host%shape
  p_eq_tempi%min = eq_tempi_host%min
  p_eq_tempi%max = eq_tempi_host%max

  eq_tempi_gpu = p_eq_tempi


  return
  end subroutine copy_eq_ftn_type

  attributes(host) &
  subroutine update_device_eq()
  use sml_module, only : sml_mype

  use eq_module, only : &
  eq_ftn_type_host => eq_ftn_type, &
  eq_axis_r_host => eq_axis_r, &
  eq_axis_z_host => eq_axis_z, &
  eq_min_r_host => eq_min_r, &
  eq_max_r_host => eq_max_r, &
  eq_min_z_host => eq_min_z, &
  eq_max_z_host => eq_max_z, &
  eq_x_psi_host => eq_x_psi, &
  eq_x_r_host => eq_x_r, &
  eq_x_z_host => eq_x_z, &
  eq_tempi_host => eq_tempi, &
  eq_tempe_host => eq_tempe, &
  eq_den_host => eq_den, &
  eq_flowi_host => eq_flowi, &
  eq_flowe_host => eq_flowe, &
  eq_x2_z_host => eq_x2_z, &
  eq_x2_r_host => eq_x2_r, &
  eq_x2_slope_host => eq_x2_slope, &
  eq_x_slope_host => eq_x_slope

  implicit none

  integer, parameter :: idebug = 0

  eq_axis_r = eq_axis_r_host
  eq_min_r = eq_min_r_host
  eq_max_r = eq_max_r_host
  eq_min_z = eq_min_z_host
  eq_max_z = eq_max_z_host
  eq_x_psi = eq_x_psi_host
  eq_x_r = eq_x_r_host
  eq_x_z = eq_x_z_host
  eq_x2_r = eq_x2_r_host
  eq_x2_z = eq_x2_z_host
  eq_x2_slope = eq_x2_slope_host
  eq_x_slope = eq_x_slope_host

  if (idebug >= 1) then
! if (sml_mype == 0) print*,'before copy_eq_ftn_type()'
  endif



  call copy_eq_ftn_type(eq_tempi, eq_tempi_host)
  call copy_eq_ftn_type(eq_tempe, eq_tempe_host)
  call copy_eq_ftn_type(eq_den, eq_den_host)
  call copy_eq_ftn_type(eq_flowi, eq_flowi_host)
  call copy_eq_ftn_type(eq_flowe, eq_flowe_host)

  end subroutine update_device_eq


end module eq_module_gpu
module itp_module_gpu
use precision_mod_gpu
implicit none

real(kind=work_p) :: itp_min_psi,itp_max_psi
attributes(device) :: itp_min_psi, itp_max_psi

contains
  attributes(host) &
  subroutine update_device_itp()
  use itp_module, only : &
    itp_min_psi_host => itp_min_psi, &
    itp_max_psi_host => itp_max_psi

  itp_min_psi = itp_min_psi_host
  itp_max_psi = itp_max_psi_host

  return
  end subroutine update_device_itp

end module itp_module_gpu
module one_d_cub_mod_gpu
use precision_mod_gpu
use cudafor

implicit none
   integer, parameter :: ndeg = 3
   integer, parameter :: max_npsi = eq_mpsi-1

   real (kind=work_p),allocatable :: one_d_cub_psi_acoef_gpu(:,:)

   real (kind=work_p) :: one_d_cub_dpsi_inv_gpu, one_d_cub_psimin_gpu

   attributes(device) :: one_d_cub_dpsi_inv_gpu, one_d_cub_psimin_gpu
   attributes(device) :: one_d_cub_psi_acoef_gpu

contains

   attributes(host) &
   subroutine update_device_one_d()
   use assert_mod
   use precision_mod_gpu
   use one_d_cub_mod, only : &
       one_d_cub_psi_acoef, one_d_cub_dpsi_inv, one_d_cub_psimin
   use cudafor

implicit none

   integer :: npsi
   integer :: icount_one_d_acoef
   integer :: ierr, i, j
   integer :: lb1,ub1, lb2,ub2
   logical :: isvalid
   integer, parameter :: idebug = 0
   logical, parameter :: use_cudaMemcpy = .true.



   if (.not.allocated(one_d_cub_psi_acoef_gpu)) then
     lb1 = lbound(one_d_cub_psi_acoef,1)
     ub1 = ubound(one_d_cub_psi_acoef,1)
     lb2 = lbound(one_d_cub_psi_acoef,2)
     ub2 = ubound(one_d_cub_psi_acoef,2)
     allocate( one_d_cub_psi_acoef_gpu( lb1:ub1, lb2:ub2 ),stat=ierr )
     call assert(ierr.eq.0,'alloc(one_d_cub_psi_acoef_gpu)',ierr)

     one_d_cub_psi_acoef_gpu = 0
   endif


   isvalid = size(one_d_cub_psi_acoef).eq.size(one_d_cub_psi_acoef_gpu)
   call assert(isvalid, &
     'invalid size(one_d_cub_psi_acoef_gpu', &
     size(one_d_cub_psi_acoef_gpu))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!

   
   one_d_cub_dpsi_inv_gpu = one_d_cub_dpsi_inv
   one_d_cub_psimin_gpu = one_d_cub_psimin

   if(use_cudaMemcpy) then
        icount_one_d_acoef = size( one_d_cub_psi_acoef )
        ierr = cudaMemcpy( one_d_cub_psi_acoef_gpu, &
                           one_d_cub_psi_acoef, &
                           icount_one_d_acoef, &
                           cudaMemcpyHostToDevice )

        call assert(ierr.eq.0, &
               'cudaMemcpy(one_d_cub_psi_acoef_gpu',ierr)
   else

     one_d_cub_psi_acoef_gpu = one_d_cub_psi_acoef

   endif
   
   
   return
   end subroutine update_device_one_d



end module one_d_cub_mod_gpu
module bicub_mod_gpu
  use precision_mod_gpu
  use cudafor
  implicit none

  integer, parameter :: ndeg = 3
  integer, parameter :: max_nr = eq_mr-1
  integer, parameter :: max_nz = eq_mz-1

  integer :: nr_gpu, nz_gpu

  real (kind=work_p) :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu

  real (kind=work_p),allocatable :: rc_cub_gpu(:)
  real (kind=work_p),allocatable :: zc_cub_gpu(:)
  real (kind=work_p),allocatable :: acoef_cub_gpu(:,:,:,:)





! real (kind=work_p), dimension(max_nr) :: rc_cub_gpu
! real (kind=work_p), dimension(max_nz) :: zc_cub_gpu
! real (kind=work_p), dimension(0:ndeg,0:ndeg,max_nr,max_nz) :: acoef_cub_gpu

  attributes(device) :: rc_cub_gpu, zc_cub_gpu, acoef_cub_gpu
  attributes(device) :: nr_gpu, nz_gpu
  attributes(device) :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu

  public :: ndeg
  public :: update_device_bicub
  public :: rc_cub_gpu, zc_cub_gpu, nr_gpu, nz_gpu
  public :: acoef_cub_gpu
  public :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu


contains

  attributes(host) &
  subroutine update_device_bicub()
  use sml_module, only : sml_mype
  use precision_mod_gpu
  use bicub_mod, only : psi_bicub

  use cudafor
  
implicit none

! type(bicub_type), intent(in) :: bicub
  integer :: icount
  integer :: ierr 
  integer :: lb1,ub1, lb2,ub2, lb3,ub3, lb4,ub4

  integer, parameter :: idebug = 0


  logical, parameter :: use_cudaMemcpy = .true.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!


  if (.not.allocated(rc_cub_gpu)) then
    lb1 = lbound(psi_bicub%rc_cub,1)
    ub1 = ubound(psi_bicub%rc_cub,1)
    allocate( rc_cub_gpu(lb1:ub1), stat=ierr )
    call assert(ierr.eq.0,'alloc(rc_cub_gpu)',ierr)

    rc_cub_gpu = 0
  endif

  if (.not.allocated(zc_cub_gpu)) then
    lb1 = lbound(psi_bicub%zc_cub,1)
    ub1 = ubound(psi_bicub%zc_cub,1)
    allocate( zc_cub_gpu(lb1:ub1),stat=ierr )
    call assert(ierr.eq.0,'alloc(zc_cub_gpu)',ierr)

    zc_cub_gpu = 0
  endif

  if (.not.allocated(acoef_cub_gpu)) then
    lb1 = lbound(psi_bicub%acoef_cub,1)
    lb2 = lbound(psi_bicub%acoef_cub,2)
    lb3 = lbound(psi_bicub%acoef_cub,3)
    lb4 = lbound(psi_bicub%acoef_cub,4)

    ub1 = ubound(psi_bicub%acoef_cub,1)
    ub2 = ubound(psi_bicub%acoef_cub,2)
    ub3 = ubound(psi_bicub%acoef_cub,3)
    ub4 = ubound(psi_bicub%acoef_cub,4)

    allocate( acoef_cub_gpu(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),stat=ierr )
    call assert(ierr.eq.0,'alloc(acoef_cub_gpu)',ierr)

    acoef_cub_gpu = 0
  endif


  
  dr_inv_gpu = psi_bicub%dr_inv
  dz_inv_gpu = psi_bicub%dz_inv
  rmin_gpu = psi_bicub%rmin
  zmin_gpu = psi_bicub%zmin 
  nr_gpu=psi_bicub%nr
  nz_gpu=psi_bicub%nz

 if(use_cudaMemcpy) then


     icount = size(rc_cub_gpu)
     ierr = cudaMemcpy( rc_cub_gpu, psi_bicub%rc_cub, &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(rc_cub_gpu)',ierr)

     icount = size(zc_cub_gpu)
     ierr = cudaMemcpy( zc_cub_gpu, psi_bicub%zc_cub, &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(zc_cub_gpu)',ierr)


     icount = size(acoef_cub_gpu)
     ierr = cudaMemcpy( acoef_cub_gpu, psi_bicub%acoef_cub, &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(acoef_cub_gpu)',ierr)

 else

     rc_cub_gpu = psi_bicub%rc_cub

     zc_cub_gpu = psi_bicub%zc_cub

     acoef_cub_gpu = psi_bicub%acoef_cub
     
 endif


    return
  end subroutine update_device_bicub

end module bicub_mod_gpu





module ptl_module_gpu
use dimensions_mod_gpu
use precision_mod_gpu
use ptl_module, only : &
   ptl_nphase_host => ptl_nphase, &
   ptl_nphase2_host => ptl_nphase2, &
   ptl_nconst_host => ptl_nconst, &
   ptl_nsp_max_host => ptl_nsp_max, &
   pir_host => pir, &
   piz_host => piz, &
   pip_host => pip, &
   pirho_host => pirho, & 
   pim_host => pim, &
   piw1_host => piw1, &
   piw2_host => piw2, &
   piw0_host => piw0

use util_mod_gpu, only : get_gpu_streamid
use cudafor
implicit none

! logical, parameter :: use_ptl_block_copy = .false.

  integer, parameter :: ptl_nphase = ptl_nphase_host
  integer, parameter :: ptl_nphase2 = ptl_nphase2_host
  integer, parameter :: ptl_nconst = ptl_nconst_host
  integer, parameter :: ptl_nsp_max = ptl_nsp_max_host

  integer, parameter :: pir = pir_host
  integer, parameter :: piz = piz_host
  integer, parameter :: pip = pip_host
  integer, parameter :: pirho = pirho_host
  integer, parameter :: pim = pim_host
  integer, parameter :: piw1 = piw1_host
  integer, parameter :: piw2 = piw2_host
  integer, parameter :: piw0 = piw0_host

! integer, parameter :: spe_num = spe_num_host
! integer, parameter :: spe_maxnum = spe_maxnum_host
! integer, parameter :: spe_type = spe_type_host

  integer (kind=work_p),allocatable :: sp_ptl_gid(:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_ptl_ph(:,:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_ptl_ct(:,:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_phase0(:,:) !!! new CPU arrays


  logical, parameter :: use_update_phase0 = .true.

  integer :: ptl_isp, ptl_nsp
  integer :: num_gpu,maxnum_gpu,type_gpu
  integer(8),allocatable :: ptl_gid_gpu(:)
  real (kind=work_p),allocatable :: ptl_ph_gpu(:,:)
  real (kind=work_p),allocatable :: ptl_ct_gpu(:,:)
  real (kind=work_p),allocatable :: phase0_gpu(:,:)
  real (kind=work_p),allocatable :: rhoi_gpu(:)

  real(kind=work_p), dimension(0:ptl_nsp_max) :: ptl_mass, ptl_c_m, ptl_c2_2m, ptl_charge
  logical,dimension(0:ptl_nsp_max) :: ptl_deltaf_sp
  integer,parameter :: idebug=0


  attributes(device) :: ptl_gid_gpu, ptl_ph_gpu, ptl_ct_gpu, phase0_gpu, rhoi_gpu
  attributes(device) :: num_gpu,maxnum_gpu,type_gpu
  attributes(constant) :: ptl_mass, ptl_c_m, ptl_c2_2m, ptl_charge, ptl_deltaf_sp
  attributes(device) :: ptl_isp, ptl_nsp

  contains

  attributes(host) &
  subroutine update_device_species_type(sp, gpu_ibegin,gpu_iend)
  use sml_module, only : sml_mype
  use ptl_module, only : &
      species_type 

implicit none

  type(species_type),intent(in) :: sp
  integer, intent(in) :: gpu_ibegin, gpu_iend !!! isp=0 electron, isp=1 ion

! integer :: lb1,ub1,lb2,ub2, llb,uub
  integer :: ierr
  integer :: i, nn, nc, np , iphase
  logical, parameter :: use_cudaMemcpy = .true.
  logical :: isvalid

  integer :: icount,streamid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! convert the ptl_type data structure !!!


  streamid = get_gpu_streamid()


  allocate(sp_ptl_ct(gpu_ibegin:gpu_iend,ptl_nconst), &
           sp_ptl_ph(gpu_ibegin:gpu_iend,ptl_nphase), &
           sp_ptl_gid(gpu_ibegin:gpu_iend), stat=ierr )
 

  do nn = gpu_ibegin, gpu_iend
     sp_ptl_gid(nn)=sp%ptl(nn)%gid
     do np = 1, ptl_nphase
        sp_ptl_ph(nn,np) = sp%ptl(nn)%ph(np)
     enddo
     do nc = 1, ptl_nconst
        sp_ptl_ct(nn,nc) = sp%ptl(nn)%ct(nc)
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!

  if (.not.allocated(ptl_gid_gpu)) then
     allocate( ptl_gid_gpu( gpu_ibegin:gpu_iend), stat=ierr)
     call assert(ierr.eq.0,'alloc(ptl_gid_gpu)',ierr)
  endif

  if (.not.allocated( ptl_ph_gpu)) then
    allocate( ptl_ph_gpu( gpu_ibegin:gpu_iend,ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(ptl_ph_gpu)',ierr)
  endif

  if (.not.allocated(ptl_ct_gpu)) then
    allocate( ptl_ct_gpu( gpu_ibegin:gpu_iend,ptl_nconst),stat=ierr)
    call assert(ierr.eq.0,'alloc(ptl_ct_gpu)',ierr)
  endif

  if (.not.allocated(phase0_gpu)) then
    allocate( phase0_gpu( gpu_ibegin:gpu_iend,ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(phase0_gpu)',ierr)
  endif



!! allocate( ptl_gid_gpu( gpu_ibegin:gpu_iend), &
!! ptl_ph_gpu( gpu_ibegin:gpu_iend,ptl_nphase), &
!! ptl_ct_gpu( gpu_ibegin:gpu_iend,ptl_nconst), &
!! phase0_gpu( gpu_ibegin:gpu_iend,ptl_nphase) )

  num_gpu = sp%num
  maxnum_gpu = sp%maxnum
  type_gpu = sp%type

     isvalid = size(ptl_gid_gpu).eq.size(sp_ptl_gid)
     call assert(isvalid,'invalid size(ptl_gid_gpu)',size(ptl_gid_gpu))

     isvalid = size(ptl_ph_gpu).eq.size(sp_ptl_ph)
     call assert(isvalid,'invalid size(ptl_ph_gpu)',size(ptl_ph_gpu))

     isvalid = size(ptl_ct_gpu).eq.size(sp_ptl_ct)
     call assert( isvalid,'invalid size(ptl_ct_gpu)',size(ptl_ct_gpu))

 if(use_cudaMemcpy) then


     icount = size(ptl_gid_gpu)
     ierr = cudaMemcpy( ptl_gid_gpu, sp_ptl_gid, &
                icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_gid_gpu)',ierr)




     icount = size(ptl_ph_gpu)

     ierr = cudaMemcpy( ptl_ph_gpu, sp_ptl_ph, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_ph_gpu)',ierr)




     icount = size(ptl_ct_gpu)
     ierr = cudaMemcpy( ptl_ct_gpu, sp_ptl_ct, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_ct_gpu)',ierr)

! if(sp%type/=0) then
! allocate(rhoi_gpu(gpu_ibegin:gpu_iend))
! ierr = cudaMemcpyAsync( rhoi_gpu, sp%rhoi, &
! icount_gid, cudaMemcpyHostToDevice,streamid)
! endif
 else
     ptl_gid_gpu(gpu_ibegin:gpu_iend) = sp_ptl_gid(gpu_ibegin:gpu_iend)
     ptl_ph_gpu(gpu_ibegin:gpu_iend,1:ptl_nphase) = sp_ptl_ph(gpu_ibegin:gpu_iend,1:ptl_nphase)
     ptl_ct_gpu(gpu_ibegin:gpu_iend,1:ptl_nconst) = sp_ptl_ct(gpu_ibegin:gpu_iend,1:ptl_nconst)

! if(type_gpu/=0) then
! allocate(rhoi_gpu(gpu_ibegin:gpu_iend))
! rhoi_gpu(gpu_ibegin:gpu_iend) = sp%rhoi(gpu_ibegin:gpu_iend)
! endif

 endif

  
  deallocate(sp_ptl_gid, stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_gid)',ierr)

  deallocate(sp_ptl_ph, stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_ph)',ierr)

  deallocate(sp_ptl_ct,stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_ct)',ierr)



! ---------------------
! update phase0 array?
! ---------------------

  if (use_update_phase0) then
    allocate( sp_phase0(gpu_ibegin:gpu_iend,1:ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(sp_phase0)',ierr)

    do nn=gpu_ibegin,gpu_iend
     sp_phase0(nn,1:ptl_nphase) = sp%phase0(1:ptl_nphase,nn)
    enddo

    if (use_cudaMemcpy) then
       icount = size(phase0_gpu)
       ierr = cudaMemcpy( phase0_gpu, sp_phase0, &
                 icount, cudaMemcpyHostToDevice)
       call assert(ierr.eq.0,'cudaMemcpy(phase0_gpu)',ierr)

    else
      phase0_gpu = sp_phase0
    endif





    deallocate( sp_phase0,stat=ierr)
    call assert(ierr.eq.0,'dealloc(sp_phase0)',ierr)
  endif


  return
  end subroutine update_device_species_type


  attributes(host) &
  subroutine update_host_species_type(sp,gpu_ibegin,gpu_iend)
  use sml_module, only : sml_mype

   use ptl_module, only : species_type, ptl_type

  integer, intent(in) :: gpu_ibegin, gpu_iend
  type(species_type),intent(inout) :: sp
  integer :: ierr, stat
! logical :: isvalid
  integer :: i, nn, np, nc
  integer :: icount_gid,icount_ph,icount_ct
  integer, parameter :: idebug = 0

  logical, parameter :: use_cudaMemcpy = .true.
  integer :: icount, streamid


  streamid = get_gpu_streamid()

  allocate(sp_ptl_ct(gpu_ibegin:gpu_iend,ptl_nconst), &
           sp_ptl_ph(gpu_ibegin:gpu_iend,ptl_nphase), &
           sp_ptl_gid(gpu_ibegin:gpu_iend), stat=ierr )
  call assert(ierr.eq.0,'alloc(sp_ptl_ct)',ierr)

  icount_ph = ptl_nphase*(gpu_iend-gpu_ibegin+1)
  icount_gid = gpu_iend-gpu_ibegin+1
  icount_ct = ptl_nconst*(gpu_iend-gpu_ibegin+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from GPU to CPU !!!
  
  sp%num = num_gpu
  sp%maxnum = maxnum_gpu
  sp%type = type_gpu


! sp%tr_save(gpu_ibegin:gpu_iend) = tr_save_gpu(gpu_ibegin:gpu_iend)
! sp%p_save(1:3,gpu_ibegin:gpu_iend) = p_save_gpu(1:3,gpu_ibegin:gpu_iend)

 if(use_cudaMemcpy) then

     icount = size(sp_ptl_gid)
     ierr = cudaMemcpy( sp_ptl_gid, ptl_gid_gpu, &
               icount, cudaMemcpyDeviceToHost)
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_gid)', ierr )



     icount = size(sp_ptl_ph)
     ierr = cudaMemcpy( sp_ptl_ph, ptl_ph_gpu, &
               icount, cudaMemcpyDeviceToHost)
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_ph)', ierr )




     icount = size(sp_ptl_ct)
     ierr = cudaMemcpy( sp_ptl_ct, ptl_ct_gpu, &
               icount, cudaMemcpyDeviceToHost )
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_ct)', ierr )

! if(sp%type/=0) then
! ierr = cudaMemcpyAsync( sp%rhoi, rhoi_gpu, &
! icount_gid, cudaMemcpyDeviceToHost,streamid)
! endif
 else
     sp_ptl_gid = ptl_gid_gpu
     sp_ptl_ph = ptl_ph_gpu
     sp_ptl_ct = ptl_ct_gpu

! if(sp%type/=0) then
! sp%rhoi(gpu_ibegin:gpu_iend) = rhoi_gpu(gpu_ibegin:gpu_iend)
! endif
 endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! convert the simple array to ptl_type data structure !!!

  do nn = gpu_ibegin, gpu_iend
     sp%ptl(nn)%gid = sp_ptl_gid(nn)
     do np = 1,ptl_nphase
        sp%ptl(nn)%ph(np) = sp_ptl_ph(nn,np)
     enddo
     do nc = 1,ptl_nconst
        sp%ptl(nn)%ct(nc) = sp_ptl_ct(nn,nc)
     enddo
  enddo

  deallocate(ptl_gid_gpu, ptl_ph_gpu, ptl_ct_gpu,stat=ierr) 
  call assert(ierr.eq.0,'dealloc(ptl_gid_gpu)',ierr)

  deallocate(sp_ptl_gid, sp_ptl_ph, sp_ptl_ct,stat=ierr )
  call assert(ierr.eq.0,'dealloc(sp_ptl_gid)',ierr)

  if (use_update_phase0) then

    allocate( sp_phase0(gpu_ibegin:gpu_iend,1:ptl_nphase), stat=ierr)
    call assert(ierr.eq.0,'alloca(sp_phase0)',ierr)


    if (use_cudaMemcpy) then
       icount = size(sp_phase0)
       ierr = cudaMemcpy( sp_phase0, phase0_gpu, &
                 icount, cudaMemcpyDeviceToHost )
       call assert(ierr.eq.0,'cudaMemcpy(sp_phase0)',ierr)
    else 
       sp_phase0 = phase0_gpu
    endif


    do nn=gpu_ibegin,gpu_iend
      sp%phase0(1:ptl_nphase, nn) = sp_phase0(nn,1:ptl_nphase)
    enddo

    deallocate( sp_phase0, stat=ierr)
    call assert(ierr.eq.0,'dealloc(sp_phase0)',ierr)

  endif

  deallocate(phase0_gpu,stat=ierr)
  call assert(ierr.eq.0,'dealloc(phase0_gpu)',ierr)

  return
  end subroutine update_host_species_type


  attributes(host) &
  subroutine update_device_ptl( )
  use ptl_module, only : &
    ptl_mass_host => ptl_mass, &
    ptl_c_m_host => ptl_c_m, &
    ptl_charge_host => ptl_charge, &
    ptl_c2_2m_host => ptl_c2_2m, &
    ptl_deltaf_sp_host => ptl_deltaf_sp, &
    ptl_isp_host => ptl_isp, &
    ptl_nsp_host => ptl_nsp


    integer :: lb,ub
    logical :: isvalid

    lb = lbound(ptl_mass_host,1)
    ub = ubound(ptl_mass_host,1)
    isvalid = (lbound(ptl_mass,1) <= lb) 
    call assert(isvalid, &
       'ptl_module_gpu: invalid lb ', lb)
    isvalid = (ub <= ubound(ptl_mass,1) )
    call assert(isvalid, &
       'ptl_module_gpu: invalid ub ', ub)

    ptl_mass(lb:ub) = ptl_mass_host(lb:ub)
    ptl_c_m(lb:ub) = ptl_c_m_host(lb:ub)
    ptl_c2_2m(lb:ub) = ptl_c2_2m_host(lb:ub)
    ptl_charge(lb:ub) = ptl_charge_host(lb:ub)
    ptl_deltaf_sp(lb:ub) = ptl_deltaf_sp_host(lb:ub)
    ptl_isp = ptl_isp_host
    ptl_nsp = ptl_nsp_host
  return
  end subroutine update_device_ptl

end module ptl_module_gpu
module grid_class_gpu
  use dimensions_mod_gpu
  use precision_mod_gpu
  use grid_class, only :grid_type

  use cudafor
  implicit none

  integer :: grid_nnode !! number of nodes in a grid system
  integer :: grid_ntriangle !! number of trianble in a grid system
  integer :: grid_iphi_offset !! toroidal angle indexing offset
  integer :: grid_nrho


  integer,allocatable :: grid_nd(:,:) !! 3 node numbers that each triangle has
  integer,allocatable :: grid_adj(:,:) !! 3 adjacent triangle
  integer,allocatable :: grid_guess_table(:,:)
  integer,allocatable :: grid_guess_list(:)
  integer,allocatable :: grid_guess_xtable(:,:)
  integer,allocatable :: grid_guess_count(:,:)

  integer, allocatable :: grid_rgn(:) !! region value for each node point
  real (kind=work_p),allocatable :: grid_x(:,:)
 
  real (kind=work_p),allocatable :: grid_mapping(:,:,:) !! shape function coefficient












  integer,dimension(2) :: grid_guess_n
  real (kind=work_p),dimension(2) :: grid_guess_min 
  real (kind=work_p),dimension(2) :: grid_guess_max
  real (kind=work_p),dimension(2) :: grid_guess_d 
  real (kind=work_p),dimension(2) :: grid_inv_guess_d




  real (kind=work_p) :: grid_delta_phi
  real (kind=work_p) :: grid_rhomax
  real (kind=work_p) :: grid_drho
  integer, parameter :: idebug = 0


  attributes(device) :: grid_nd, grid_adj 
  attributes(device) :: grid_guess_table, grid_guess_list
  attributes(device) :: grid_guess_xtable, grid_guess_count
  attributes(device) :: grid_rgn, grid_x, grid_mapping

  attributes(device) :: grid_guess_n
  attributes(device) :: grid_iphi_offset, grid_nrho 
  attributes(device) :: grid_nnode, grid_ntriangle, grid_guess_min, &
                        grid_guess_max, grid_guess_d, grid_inv_guess_d, &
                        grid_delta_phi, grid_rhomax, grid_drho 

  contains

  attributes(host) &
  subroutine update_device_grid_type(grid )
  use assert_mod
  use sml_module, only : sml_mype
  use grid_class, only : grid_type
  type(grid_type), intent(in) :: grid
  logical, parameter :: use_cudaMemcpy = .true.
  integer :: icount_nd, icount_adj, icount_gn, icount_map, icount_gt, icount_gx, icount_gc, &
             icount_list, icount_gmin, icount_gmax, icount_gd, icount_igd, icount_rg, icount_x, ierr
  integer, parameter :: iundefined = huge(0)

  integer :: icount 
  integer :: lb1,ub1, lb2,ub2, lb3,ub3
  logical :: isvalid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy grid data from CPU to GPU !!!

! icount_nd = 3*ntriangle_dim
! icount_adj = 3*ntriangle_dim
! icount_map = 2*3*ntriangle_dim
! icount_gt = guess_table_n1_dim*guess_table_n2_dim
! icount_list = guess_list_size
! icount_gx = guess_table_n1_dim*guess_table_n2_dim 
! icount_gc = guess_table_n1_dim*guess_table_n2_dim
   
!! allocate(grid_nd(3,ntriangle_dim), &
!! grid_adj(3,ntriangle_dim),
!! grid_guess_table(guess_table_n1_dim,guess_table_n2_dim), &
!! grid_guess_list(guess_list_size), &
!! grid_guess_xtable(guess_table_n1_dim,guess_table_n2_dim), &
!! grid_guess_count(guess_table_n1_dim,guess_table_n2_dim), &
!! grid_mapping(2,3,ntriangle_dim), &

  
! =====================
! allocate grid_nd(:,:)
! =====================

  if (.not.allocated(grid_nd)) then
    lb1 = lbound(grid%nd,1)
    lb2 = lbound(grid%nd,2)

    ub1 = ubound(grid%nd,1)
    ub2 = ubound(grid%nd,2)

    allocate( grid_nd( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_ND)',ierr)

  endif

  isvalid = size(grid_nd).eq.size(grid%nd)
  call assert(isvalid,'invalid size(GRID_ND)',size(grid_nd))



! =====================
! allocate grid_adj(:,:)
! =====================

  if (.not.allocated(grid_adj)) then
    lb1 = lbound(grid%adj,1)
    lb2 = lbound(grid%adj,2)

    ub1 = ubound(grid%adj,1)
    ub2 = ubound(grid%adj,2)

    allocate( grid_adj( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_ADJ)',ierr)

  endif

  isvalid = size(grid_adj).eq.size(grid%adj)
  call assert(isvalid,'invalid size(GRID_ADJ)',size(grid_adj))



! =====================
! allocate grid_guess_table(:,:)
! =====================

  if (.not.allocated(grid_guess_table)) then
    lb1 = lbound(grid%guess_table,1)
    lb2 = lbound(grid%guess_table,2)

    ub1 = ubound(grid%guess_table,1)
    ub2 = ubound(grid%guess_table,2)

    allocate( grid_guess_table( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_TABLE)',ierr)

  endif

  isvalid = size(grid_guess_table).eq.size(grid%guess_table)
  call assert(isvalid,'invalid size(GRID_GUESS_TABLE)',size(grid_guess_table))



! =====================
! allocate grid_guess_xtable(:,:)
! =====================

  if (.not.allocated(grid_guess_xtable)) then
    lb1 = lbound(grid%guess_xtable,1)
    lb2 = lbound(grid%guess_xtable,2)

    ub1 = ubound(grid%guess_xtable,1)
    ub2 = ubound(grid%guess_xtable,2)

    allocate( grid_guess_xtable( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_XTABLE)',ierr)

  endif

  isvalid = size(grid_guess_xtable).eq.size(grid%guess_xtable)
  call assert(isvalid,'invalid size(GRID_GUESS_XTABLE)',size(grid_guess_xtable))




! =====================
! allocate grid_guess_count(:,:)
! =====================

  if (.not.allocated(grid_guess_count)) then
    lb1 = lbound(grid%guess_count,1)
    lb2 = lbound(grid%guess_count,2)

    ub1 = ubound(grid%guess_count,1)
    ub2 = ubound(grid%guess_count,2)

    allocate( grid_guess_count( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_COUNT)',ierr)

  endif

  isvalid = size(grid_guess_count).eq.size(grid%guess_count)
  call assert(isvalid,'invalid size(GRID_GUESS_COUNT)',size(grid_guess_count))

! =====================
! allocate grid_guess_list(:)
! =====================

  if (.not.allocated(grid_guess_list)) then
    lb1 = lbound(grid%guess_list,1)

    ub1 = ubound(grid%guess_list,1)

    allocate( grid_guess_list( lb1:ub1 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_LIST)',ierr)

  endif

  isvalid = size(grid_guess_list).eq.size(grid%guess_list)
  call assert(isvalid,'invalid size(GRID_GUESS_LIST)',size(grid_guess_list))


! =====================
! allocate grid_mapping(:,:,:)
! =====================

  if (.not.allocated(grid_mapping)) then
    lb1 = lbound(grid%mapping,1)
    lb2 = lbound(grid%mapping,2)
    lb3 = lbound(grid%mapping,3)

    ub1 = ubound(grid%mapping,1)
    ub2 = ubound(grid%mapping,2)
    ub3 = ubound(grid%mapping,3)

    allocate( grid_mapping( lb1:ub1, lb2:ub2, lb3:ub3 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_MAPPING)',ierr)

  endif

  isvalid = size(grid_mapping).eq.size(grid%mapping)
  call assert(isvalid,'invalid size(GRID_MAPPING)',size(grid_mapping))







! =====================
! allocate grid_rgn(:)
! =====================

  if (.not.allocated(grid_rgn)) then
    lb1 = lbound(grid%rgn,1)

    ub1 = ubound(grid%rgn,1)

    allocate( grid_rgn( lb1:ub1 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_RGN)',ierr)

  endif

  isvalid = size(grid_rgn).eq.size(grid%rgn)
  call assert(isvalid,'invalid size(GRID_RGN)',size(grid_rgn))

! =====================
! allocate grid_x(:)
! =====================

  if (.not.allocated(grid_x)) then
    lb1 = lbound(grid%x,1)
    lb2 = lbound(grid%x,2)

    ub1 = ubound(grid%x,1)
    ub2 = ubound(grid%x,2)

    allocate( grid_x( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_X)',ierr)

  endif

  isvalid = size(grid_x).eq.size(grid%x)
  call assert(isvalid,'invalid size(GRID_X)',size(grid_x))




  grid_nnode = grid%nnode
  grid_ntriangle = grid%ntriangle
  grid_iphi_offset = grid%iphi_offset
  grid_nrho = grid%nrho
  grid_rhomax = grid%rhomax
  grid_drho = grid%drho
  grid_delta_phi = grid%delta_phi


  if(use_cudaMemcpy) then


     icount_nd = size(grid_nd)
     ierr = cudaMemcpy( grid_nd, grid%nd, &
               icount_nd, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_ND)',ierr)


     icount_adj = size(grid_adj)
     ierr = cudaMemcpy( grid_adj, grid%adj, &
               icount_adj, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_ADJ)',ierr)

     icount_map = size(grid_mapping)
     ierr = cudaMemcpy( grid_mapping, grid%mapping, &
              icount_map, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_MAPPING)',ierr)


     icount_rg = size(grid_rgn)
     ierr = cudaMemcpy( grid_rgn, grid%rgn, &
               icount_rg, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_RGN)',ierr)

     icount_x = size(grid_x)
     ierr = cudaMemcpy( grid_x, grid%x, &
               icount_x, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_X)',ierr)


     icount_gt = size(grid_guess_table)
     ierr = cudaMemcpy( grid_guess_table, grid%guess_table, &
               icount_gt, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_TABLE)',ierr)


     icount_gx = size(grid_guess_xtable)
     ierr = cudaMemcpy( grid_guess_xtable, grid%guess_xtable, &
               icount_gx, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_XTABLE)',ierr)


     icount_gc = size(grid_guess_count)
     ierr = cudaMemcpy( grid_guess_count, grid%guess_count, &
               icount_gc, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_COUNT)',ierr)


     icount_list = size(grid_guess_list)
     ierr = cudaMemcpy( grid_guess_list, grid%guess_list, &
               icount_list, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_LIST)',ierr)



     icount = size(grid_guess_n)
     ierr = cudaMemcpy( grid_guess_n, grid%guess_n, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_n)',ierr)

     icount = size(grid_guess_min)
     ierr = cudaMemcpy( grid_guess_min, grid%guess_min, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_min)',ierr)


     icount = size(grid_guess_max)
     ierr = cudaMemcpy( grid_guess_max, grid%guess_max, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_max)',ierr)


     icount = size(grid_guess_d)
     ierr = cudaMemcpy( grid_guess_d, grid%guess_d, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_d)',ierr)


     icount = size(grid_inv_guess_d)
     ierr = cudaMemcpy( grid_inv_guess_d, grid%inv_guess_d, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_inv_guess_d)',ierr)


  else
     grid_nd = grid%nd
     grid_adj = grid%adj

     grid_mapping = grid%mapping
     grid_x = grid%x
     grid_rgn = grid%rgn

     grid_guess_table = grid%guess_table
     grid_guess_xtable = grid%guess_xtable
     grid_guess_count = grid%guess_count
     grid_guess_list = grid%guess_list


     grid_guess_n = grid%guess_n
     grid_guess_min = grid%guess_min
     grid_guess_max = grid%guess_max
     grid_guess_d = grid%guess_d
     grid_inv_guess_d = grid%inv_guess_d


  endif


  return
  end subroutine update_device_grid_type


  attributes(host) &
  subroutine update_host_grid_type()

  integer :: ierr

  if (allocated(grid_nd)) then
    deallocate(grid_nd,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_ND)',ierr)
  endif


  if (allocated(grid_adj)) then
    deallocate(grid_adj,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_ADJ)',ierr)
  endif



  if (allocated(grid_guess_table)) then
    deallocate(grid_guess_table,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_TABLE)',ierr)
  endif


  if (allocated(grid_guess_xtable)) then
    deallocate(grid_guess_xtable,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_XTABLE)',ierr)
  endif


  if (allocated(grid_guess_count)) then
    deallocate(grid_guess_count,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_COUNT)',ierr)
  endif



  if (allocated(grid_guess_list)) then
    deallocate(grid_guess_list,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_LIST)',ierr)
  endif

  if (allocated(grid_mapping)) then
    deallocate(grid_mapping,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_MAPPING)',ierr)
  endif


  if (allocated(grid_rgn)) then
    deallocate(grid_rgn,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_RGN)',ierr)
  endif


  if (allocated(grid_x)) then
    deallocate(grid_x,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_X)',ierr)
  endif


  return
  end subroutine update_host_grid_type


end module grid_class_gpu

















module psn_class_gpu
use dimensions_mod_gpu
use precision_mod_gpu

use cudafor
!use boundary_class_gpu, only : boundary2_type

  real (kind=work_p),allocatable :: psn_pot_rho_ff(:,:,:) 
  real (kind=work_p),allocatable :: psn_E_rho_ff(:,:,:,:)
  real (kind=work_p),allocatable :: psn_ddpotdt_phi(:,:,:)
  real (kind=work_p),allocatable :: psn_pot0(:)
  real (kind=work_p),allocatable :: psn_sheath_lost(:)
  real (kind=work_p),allocatable :: psn_sheath_pot(:)
  integer, allocatable :: psn_node_to_wall(:)
  integer, allocatable :: psn_wall_nodes(:)


  real (kind=work_p),allocatable :: psn_E_phi_ff(:,:,:,:)

  attributes(device) :: psn_E_phi_ff
 
  attributes(device) :: psn_pot_rho_ff, psn_E_rho_ff, &
                        psn_ddpotdt_phi, psn_pot0, psn_sheath_lost, psn_sheath_pot, &
                        psn_node_to_wall, psn_wall_nodes, psn_nwall !, psn_E00_ff, psn_pot_phi_ff,psn_ddpotdt

  contains

  attributes(host) &
  subroutine update_device_psn_type(psn)
  use sml_module, only : sml_mype
  use psn_class, only : psn_type
  use cudafor
  implicit none

  type(psn_type) :: psn
! real (kind=work_p), parameter :: undefined = huge(0.0_work_p)
! integer :: i3,i4
! integer :: lb1,lb2,lb3,lb4, ub1,ub2,ub3,ub4
  integer :: icount_pr,icount_er, icount_e0, icount_pp, icount_ep, icount_dd
  integer :: ierr
  logical :: use_cudaMemcpy = .true.
  integer, parameter :: idebug = 0 

  integer :: lb1,ub1, lb2,ub2, lb3,ub3, lb4,ub4
  logical :: isvalid
  integer :: icount 


  psn_nwall = psn%nwall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy psn data from CPU to GPU !!!
  
!! allocate(psn_pot_rho_ff(0:1,0:nrho_dim,nnode_dim), 
!! psn_E_rho_ff(3,0:1,0:nrho_dim,nnode_dim), 
!! psn_ddpotdt(nnode_dim,0:1), &
!! psn_E00_ff(2,0:1,nnode_dim), 
!! psn_pot_phi_ff(0:1,nnode_dim,0:(nphi_dim-1)), &
!! psn_E_phi_ff(3,0:1,nnode_dim,0:(nphi_dim-1)),
!! psn_ddpotdt_phi(nnode_dim,0:1,0:(nphi_dim-1)))


! ========================
! allocate psn_pot_rho_ff(:,:,:)
! ========================
  if (.not.allocated(psn_pot_rho_ff)) then

    lb1 = lbound(psn%pot_rho_ff,1)
    lb2 = lbound(psn%pot_rho_ff,2)
    lb3 = lbound(psn%pot_rho_ff,3)

    ub1 = ubound(psn%pot_rho_ff,1)
    ub2 = ubound(psn%pot_rho_ff,2)
    ub3 = ubound(psn%pot_rho_ff,3)

    allocate( psn_pot_rho_ff(lb1:ub1,lb2:ub2,lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_pot_rho_ff)',ierr)
    psn_pot_rho_ff = 0
  endif

  isvalid = size(psn%pot_rho_ff).eq.size(psn_pot_rho_ff)
  call assert(isvalid,'invalid size(psn_pot_rho_ff)',size(psn_pot_rho_ff))

! ======================
! allocate psn_E_rho_ff(:,:,:,:)
! ======================
  if (.not.allocated(psn_E_rho_ff)) then

    lb1 = lbound(psn%E_rho_ff,1)
    lb2 = lbound(psn%E_rho_ff,2)
    lb3 = lbound(psn%E_rho_ff,3)
    lb4 = lbound(psn%E_rho_ff,4)

    ub1 = ubound(psn%E_rho_ff,1)
    ub2 = ubound(psn%E_rho_ff,2)
    ub3 = ubound(psn%E_rho_ff,3)
    ub4 = ubound(psn%E_rho_ff,4)

    allocate( psn_E_rho_ff(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_E_rho_ff)',ierr)

    psn_E_rho_ff = 0
  endif

  isvalid = size(psn_E_rho_ff).eq.size(psn%E_rho_ff)
  call assert(isvalid,'invalid size(psn_E_rho_ff)',size(psn_E_rho_ff))

! ======================
! allocate psn_pot0(:)
! ======================
  if (.not.allocated(psn_pot0)) then

    lb1 = lbound(psn%pot0,1)

    ub1 = ubound(psn%pot0,1)

    allocate( psn_pot0(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_pot0)',ierr)

    psn_pot0 = 0
  endif

  isvalid = size(psn_pot0).eq.size(psn%pot0)
  call assert(isvalid,'invalid size(psn_pot0)',size(psn_pot0))

! ======================
! allocate psn_sheath_lost(:)
! ======================
  if (.not.allocated(psn_sheath_lost)) then

    lb1 = lbound(psn%sheath_lost(:,1),1)

    ub1 = ubound(psn%sheath_lost(:,1),1)
!print *, 'sheath_lost', lb1, ub1, psn%nwall, size(psn%sheath_lost)
    if(ub1>=lb1) allocate( psn_sheath_lost(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_sheath_lost)',ierr)

    psn_sheath_lost = 0
  endif

  isvalid = size(psn_sheath_lost).eq.size(psn%sheath_lost(:,1))
  call assert(isvalid,'invalid size(psn_sheath_lost)',size(psn_sheath_lost))

! ======================
! allocate psn_sheath_pot(:)
! ======================
  if (.not.allocated(psn_sheath_pot)) then

    lb1 = lbound(psn%sheath_pot,1)

    ub1 = ubound(psn%sheath_pot,1)
!print *, 'sheath_lost', lb1, ub1, psn%nwall, size(psn%sheath_lost)
    if(ub1>=lb1) allocate( psn_sheath_pot(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_sheath_pot)',ierr)

    psn_sheath_pot = 0
  endif

  isvalid = size(psn_sheath_pot).eq.size(psn%sheath_pot)
  call assert(isvalid,'invalid size(psn_sheath_pot)',size(psn_sheath_pot))

! ======================
! allocate psn_wall_nodes(:)
! ======================
  if (.not.allocated(psn_wall_nodes)) then

    lb1 = lbound(psn%wall_nodes,1)

    ub1 = ubound(psn%wall_nodes,1)

    if(ub1>=lb1) allocate( psn_wall_nodes(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_wall_nodes)',ierr)

    psn_wall_nodes = 0
  endif

  isvalid = size(psn_wall_nodes).eq.size(psn%wall_nodes)
  call assert(isvalid,'invalid size(psn_wall_nodes)',size(psn_wall_nodes))

! ======================
! allocate psn_node_to_wall(:)
! ======================
  if (.not.allocated(psn_node_to_wall)) then

    lb1 = lbound(psn%node_to_wall,1)

    ub1 = ubound(psn%node_to_wall,1)

    if(ub1>=lb1) allocate( psn_node_to_wall(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_node_to_wall)',ierr)

    psn_node_to_wall = 0
  endif

  isvalid = size(psn_node_to_wall).eq.size(psn%node_to_wall)
  call assert(isvalid,'invalid size(psn_node_to_wall)',size(psn_node_to_wall))


! =====================
! allocate psn_ddpotdt(:,:)
! =====================

! if (.not.allocated(psn_ddpotdt)) then
  
! lb1 = lbound(psn%ddpotdt,1)
! lb2 = lbound(psn%ddpotdt,2)

! ub1 = ubound(psn%ddpotdt,1)
! ub2 = ubound(psn%ddpotdt,2)

! allocate( psn_ddpotdt(lb1:ub1, lb2:ub2),stat=ierr)
! call assert(ierr.eq.0,'alloc(psn_ddpotdt)',ierr)

! psn_ddpotdt = 0
! endif

! isvalid = size(psn%ddpotdt).eq.size(psn_ddpotdt)
! call assert(isvalid,'invalid size(psn_ddpotdt)',size(psn_ddpotdt))


! =====================
! allocate psn_E00_ff(:,:,:)
! =====================

! if (.not.allocated(psn_E00_ff)) then

! lb1 = lbound(psn%E00_ff,1)
! lb2 = lbound(psn%E00_ff,2)
! lb3 = lbound(psn%E00_ff,3)

! ub1 = ubound(psn%E00_ff,1)
! ub2 = ubound(psn%E00_ff,2)
! ub3 = ubound(psn%E00_ff,3)

! allocate( psn_E00_ff(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
! call assert(ierr.eq.0,'alloc(psn_E00_ff)',ierr)

! psn_E00_ff = 0
! endif

! isvalid = size(psn%E00_ff).eq.size(psn_E00_ff)
! call assert(isvalid,'invalid size(psn_E00_ff)',size(psn_E00_ff))


! =====================
! allocate psn_pot_phi_ff(:,:,:)
! =====================

! if (.not.allocated(psn_pot_phi_ff)) then

! lb1 = lbound(psn%pot_phi_ff,1)
! lb2 = lbound(psn%pot_phi_ff,2)
! lb3 = lbound(psn%pot_phi_ff,3)

! ub1 = ubound(psn%pot_phi_ff,1)
! ub2 = ubound(psn%pot_phi_ff,2)
! ub3 = ubound(psn%pot_phi_ff,3)

! allocate( psn_pot_phi_ff(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
! call assert(ierr.eq.0,'alloc(psn_pot_phi_ff)',ierr)

! psn_pot_phi_ff = 0
! endif

! isvalid = size(psn%pot_phi_ff).eq.size(psn_pot_phi_ff)
! call assert(isvalid,'invalid size(psn_pot_phi_ff)',size(psn_pot_phi_ff))


! =====================
! allocate psn_E_phi_ff(:,:,:,:)
! =====================

  if (.not.allocated(psn_E_phi_ff)) then

    lb1 = lbound(psn%E_phi_ff,1)
    lb2 = lbound(psn%E_phi_ff,2)
    lb3 = lbound(psn%E_phi_ff,3)
    lb4 = lbound(psn%E_phi_ff,4)

    ub1 = ubound(psn%E_phi_ff,1)
    ub2 = ubound(psn%E_phi_ff,2)
    ub3 = ubound(psn%E_phi_ff,3)
    ub4 = ubound(psn%E_phi_ff,4)

    allocate( psn_E_phi_ff(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_E_phi_ff)',ierr)

    psn_E_phi_ff = 0
  endif

  isvalid = size(psn%E_phi_ff).eq.size(psn_E_phi_ff)
  call assert(isvalid,'invalid size(psn_E_phi_ff)',size(psn_E_phi_ff))




    
! =====================
! allocate psn_ddpotdt_phi(:,:,:)
! =====================

  if (.not.allocated(psn_ddpotdt_phi)) then

    lb1 = lbound(psn%ddpotdt_phi,1)
    lb2 = lbound(psn%ddpotdt_phi,2)
    lb3 = lbound(psn%ddpotdt_phi,3)

    ub1 = ubound(psn%ddpotdt_phi,1)
    ub2 = ubound(psn%ddpotdt_phi,2)
    ub3 = ubound(psn%ddpotdt_phi,3)

    allocate( psn_ddpotdt_phi(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ddpotdt_phi)',ierr)

    psn_ddpotdt_phi = 0
  endif

  isvalid = size(psn%ddpotdt_phi).eq.size(psn_ddpotdt_phi)
  call assert(isvalid,'invalid size(psn_ddpotdt_phi)',size(psn_ddpotdt_phi))







!! icount_pr = 2*(nrho_dim+1)*nnode_dim
!! icount_er = 3*2*(nrho_dim+1)*nnode_dim
!! icount_e0 = nnode_dim*2
!! icount_pp = 2*nnode_dim*nphi_dim
!! icount_ep = 3*2*nnode_dim*nphi_dim
!! icount_dd = nnode_dim*2*nphi_dim
  
  if(use_cudaMemcpy) then
     icount = size(psn_pot_rho_ff)
     ierr = cudaMemcpy( psn_pot_rho_ff, psn%pot_rho_ff, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_pot_rho_ff)',ierr)

     icount = size(psn_E_rho_ff)
     ierr = cudaMemcpy( psn_E_rho_ff, psn%E_rho_ff, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_E_rho_ff)',ierr)

     icount = size(psn_pot0)
     ierr = cudaMemcpy( psn_pot0, psn%pot0, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_pot0)',ierr)

     icount = size(psn_wall_nodes)
     ierr = cudaMemcpy( psn_wall_nodes, psn%wall_nodes, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_wall_nodes)',ierr)

! icount = size(psn_sheath_lost)
! ierr = cudaMemcpy( psn_sheath_lost, psn%sheath_lost, &
! icount, cudaMemcpyHostToDevice)
! call assert(ierr.eq.0, 'cudaMemcpy(psn_sheath_lost)',ierr)

     icount = size(psn_sheath_pot)
     ierr = cudaMemcpy( psn_sheath_pot, psn%sheath_pot, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_sheath_pot)',ierr)

     icount = size(psn_node_to_wall)
     ierr = cudaMemcpy( psn_node_to_wall, psn%node_to_wall, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_node_wall)',ierr)


! icount = size(psn_ddpotdt)
! ierr = cudaMemcpy( psn_ddpotdt, psn%ddpotdt, &
! icount, cudaMemcpyHostToDevice)
! call assert(ierr.eq.0, 'cudaMemcpy(psn_ddpotdt)',ierr)

! icount = size(psn_E00_ff)
! ierr = cudaMemcpy( psn_E00_ff, psn%E00_ff, &
! icount, cudaMemcpyHostToDevice)
! call assert(ierr.eq.0, 'cudaMemcpy(psn_E00_ff)',ierr)



! icount = size(psn_pot_phi_ff)
! ierr = cudaMemcpy( psn_pot_phi_ff, psn%pot_phi_ff, &
! icount, cudaMemcpyHostToDevice)
! call assert(ierr.eq.0,'cudaMemcpy(psn_pot_phi_ff)',ierr)


     icount = size(psn_E_phi_ff)
     ierr = cudaMemcpy( psn_E_phi_ff, psn%E_phi_ff, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_E_phi_ff)',ierr)




     icount = size(psn_ddpotdt_phi)
     ierr = cudaMemcpy( psn_ddpotdt_phi, psn%ddpotdt_phi, &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ddpotdt_phi)',ierr)
  else
     psn_pot_rho_ff = psn%pot_rho_ff
     psn_E_rho_ff = psn%E_rho_ff
! psn_ddpotdt = psn%ddpotdt
! psn_E00_ff = psn%E00_ff

! psn_pot_phi_ff = psn%pot_phi_ff

     psn_E_phi_ff = psn%E_phi_ff
     psn_ddpotdt_phi = psn%ddpotdt_phi
  endif
 


  return
  end subroutine update_device_psn_type



  attributes(host) &
  subroutine update_host_psn_type(psn)
  use psn_class, only : psn_type
  use cudafor
  implicit none
  integer, parameter :: idebug = 0
  integer :: ierr, icount 
  type(psn_type):: psn
  real (8) :: tmp(psn%nwall)
  logical, parameter :: use_cudaMemcpy = .true.

  if(use_cudaMemcpy) then

     !icount = size(psn%sheath_lost)
     icount = size(psn_sheath_lost)
     ierr = cudaMemcpy( tmp, psn_sheath_lost, &
               icount, cudaMemcpyDeviceToHost)
     call assert(ierr.eq.0,'cudaMemcpy(psn%sheath_lost)',ierr)

  else

     tmp = psn_sheath_lost

  endif
  psn%sheath_lost(:,1) = psn%sheath_lost(:,1) + tmp

  if(allocated(psn_pot_rho_ff)) then
    deallocate( psn_pot_rho_ff,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_pot_rho_ff)',ierr)
  endif

  if(allocated(psn_E_rho_ff)) then
    deallocate( psn_E_rho_ff,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_E_rho_ff)',ierr)
  endif

  if(allocated(psn_ddpotdt_phi)) then
    deallocate( psn_ddpotdt_phi,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ddpotdt_phi)',ierr)

  endif

  if(allocated(psn_node_to_wall)) then
    deallocate( psn_node_to_wall,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_node_to_wall)',ierr)
  endif

  if(allocated(psn_wall_nodes)) then
    deallocate( psn_wall_nodes,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_wall_nodes)',ierr)
  endif

  if(allocated(psn_sheath_lost)) then
    deallocate( psn_sheath_lost,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_sheath_lost)',ierr)
  endif

  if(allocated(psn_sheath_pot)) then
    deallocate( psn_sheath_pot,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_sheath_pot)',ierr)
  endif

  if(allocated(psn_pot0)) then
    deallocate( psn_pot0,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_pot0)',ierr)
  endif


  if(allocated(psn_E_phi_ff)) then
    deallocate( psn_E_phi_ff,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_E_phi_ff)',ierr)
  endif


  return
  end subroutine update_host_psn_type


end module psn_class_gpu
module boundary_class_gpu
use dimensions_mod_gpu
implicit none


  type range_type
     integer :: start,end 
  end type range_type

  type boundary2_type
     integer :: nseg
     type(range_type) :: iseg(iseg_dim)
  end type boundary2_type


end module boundary_class_gpu
!! Diagnosis module
module diag_module_gpu
   use dimensions_mod_gpu, only : nthreads_dim
   use ptl_module_gpu,only : ptl_nsp_max

   use util_mod_gpu 
   use precision_mod_gpu
   use diag_module, only : &
      diag_max_sp_num_host => diag_max_sp_num, &
      diag_1d_npv1_host => diag_1d_npv1, &
      diag_heat_nvar_host => diag_heat_nvar
   implicit none

   integer, parameter :: diag_max_sp_num=diag_max_sp_num_host ! two species

  ! tracer 
! integer :: diag_tracer_period !! Period for tracer
! integer :: diag_tracer_n !! Particle index for tracer routine
! integer :: diag_tracer_sp ! Species index for tracer 

  ! for 1D
   integer, parameter :: diag_1d_npv1=diag_1d_npv1_host

! Integer parameters for heat load diagnosis
   integer, parameter :: diag_heat_nvar=diag_heat_nvar_host

   integer :: diag_heat_nr, diag_heat_nz, diag_heat_npsi, diag_heat_nsection 

! logical :: diag_1d_on
   logical :: diag_eflux_on
! integer :: diag_1d_period
   integer :: diag_1d_npsi !
   integer :: diag_1d_ne
   real (kind=work_p) :: diag_1d_pin 
! real (8) :: diag_1d_pout
! real (8) :: diag_1d_dp 
   real (kind=work_p) :: diag_1d_dp_inv
   real (kind=work_p) :: diag_1d_emin
   real (kind=work_p) :: diag_1d_emax
! real (8),allocatable :: diag_1d_vol(:)
! integer :: diag_1d_isp,diag_1d_nsp ! same as ptl_isp/nsp
! logical :: diag_tavg_on
! 
! integer, parameter :: ptl_isp = 0
   integer, parameter :: diag_1d_isp = 0

   ! ------------------------------------------
   ! value for diag_1d_npsi is set in setup.F90
   ! ------------------------------------------
   integer, parameter :: diag_1d_npsi_dim = 80 
   integer, parameter :: diag_1d_nsp = ptl_nsp_max

   integer, parameter :: np = diag_1d_npsi_dim
   integer, parameter :: nsp = diag_1d_nsp
   integer, parameter :: isp = diag_1d_isp
   logical :: diag_heat_on 

   private :: diag_1d_isp, diag_1d_npsi_dim, diag_1d_nsp
   private :: np,nsp,isp

   integer, parameter :: diag_1d_dim4 = 16
   real(kind=work_p), allocatable :: diag_1d_f_pv1(:,:,:,:)
   real(kind=work_p), allocatable :: diag_1d_df_pv1(:,:,:,:)
   real(kind=work_p), allocatable :: diag_1d_eflux_pv(:,:,:,:,:)
   real(kind=work_p), allocatable :: diag_2d_dflux_pv(:,:,:,:,:)

! for heat load diagnosis
   real(kind=work_p), allocatable :: diag_heat_pv(:,:,:,:,:,:), diag_heat_pv_psi(:,:,:,:,:)
   real(kind=work_p), allocatable :: diag_heat_rmax(:), diag_heat_rmin(:), diag_heat_zmax(:), diag_heat_zmin(:), &
                                     diag_heat_dr(:), diag_heat_dz(:), diag_heat_pmax(:), diag_heat_pmin(:), &
                                     diag_heat_dp(:)
! real (8), allocatable :: diag_1d_tavg_f_pv1(:,:,:,:)
! real (8), allocatable :: diag_1d_tavg_df_pv1(:,:,:,:)
!
! logical :: diag_3d_on
! integer :: diag_3d_period 

   attributes(device) :: diag_1d_npsi
   attributes(device) :: diag_1d_pin, diag_1d_dp_inv
   attributes(device) :: diag_1d_f_pv1, diag_1d_df_pv1
   attributes(device) :: diag_heat_nr, diag_heat_nz, diag_heat_npsi, diag_heat_nsection , diag_heat_on
! attributes(device) :: diag_heat_rmax1, diag_heat_rmax2, diag_heat_rmax3, diag_heat_rmin1, diag_heat_rmin2, &
! diag_heat_rmin2, diag_heat_rmin3, diag_heat_zmax1, diag_heat_zmax2, diag_heat_zmax3, &
! diag_heat_zmin1, diag_heat_zmin2, diag_heat_zmin3
   attributes(device) :: diag_heat_rmax, diag_heat_rmin, diag_heat_zmax, diag_heat_zmin, diag_heat_dr, &
                         diag_heat_dz, diag_heat_pmax, diag_heat_pmin, diag_heat_dp
   attributes(device) :: diag_heat_pv, diag_heat_pv_psi
   attributes(device) :: diag_1d_eflux_pv, diag_2d_dflux_pv
   attributes(device) :: diag_1d_ne, diag_eflux_on, diag_1d_emin, diag_1d_emax
   contains



   attributes(host) &
   subroutine update_device_diag()
   use assert_mod
   use sml_module, only : sml_mype

   use diag_module, only : &
      diag_1d_npsi_host => diag_1d_npsi, &
      diag_1d_dp_inv_host => diag_1d_dp_inv, &
      diag_1d_pin_host => diag_1d_pin, &
      diag_1d_df_pv1_host => diag_1d_df_pv1, &
      diag_1d_f_pv1_host => diag_1d_f_pv1, &
      diag_heat_nr_host => diag_heat_nr, &
      diag_heat_nz_host => diag_heat_nz, &
      diag_heat_npsi_host => diag_heat_npsi, &
      diag_heat_nsection_host => diag_heat_nsection, &
      diag_heat_on_host => diag_heat_on, &
! diag_heat_rmax1_host => diag_heat_rmax1, &
! diag_heat_rmax2_host => diag_heat_rmax2, &
! diag_heat_rmax3_host => diag_heat_rmax3, &
! diag_heat_rmin1_host => diag_heat_rmin1, &
! diag_heat_rmin2_host => diag_heat_rmin2, &
! diag_heat_rmin3_host => diag_heat_rmin3, &
! diag_heat_zmax1_host => diag_heat_zmax1, &
! diag_heat_zmax2_host => diag_heat_zmax2, &
! diag_heat_zmax3_host => diag_heat_zmax3, &
! diag_heat_zmin1_host => diag_heat_zmin1, &
! diag_heat_zmin2_host => diag_heat_zmin2, &
! diag_heat_zmin3_host => diag_heat_zmin3, &
      diag_heat_rmax_host => diag_heat_rmax, &
      diag_heat_rmin_host => diag_heat_rmin, &
      diag_heat_zmax_host => diag_heat_zmax, &
      diag_heat_zmin_host => diag_heat_zmin, &
      diag_heat_dr_host => diag_heat_dr, &
      diag_heat_dz_host => diag_heat_dz, &
      diag_heat_pmax_host => diag_heat_pmax, &
      diag_heat_pmin_host => diag_heat_pmin, &
      diag_heat_dp_host => diag_heat_dp, &
      diag_heat_pv_host => diag_heat_pv, &
      diag_heat_pv_psi_host => diag_heat_pv_psi, &
      diag_eflux_on_host => diag_eflux_on, &
      diag_1d_ne_host => diag_1d_ne, &
      diag_1d_emin_host => diag_1d_emin, &
      diag_1d_emax_host => diag_1d_emax, &
      diag_1d_eflux_pv_host => diag_1d_eflux_pv, &
      diag_2d_dflux_pv_host => diag_2d_dflux_pv

   integer, parameter :: idebug = 2

   integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb,ub
   integer :: ierr
   logical :: need_diag_1d_f_pv1
   logical :: need_diag_1d_df_pv1
   logical :: need_diag_heat_pv
   logical :: need_diag_heat_pv_psi
   logical :: need_diag_1d_eflux_pv
   logical :: need_diag_2d_dflux_pv

   diag_1d_npsi = diag_1d_npsi_host
   diag_1d_dp_inv = diag_1d_dp_inv_host
   diag_1d_pin = diag_1d_pin_host

   diag_heat_npsi = diag_heat_npsi_host
   diag_heat_nr = diag_heat_nr_host
   diag_heat_nz = diag_heat_nz_host
   diag_heat_nsection = diag_heat_nsection_host
   diag_heat_on = diag_heat_on_host

! diag_heat_rmax1 = diag_heat_rmax1_host
! diag_heat_rmax2 = diag_heat_rmax2_host
! diag_heat_rmax3 = diag_heat_rmax3_host
! diag_heat_rmin1 = diag_heat_rmin1_host
! diag_heat_rmin2 = diag_heat_rmin2_host
! diag_heat_rmin3 = diag_heat_rmin3_host
! diag_heat_zmax1 = diag_heat_zmax1_host
! diag_heat_zmax2 = diag_heat_zmax2_host
! diag_heat_zmax3 = diag_heat_zmax3_host
! diag_heat_zmin1 = diag_heat_zmin1_host
! diag_heat_zmin2 = diag_heat_zmin2_host
! diag_heat_zmin3 = diag_heat_zmin3_host
   diag_eflux_on =diag_eflux_on_host
   diag_1d_ne =diag_1d_ne_host
   diag_1d_emin =diag_1d_emin_host
   diag_1d_emax =diag_1d_emax_host

   need_diag_1d_f_pv1 = allocated(diag_1d_f_pv1_host)
   if (need_diag_1d_f_pv1) then

! ------------------
! clear previous storage
! ------------------
     if (allocated(diag_1d_f_pv1)) then
      deallocate( diag_1d_f_pv1, stat=ierr)
      call assert(ierr.eq.0,'dealloc(diag_1d_f_pv1),ierr',ierr)
     endif

   lb1 = lbound(diag_1d_f_pv1_host,1)
   lb2 = lbound(diag_1d_f_pv1_host,2)
   lb3 = lbound(diag_1d_f_pv1_host,3)

   ub1 = ubound(diag_1d_f_pv1_host,1)
   ub2 = ubound(diag_1d_f_pv1_host,2)
   ub3 = ubound(diag_1d_f_pv1_host,3)

    
     allocate( diag_1d_f_pv1(lb1:ub1,lb2:ub2,lb3:ub3,diag_1d_dim4), &
             stat=ierr)
     call assert(ierr.eq.0,'alloc(diag_1d_f_pv1),diag_1d_dim4', &
             diag_1d_dim4)

   endif

   need_diag_1d_df_pv1 = allocated(diag_1d_df_pv1_host)
   if (need_diag_1d_df_pv1) then
! ------------------
! clear previous storage
! ------------------
     if (allocated(diag_1d_df_pv1)) then
       deallocate( diag_1d_df_pv1,stat=ierr)
       call assert(ierr.eq.0,'dealloc(diag_1d_df_pv1,ierr=',ierr)
     endif

   lb1 = lbound(diag_1d_df_pv1_host,1)
   lb2 = lbound(diag_1d_df_pv1_host,2)
   lb3 = lbound(diag_1d_df_pv1_host,3)

   ub1 = ubound(diag_1d_df_pv1_host,1)
   ub2 = ubound(diag_1d_df_pv1_host,2)
   ub3 = ubound(diag_1d_df_pv1_host,3)


     allocate( diag_1d_df_pv1(lb1:ub1,lb2:ub2,lb3:ub3,diag_1d_dim4), &
             stat=ierr)
     call assert(ierr.eq.0,'alloc(diag_1d_df_pv1),diag_1d_dim4', &
             diag_1d_dim4)
   endif


   if (allocated(diag_1d_f_pv1)) then
      diag_1d_f_pv1 = 0
   endif

   if (allocated(diag_1d_df_pv1)) then
      diag_1d_df_pv1 = 0
   endif
   
   if (diag_heat_on_host) then
     if (allocated(diag_heat_rmax)) then
        deallocate(diag_heat_rmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_rmax),ierr=',ierr) 
     endif
     lb = lbound(diag_heat_rmax_host,1)
     ub = ubound(diag_heat_rmax_host,1)
     allocate( diag_heat_rmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_rmax),ierr=',ierr) 
     diag_heat_rmax=diag_heat_rmax_host


     if (allocated(diag_heat_rmin)) then 
        deallocate(diag_heat_rmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_rmin),ierr=',ierr)
     endif 
     lb = lbound(diag_heat_rmin_host,1)
     ub = ubound(diag_heat_rmin_host,1)
     allocate( diag_heat_rmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_rmin),ierr=',ierr)
     diag_heat_rmin=diag_heat_rmin_host


     if (allocated(diag_heat_zmax)) then
        deallocate(diag_heat_zmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_zmax),ierr=',ierr)
     endif
     lb = lbound(diag_heat_zmax_host,1)
     ub = ubound(diag_heat_zmax_host,1)
     allocate( diag_heat_zmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_zmax),ierr=',ierr)
     diag_heat_zmax=diag_heat_zmax_host


     if (allocated(diag_heat_zmin)) then
        deallocate(diag_heat_zmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_zmin),ierr=',ierr)
     endif
     lb = lbound(diag_heat_zmin_host,1)
     ub = ubound(diag_heat_zmin_host,1)
     allocate( diag_heat_zmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_zmin),ierr=',ierr)
     diag_heat_zmin=diag_heat_zmin_host


     if (allocated(diag_heat_dr)) then
        deallocate(diag_heat_dr, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dr),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dr_host,1)
     ub = ubound(diag_heat_dr_host,1)
     allocate( diag_heat_dr(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dr),ierr=',ierr)
     diag_heat_dr=diag_heat_dr_host


     if (allocated(diag_heat_dz)) then
        deallocate(diag_heat_dz, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dz),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dz_host,1)
     ub = ubound(diag_heat_dz_host,1)
     allocate( diag_heat_dz(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dz),ierr=',ierr)
    diag_heat_dz=diag_heat_dz_host


     if (allocated(diag_heat_pmax)) then
        deallocate(diag_heat_pmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_pmax),ierr=',ierr)
     endif
     lb = lbound(diag_heat_pmax_host,1)
     ub = ubound(diag_heat_pmax_host,1)
     allocate( diag_heat_pmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_pmax),ierr=',ierr)
     diag_heat_pmax=diag_heat_pmax_host


     if (allocated(diag_heat_pmin)) then
        deallocate(diag_heat_pmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_pmin),ierr=',ierr)
     endif
     lb = lbound(diag_heat_pmin_host,1)
     ub = ubound(diag_heat_pmin_host,1)
     allocate( diag_heat_pmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_pmin),ierr=',ierr)
     diag_heat_pmin=diag_heat_pmin_host


     if (allocated(diag_heat_dp)) then
        deallocate(diag_heat_dp, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dp),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dp_host,1)
     ub = ubound(diag_heat_dp_host,1)
     allocate( diag_heat_dp(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dp),ierr=',ierr)
     diag_heat_dp=diag_heat_dp_host


     need_diag_heat_pv = allocated(diag_heat_pv_host)
     if (need_diag_heat_pv) then
        if (allocated(diag_heat_pv)) then
           deallocate( diag_heat_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_heat_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_heat_pv_host,1)
        lb2 = lbound(diag_heat_pv_host,2)
        lb3 = lbound(diag_heat_pv_host,3)
        lb4 = lbound(diag_heat_pv_host,4)
        lb = lbound(diag_heat_pv_host,5)

        ub1 = ubound(diag_heat_pv_host,1)
        ub2 = ubound(diag_heat_pv_host,2)
        ub3 = ubound(diag_heat_pv_host,3)
        ub4 = ubound(diag_heat_pv_host,4)
        ub = ubound(diag_heat_pv_host,5)
        allocate( diag_heat_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb:ub,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_heat_pv),diag_1d_dim4', &
             diag_1d_dim4,size(diag_heat_pv_host))
     endif


     need_diag_heat_pv_psi = allocated(diag_heat_pv_psi_host)
     if (need_diag_heat_pv_psi) then
        if (allocated(diag_heat_pv_psi)) then
           deallocate( diag_heat_pv_psi,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_heat_pv_psi,ierr=',ierr)
        endif
        lb1 = lbound(diag_heat_pv_psi_host,1)
        lb2 = lbound(diag_heat_pv_psi_host,2)
        lb3 = lbound(diag_heat_pv_psi_host,3)
        lb4 = lbound(diag_heat_pv_psi_host,4)

        ub1 = ubound(diag_heat_pv_psi_host,1)
        ub2 = ubound(diag_heat_pv_psi_host,2)
        ub3 = ubound(diag_heat_pv_psi_host,3)
        ub4 = ubound(diag_heat_pv_psi_host,4)
        allocate( diag_heat_pv_psi(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_heat_pv_psi),diag_1d_dim4', &
             diag_1d_dim4,size(diag_heat_pv_psi_host))
     endif

     if (allocated(diag_heat_pv)) then
        diag_heat_pv = 0_work_p
     endif

     if (allocated(diag_heat_pv_psi)) then
        diag_heat_pv_psi = 0_work_p
     endif

   endif

   if(diag_eflux_on_host) then

     need_diag_1d_eflux_pv = allocated(diag_1d_eflux_pv_host)
  
     if (need_diag_1d_eflux_pv) then
        if (allocated(diag_1d_eflux_pv)) then
           deallocate( diag_1d_eflux_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_1d_eflux_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_1d_eflux_pv_host,1)
        lb2 = lbound(diag_1d_eflux_pv_host,2)
        lb3 = lbound(diag_1d_eflux_pv_host,3)
        lb4 = lbound(diag_1d_eflux_pv_host,4)

        ub1 = ubound(diag_1d_eflux_pv_host,1)
        ub2 = ubound(diag_1d_eflux_pv_host,2)
        ub3 = ubound(diag_1d_eflux_pv_host,3)
        ub4 = ubound(diag_1d_eflux_pv_host,4)
        allocate( diag_1d_eflux_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_1d_eflux_pv),diag_1d_dim4', &
             diag_1d_dim4,size(diag_1d_eflux_pv_host))
     endif

     need_diag_2d_dflux_pv = allocated(diag_2d_dflux_pv_host)
     if (need_diag_2d_dflux_pv) then
        if (allocated(diag_2d_dflux_pv)) then
           deallocate( diag_2d_dflux_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_2d_dflux_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_2d_dflux_pv_host,1)
        lb2 = lbound(diag_2d_dflux_pv_host,2)
        lb3 = lbound(diag_2d_dflux_pv_host,3)
        lb4 = lbound(diag_2d_dflux_pv_host,4)

        ub1 = ubound(diag_2d_dflux_pv_host,1)
        ub2 = ubound(diag_2d_dflux_pv_host,2)
        ub3 = ubound(diag_2d_dflux_pv_host,3)
        ub4 = ubound(diag_2d_dflux_pv_host,4)
        allocate( diag_2d_dflux_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_2d_dflux_pv),diag_1d_dim4', &
             diag_1d_dim4,size(diag_2d_dflux_pv))
     endif

     if (allocated(diag_1d_eflux_pv)) then
        diag_1d_eflux_pv = 0_work_p
     endif

     if (allocated(diag_2d_dflux_pv)) then
        diag_2d_dflux_pv = 0_work_p
     endif

   endif ! diag_eflux_on_host

   return
   end subroutine update_device_diag



   attributes(host) &
   subroutine update_host_diag()
   use sml_module, only : sml_mype
   use precision_mod_gpu
   use diag_module, only : &
      diag_1d_npsi_host => diag_1d_npsi, &
      diag_1d_dp_inv_host => diag_1d_dp_inv, &
      diag_1d_pin_host => diag_1d_pin, &
      diag_1d_df_pv1_host => diag_1d_df_pv1, &
      diag_1d_f_pv1_host => diag_1d_f_pv1, &
      diag_heat_pv_host => diag_heat_pv, &
      diag_heat_pv_psi_host => diag_heat_pv_psi, &
      diag_1d_eflux_pv_host => diag_1d_eflux_pv, &
      diag_2d_dflux_pv_host => diag_2d_dflux_pv, &
      diag_eflux_on_host => diag_eflux_on, &
      diag_heat_on_host => diag_heat_on

   integer, parameter :: idebug = 2

   real (kind=work_p), allocatable, dimension(:,:,:,:) :: p_diag_1d_f_pv1
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: p_diag_1d_df_pv1

   real (kind=work_p), allocatable, dimension(:,:,:,:,:,:) :: p_diag_heat_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_heat_pv_psi

   real (kind=work_p), allocatable, dimension(:,:,:) :: dsum_f, dsum_df
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: dsum_heat_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_heat_psi
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_eflux_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_dflux_pv

   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_1d_eflux_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_2d_dflux_pv

   integer :: i1,i2,i3,i4,i5,i6
   integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4
   integer :: llb1,uub1,llb2,uub2,llb3,uub3,llb4,uub4

   integer :: lb1_df,ub1_df,lb2_df,ub2_df,lb3_df,ub3_df,lb4_df,ub4_df
   integer :: llb1_df,uub1_df,llb2_df,uub2_df,llb3_df,uub3_df,llb4_df,uub4_df

   integer :: lb1_pv,ub1_pv,lb2_pv,ub2_pv,lb3_pv,ub3_pv,lb4_pv,ub4_pv,lb5_pv,ub5_pv,lb6_pv,ub6_pv
   integer :: lb1_psi,ub1_psi,lb2_psi,ub2_psi,lb3_psi,ub3_psi,lb4_psi,ub4_psi,lb5_psi,ub5_psi

   integer :: lb1_ef,ub1_ef,lb2_ef,ub2_ef,lb3_ef,ub3_ef,lb4_ef,ub4_ef,lb5_ef,ub5_ef
   integer :: lb1_2df,ub1_2df,lb2_2df,ub2_2df,lb3_2df,ub3_2df,lb4_2df,ub4_2df,lb5_2df,ub5_2df

   integer :: ierr

   logical,parameter :: use_dsum = .true.
   logical,parameter :: use_p_diag = .true.
   logical :: need_diag_1d_df_pv1
   logical :: need_diag_heat_pv
   logical :: need_diag_heat_pv_psi
   logical :: need_diag_1d_eflux_pv
   logical :: need_diag_2d_dflux_pv
   integer :: n1,n2,n3,n4,n5,n6
   integer :: l1,u1,l2,u2,l3,u3,l4,u4,l5,u5,l6,u6

   diag_1d_npsi_host = diag_1d_npsi
   diag_1d_dp_inv_host = diag_1d_dp_inv
   diag_1d_pin_host = diag_1d_pin

   need_diag_1d_df_pv1 = allocated(diag_1d_df_pv1_host)
   need_diag_heat_pv = allocated(diag_heat_pv_host) 
   need_diag_heat_pv_psi = allocated(diag_heat_pv_psi_host)
   need_diag_1d_eflux_pv = allocated(diag_1d_eflux_pv_host)
   need_diag_2d_dflux_pv = allocated(diag_2d_dflux_pv_host)

   lb1 = lbound(diag_1d_f_pv1_host,1)
   lb2 = lbound(diag_1d_f_pv1_host,2)
   lb3 = lbound(diag_1d_f_pv1_host,3)
   lb4 = lbound(diag_1d_f_pv1_host,4)

   ub1 = ubound(diag_1d_f_pv1_host,1)
   ub2 = ubound(diag_1d_f_pv1_host,2)
   ub3 = ubound(diag_1d_f_pv1_host,3)
   ub4 = ubound(diag_1d_f_pv1_host,4)


   lb1_df = lbound(diag_1d_df_pv1_host,1)
   lb2_df = lbound(diag_1d_df_pv1_host,2)
   lb3_df = lbound(diag_1d_df_pv1_host,3)
   lb4_df = lbound(diag_1d_df_pv1_host,4)

   ub1_df = ubound(diag_1d_df_pv1_host,1)
   ub2_df = ubound(diag_1d_df_pv1_host,2)
   ub3_df = ubound(diag_1d_df_pv1_host,3)
   ub4_df = ubound(diag_1d_df_pv1_host,4)

   llb1 = lbound(diag_1d_f_pv1,1)
   llb2 = lbound(diag_1d_f_pv1,2)
   llb3 = lbound(diag_1d_f_pv1,3)
   llb4 = lbound(diag_1d_f_pv1,4)

   uub1 = ubound(diag_1d_f_pv1,1)
   uub2 = ubound(diag_1d_f_pv1,2)
   uub3 = ubound(diag_1d_f_pv1,3)
   uub4 = ubound(diag_1d_f_pv1,4)

   if (idebug >= 2) then
   call assert( llb1 <= lb1,'diag_module_gpu: invalid llb1 ',lb1,llb1)
   call assert( llb2 <= lb2,'diag_module_gpu: invalid llb2 ',lb2,llb2)
   call assert( llb3 <= lb3,'diag_module_gpu: invalid llb3 ',lb3,llb3)

   call assert( ub1 <= uub1,'diag_module_gpu: invalid uub1 ',ub1,uub1)
   call assert( ub2 <= uub2,'diag_module_gpu: invalid uub2 ',ub2,uub2)
   call assert( ub3 <= uub3,'diag_module_gpu: invalid uub3 ',ub3,uub3)
   endif

   if (use_p_diag) then

   allocate(p_diag_1d_f_pv1(llb1:uub1,llb2:uub2,llb3:uub3,llb4:uub4), &
            stat=ierr)
   if (idebug >= 1) then
     call assert(ierr.eq.0,'alloc(p_diag_1d_f_pv1) ',ierr)
   endif

   p_diag_1d_f_pv1 = diag_1d_f_pv1

   endif

   if (need_diag_1d_df_pv1) then
     llb1_df = lbound(diag_1d_df_pv1,1)
     llb2_df = lbound(diag_1d_df_pv1,2)
     llb3_df = lbound(diag_1d_df_pv1,3)
     llb4_df = lbound(diag_1d_df_pv1,4)

     uub1_df = ubound(diag_1d_df_pv1,1)
     uub2_df = ubound(diag_1d_df_pv1,2)
     uub3_df = ubound(diag_1d_df_pv1,3)
     uub4_df = ubound(diag_1d_df_pv1,4)

     if (use_p_diag) then

     allocate(p_diag_1d_df_pv1(llb1_df:uub1_df, &
                llb2_df:uub2_df,llb3_df:uub3_df,llb4_df:uub4_df),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_1d_df_pv1) ',ierr)
     endif
     p_diag_1d_df_pv1 = 0

     p_diag_1d_df_pv1 = diag_1d_df_pv1

     endif
   endif

  if(diag_heat_on_host) then 
   if (need_diag_heat_pv) then
     lb1_pv = lbound(diag_heat_pv,1)
     lb2_pv = lbound(diag_heat_pv,2)
     lb3_pv = lbound(diag_heat_pv,3)
     lb4_pv = lbound(diag_heat_pv,4)
     lb5_pv = lbound(diag_heat_pv,5)
     lb6_pv = lbound(diag_heat_pv,6)

     ub1_pv = ubound(diag_heat_pv,1)
     ub2_pv = ubound(diag_heat_pv,2)
     ub3_pv = ubound(diag_heat_pv,3)
     ub4_pv = ubound(diag_heat_pv,4)
     ub5_pv = ubound(diag_heat_pv,5)
     ub6_pv = ubound(diag_heat_pv,6)

     if (use_p_diag) then

     allocate(p_diag_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv, &
            lb4_pv:ub4_pv,lb5_pv:ub5_pv,lb6_pv:ub6_pv),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_heat_pv) ',ierr)
     endif
     p_diag_heat_pv = 0
     p_diag_heat_pv = diag_heat_pv
     endif
   endif

   if (need_diag_heat_pv_psi) then
     lb1_psi = lbound(diag_heat_pv_psi,1)
     lb2_psi = lbound(diag_heat_pv_psi,2)
     lb3_psi = lbound(diag_heat_pv_psi,3)
     lb4_psi = lbound(diag_heat_pv_psi,4)
     lb5_psi = lbound(diag_heat_pv_psi,5)

     ub1_psi = ubound(diag_heat_pv_psi,1)
     ub2_psi = ubound(diag_heat_pv_psi,2)
     ub3_psi = ubound(diag_heat_pv_psi,3)
     ub4_psi = ubound(diag_heat_pv_psi,4)
     ub5_psi = ubound(diag_heat_pv_psi,5)

     if (use_p_diag) then

     allocate(p_diag_heat_pv_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi, &
            lb3_psi:ub3_psi,lb4_psi:ub4_psi,lb5_psi:ub5_psi),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_heat_pv_psi) ',ierr)
     endif
     p_diag_heat_pv_psi = 0
     p_diag_heat_pv_psi = diag_heat_pv_psi

     endif
   endif
 endif ! diag_heat_on

  if(diag_eflux_on_host) then
   if (need_diag_1d_eflux_pv) then
     lb1_ef = lbound(diag_1d_eflux_pv,1)
     lb2_ef = lbound(diag_1d_eflux_pv,2)
     lb3_ef = lbound(diag_1d_eflux_pv,3)
     lb4_ef = lbound(diag_1d_eflux_pv,4)
     lb5_ef = lbound(diag_1d_eflux_pv,5)

     ub1_ef = ubound(diag_1d_eflux_pv,1)
     ub2_ef = ubound(diag_1d_eflux_pv,2)
     ub3_ef = ubound(diag_1d_eflux_pv,3)
     ub4_ef = ubound(diag_1d_eflux_pv,4)
     ub5_ef = ubound(diag_1d_eflux_pv,5)

     if (use_p_diag) then

     allocate(p_diag_1d_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef, &
            lb4_ef:ub4_ef,lb5_ef:ub5_ef),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_1d_eflux_pv) ',ierr)
     endif
     p_diag_1d_eflux_pv = 0
     p_diag_1d_eflux_pv = diag_1d_eflux_pv

     endif
   endif

   if (need_diag_2d_dflux_pv) then
     lb1_2df = lbound(diag_2d_dflux_pv,1)
     lb2_2df = lbound(diag_2d_dflux_pv,2)
     lb3_2df = lbound(diag_2d_dflux_pv,3)
     lb4_2df = lbound(diag_2d_dflux_pv,4)
     lb5_2df = lbound(diag_2d_dflux_pv,5)

     ub1_2df = ubound(diag_2d_dflux_pv,1)
     ub2_2df = ubound(diag_2d_dflux_pv,2)
     ub3_2df = ubound(diag_2d_dflux_pv,3)
     ub4_2df = ubound(diag_2d_dflux_pv,4)
     ub5_2df = ubound(diag_2d_dflux_pv,5)

     if (use_p_diag) then

     allocate(p_diag_2d_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df, &
            lb3_2df:ub3_2df,lb4_2df:ub4_2df,lb5_2df:ub5_2df),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_2d_dflux_pv) ',ierr)
     endif
     p_diag_2d_dflux_pv = 0
     p_diag_2d_dflux_pv = diag_2d_dflux_pv

     endif
   endif
 endif ! diag_eflux_on

! ------------------
! perform sum reduce
! into diag_1d_f_pv1(:,:,:,1) so
! that it is independent of sml_nthreads
! ------------------
   if (use_dsum) then
     allocate( dsum_f(llb1:uub1,llb2:uub2,llb3:uub3), stat=ierr)
     call assert(ierr.eq.0,'alloc(dsum_f) ',ierr)

     dsum_f = 0

     if (use_p_diag) then
     do i4=llb4,uub4
!$omp parallel do private(i1,i2,i3)
     do i3=lb3,ub3
     do i2=lb2,ub2
     do i1=lb1,ub1
       dsum_f(i1,i2,i3) = dsum_f(i1,i2,i3) + &
          p_diag_1d_f_pv1(i1,i2,i3,i4)
     enddo
     enddo
     enddo
     enddo

     else
       n1 = size(diag_1d_f_pv1,1)
       n2 = size(diag_1d_f_pv1,2)
       n3 = size(diag_1d_f_pv1,3)
       n4 = size(diag_1d_f_pv1,4)

       l1 = lbound(diag_1d_f_pv1,1)
       l2 = lbound(diag_1d_f_pv1,2)
       l3 = lbound(diag_1d_f_pv1,3)
       l4 = lbound(diag_1d_f_pv1,4)

       u1 = ubound(diag_1d_f_pv1,1)
       u2 = ubound(diag_1d_f_pv1,2)
       u3 = ubound(diag_1d_f_pv1,3)
       u4 = ubound(diag_1d_f_pv1,4)
       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_f_pv1, dsum_f )
       call sum2_gpu( n1*n2*n3, n4, diag_1d_f_pv1, dsum_f )
     endif

     i4 = 1
!$omp parallel do private(i1,i2,i3)
     do i3=lb3,ub3
     do i2=lb2,ub2
     do i1=lb1,ub1
       diag_1d_f_pv1_host(i1,i2,i3,1) = &
         diag_1d_f_pv1_host(i1,i2,i3,1) + dsum_f(i1,i2,i3)
     enddo
     enddo
     enddo

     if (need_diag_1d_df_pv1) then
       allocate( dsum_df(llb1:uub1,llb2:uub2,llb3:uub3), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_df) ',ierr)
       endif

       dsum_df = 0


       if (use_p_diag) then

       do i4=llb4_df,uub4_df
!$omp parallel do private(i1,i2,i3)
       do i3=lb3_df,ub3_df
       do i2=lb2_df,ub2_df
       do i1=lb1_df,ub1_df
        dsum_df(i1,i2,i3) = dsum_df(i1,i2,i3) + &
          p_diag_1d_df_pv1(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo

       else


       n1 = size(diag_1d_df_pv1,1)
       n2 = size(diag_1d_df_pv1,2)
       n3 = size(diag_1d_df_pv1,3)
       n4 = size(diag_1d_df_pv1,4)

       l1 = lbound(diag_1d_df_pv1,1)
       l2 = lbound(diag_1d_df_pv1,2)
       l3 = lbound(diag_1d_df_pv1,3)
       l4 = lbound(diag_1d_df_pv1,4)

       u1 = ubound(diag_1d_df_pv1,1)
       u2 = ubound(diag_1d_df_pv1,2)
       u3 = ubound(diag_1d_df_pv1,3)
       u4 = ubound(diag_1d_df_pv1,4)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3, n4, diag_1d_df_pv1, dsum_df )

       endif


       i4 = 1
!$omp parallel do private(i1,i2,i3)
       do i3=lb3_df,ub3_df
       do i2=lb2_df,ub2_df
       do i1=lb1_df,ub1_df
         diag_1d_df_pv1_host(i1,i2,i3,1) = &
           diag_1d_df_pv1_host(i1,i2,i3,1) + dsum_df(i1,i2,i3)
       enddo
       enddo
       enddo
     endif

    if(diag_heat_on_host) then
     if (need_diag_heat_pv) then
       allocate( dsum_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv, &
                 lb4_pv:ub4_pv,lb5_pv:ub5_pv), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_heat_pv) ',ierr)
       endif

       dsum_heat_pv = 0


       if (use_p_diag) then

       do i6=lb6_pv,ub6_pv
!$omp parallel do private(i1,i2,i3,i4,i5)
       do i5=lb5_pv,ub5_pv
       do i4=lb4_pv,ub4_pv
       do i3=lb3_pv,ub3_pv
       do i2=lb2_pv,ub2_pv
       do i1=lb1_pv,ub1_pv
        dsum_heat_pv(i1,i2,i3,i4,i5) = dsum_heat_pv(i1,i2,i3,i4,i5) + &
          p_diag_heat_pv(i1,i2,i3,i4,i5,i6)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1 = size(diag_heat_pv,1)
       n2 = size(diag_heat_pv,2)
       n3 = size(diag_heat_pv,3)
       n4 = size(diag_heat_pv,4)
       n5 = size(diag_heat_pv,5)
       n6 = size(diag_heat_pv,6)

       l1 = lbound(diag_heat_pv,1)
       l2 = lbound(diag_heat_pv,2)
       l3 = lbound(diag_heat_pv,3)
       l4 = lbound(diag_heat_pv,4)
       l5 = lbound(diag_heat_pv,5)
       l6 = lbound(diag_heat_pv,6)

       u1 = ubound(diag_heat_pv,1)
       u2 = ubound(diag_heat_pv,2)
       u3 = ubound(diag_heat_pv,3)
       u4 = ubound(diag_heat_pv,4)
       u5 = ubound(diag_heat_pv,5)
       u6 = ubound(diag_heat_pv,6)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4*n5, n6, diag_heat_pv, dsum_heat_pv )

       endif

       i6 = 1
!$omp parallel do private(i1,i2,i3,i4,i5)
       do i5=lb5_pv,ub5_pv
       do i4=lb4_pv,ub4_pv
       do i3=lb3_pv,ub3_pv
       do i2=lb2_pv,ub2_pv
       do i1=lb1_pv,ub1_pv
         diag_heat_pv_host(i1,i2,i3,i4,i5,1) = &
           diag_heat_pv_host(i1,i2,i3,i4,i5,1) + dsum_heat_pv(i1,i2,i3,i4,i5)
       !if(diag_heat_pv_host(i1,i2,i3,i4,0,1)/=0) print *, "nonzero diag_heat_pv_host(i1,i2,i3,i4,0,1)"
       enddo
       enddo
       enddo
       enddo
       enddo
     endif


     if (need_diag_heat_pv_psi) then
       allocate( dsum_heat_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi, &
                 lb4_psi:ub4_psi), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_heat_psi) ',ierr)
       endif

       dsum_heat_psi = 0

       if (use_p_diag) then
       do i5=lb5_psi,ub5_psi
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_psi,ub4_psi
       do i3=lb3_psi,ub3_psi
       do i2=lb2_psi,ub2_psi
       do i1=lb1_psi,ub1_psi
        dsum_heat_psi(i1,i2,i3,i4) = dsum_heat_psi(i1,i2,i3,i4) + &
          p_diag_heat_pv_psi(i1,i2,i3,i4,i5)
          !if(p_diag_heat_pv_psi(i1,i2,i3,0,i5)/=0) print *, "nonzero p_diag_heat_pv_psi(i1,i2,i3,i4,i5)"
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1 = size(diag_heat_pv_psi,1)
       n2 = size(diag_heat_pv_psi,2)
       n3 = size(diag_heat_pv_psi,3)
       n4 = size(diag_heat_pv_psi,4)
       n5 = size(diag_heat_pv_psi,5)

       l1 = lbound(diag_heat_pv_psi,1)
       l2 = lbound(diag_heat_pv_psi,2)
       l3 = lbound(diag_heat_pv_psi,3)
       l4 = lbound(diag_heat_pv_psi,4)
       l5 = lbound(diag_heat_pv_psi,5)

       u1 = ubound(diag_heat_pv_psi,1)
       u2 = ubound(diag_heat_pv_psi,2)
       u3 = ubound(diag_heat_pv_psi,3)
       u4 = ubound(diag_heat_pv_psi,4)
       u5 = ubound(diag_heat_pv_psi,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_heat_pv_psi, dsum_heat_psi )

       endif

       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_psi,ub4_psi
       do i3=lb3_psi,ub3_psi
       do i2=lb2_psi,ub2_psi
       do i1=lb1_psi,ub1_psi
         diag_heat_pv_psi_host(i1,i2,i3,i4,1) = &
           diag_heat_pv_psi_host(i1,i2,i3,i4,1) + dsum_heat_psi(i1,i2,i3,i4)
           !if(abs(diag_heat_pv_psi_host(i1,i2,i3,0,1))>0) print *,"nonzero diag_heat_psi_cpu"
       enddo
       enddo
       enddo
       enddo
     endif
    endif ! diag_heat_on

    if(diag_eflux_on_host) then
     if (need_diag_1d_eflux_pv) then
       allocate( dsum_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef, &
                 lb4_ef:ub4_ef), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_eflux_pv) ',ierr)
       endif

       dsum_eflux_pv = 0

       if (use_p_diag) then

       do i5=lb5_ef,ub5_ef
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_ef,ub4_ef
       do i3=lb3_ef,ub3_ef
       do i2=lb2_ef,ub2_ef
       do i1=lb1_ef,ub1_ef
        dsum_eflux_pv(i1,i2,i3,i4) = dsum_eflux_pv(i1,i2,i3,i4) + &
          p_diag_1d_eflux_pv(i1,i2,i3,i4,i5)
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1 = size(diag_1d_eflux_pv,1)
       n2 = size(diag_1d_eflux_pv,2)
       n3 = size(diag_1d_eflux_pv,3)
       n4 = size(diag_1d_eflux_pv,4)
       n5 = size(diag_1d_eflux_pv,5)

       l1 = lbound(diag_1d_eflux_pv,1)
       l2 = lbound(diag_1d_eflux_pv,2)
       l3 = lbound(diag_1d_eflux_pv,3)
       l4 = lbound(diag_1d_eflux_pv,4)
       l5 = lbound(diag_1d_eflux_pv,5)

       u1 = ubound(diag_1d_eflux_pv,1)
       u2 = ubound(diag_1d_eflux_pv,2)
       u3 = ubound(diag_1d_eflux_pv,3)
       u4 = ubound(diag_1d_eflux_pv,4)
       u5 = ubound(diag_1d_eflux_pv,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_1d_eflux_pv, dsum_eflux_pv )

       endif
       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_ef,ub4_ef
       do i3=lb3_ef,ub3_ef
       do i2=lb2_ef,ub2_ef
       do i1=lb1_ef,ub1_ef
         diag_1d_eflux_pv_host(i1,i2,i3,i4,1) = &
           diag_1d_eflux_pv_host(i1,i2,i3,i4,1) + dsum_eflux_pv(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo
     endif

     if (need_diag_2d_dflux_pv) then
       allocate( dsum_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df, &
                 lb4_2df:ub4_2df), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_dflux_pv) ',ierr)
       endif

       dsum_dflux_pv = 0

       if (use_p_diag) then

       do i5=lb5_2df,ub5_2df
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_2df,ub4_2df
       do i3=lb3_2df,ub3_2df
       do i2=lb2_2df,ub2_2df
       do i1=lb1_2df,ub1_2df
        dsum_dflux_pv(i1,i2,i3,i4) = dsum_dflux_pv(i1,i2,i3,i4) + &
          p_diag_2d_dflux_pv(i1,i2,i3,i4,i5)
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1 = size(diag_2d_dflux_pv,1)
       n2 = size(diag_2d_dflux_pv,2)
       n3 = size(diag_2d_dflux_pv,3)
       n4 = size(diag_2d_dflux_pv,4)
       n5 = size(diag_2d_dflux_pv,5)

       l1 = lbound(diag_2d_dflux_pv,1)
       l2 = lbound(diag_2d_dflux_pv,2)
       l3 = lbound(diag_2d_dflux_pv,3)
       l4 = lbound(diag_2d_dflux_pv,4)
       l5 = lbound(diag_2d_dflux_pv,5)

       u1 = ubound(diag_2d_dflux_pv,1)
       u2 = ubound(diag_2d_dflux_pv,2)
       u3 = ubound(diag_2d_dflux_pv,3)
       u4 = ubound(diag_2d_dflux_pv,4)
       u5 = ubound(diag_2d_dflux_pv,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_2d_dflux_pv, dsum_dflux_pv )

       endif
       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_2df,ub4_2df
       do i3=lb3_2df,ub3_2df
       do i2=lb2_2df,ub2_2df
       do i1=lb1_2df,ub1_2df
         diag_2d_dflux_pv_host(i1,i2,i3,i4,1) = &
           diag_2d_dflux_pv_host(i1,i2,i3,i4,1) + dsum_dflux_pv(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo
     endif
    endif ! diag_eflux_on

     if (allocated(dsum_f)) then
       deallocate( dsum_f, stat=ierr )
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_f) ',ierr)
       endif
     endif

     if (allocated(dsum_df)) then
       deallocate( dsum_df, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_df) ',ierr)
       endif
     endif

   if(diag_heat_on_host) then
     if (allocated(dsum_heat_pv)) then
       deallocate( dsum_heat_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_heat_pv) ',ierr)
       endif
     endif

     if (allocated(dsum_heat_psi)) then
       deallocate( dsum_heat_psi, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_heat_psi) ',ierr)
       endif
     endif
   endif ! diag_heat_on

   if(diag_eflux_on_host) then
     if (allocated(dsum_eflux_pv)) then
       deallocate( dsum_eflux_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_eflux_pv) ',ierr)
       endif
     endif

     if (allocated(dsum_dflux_pv)) then
       deallocate( dsum_dflux_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_dflux_psi) ',ierr)
       endif
     endif
   endif ! diag_eflux_on

   else
   do i4=llb4,uub4
     diag_1d_f_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) = &
       diag_1d_f_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) + &
         p_diag_1d_f_pv1(lb1:ub1,lb2:ub2,lb3:ub3, i4 )
   enddo

    if (need_diag_1d_df_pv1) then
     do i4=llb4_df,uub4_df
      diag_1d_df_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) = &
       diag_1d_df_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) + &
         p_diag_1d_df_pv1(lb1:ub1,lb2:ub2,lb3:ub3, i4 )
     enddo
    endif

   if(diag_heat_on_host) then
    if (need_diag_heat_pv) then
     do i6=lb6_pv,ub6_pv
      diag_heat_pv_host(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,1) = &
       diag_heat_pv_host(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,1) + &
         p_diag_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,i6 )
     enddo
    endif

    if (need_diag_heat_pv_psi) then
     do i5=lb5_psi,ub5_psi
      diag_heat_pv_psi_host(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,1) = &
       diag_heat_pv_psi_host(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,1) + &
       p_diag_heat_pv_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,i5 )
     enddo
    endif
   endif ! diag_heat_on

   if(diag_eflux_on_host) then
    if (need_diag_1d_eflux_pv) then
     do i5=lb5_ef,ub5_ef
      diag_1d_eflux_pv_host(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,1) = &
       diag_1d_eflux_pv_host(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,1) + &
         p_diag_1d_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,i5 )
     enddo
    endif

    if (need_diag_2d_dflux_pv) then
     do i5=lb5_2df,ub5_2df
      diag_2d_dflux_pv_host(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,1) = &
       diag_2d_dflux_pv_host(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,1) + &
       p_diag_2d_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,i5 )
     enddo
    endif

   endif ! diag_eflux_on

   endif

! -------
! clean up
! -------

   if (allocated(p_diag_1d_f_pv1)) then
     deallocate( p_diag_1d_f_pv1, stat=ierr )
     call assert( ierr.eq.0,'dealloc(p_diag_1d_f_pv1 ',ierr)
   endif


   if (allocated(p_diag_1d_df_pv1)) then
     deallocate( p_diag_1d_df_pv1, stat=ierr )
     call assert( ierr.eq.0,'dealloc(p_diag_1d_df_pv1 ',ierr)
   endif

   if (allocated(diag_1d_f_pv1)) then
     deallocate(diag_1d_f_pv1,stat=ierr)
     call assert(ierr.eq.0,'final dealloc(diag_1d_f_pv1),ierr=',ierr)
   endif

   if (allocated(diag_1d_df_pv1)) then
     deallocate( diag_1d_df_pv1,stat=ierr)
     call assert(ierr.eq.0,'final dealloc(diag_1d_df_pv1),ierr=',ierr)
   endif

   if(diag_eflux_on_host)then
     if (allocated(diag_1d_eflux_pv)) then
       deallocate(diag_1d_eflux_pv,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(diag_1d_eflux_pv),ierr=',ierr)
     endif
     if (allocated(diag_2d_dflux_pv)) then
       deallocate(diag_2d_dflux_pv,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(diag_2d_dflux_pv),ierr=',ierr)
     endif
   endif ! diag_eflux_on

   if(diag_heat_on_host)then
     if (allocated(diag_heat_pv)) then
        deallocate(diag_heat_pv,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(diag_heat_pv),ierr=',ierr)
     endif
     if (allocated(diag_heat_pv_psi)) then
        deallocate(diag_heat_pv_psi,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(diag_heat_pv_psi),ierr=',ierr)
     endif
     if (allocated(p_diag_heat_pv)) then
        deallocate(p_diag_heat_pv,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(p_diag_heat_pv),ierr=',ierr)
     endif
     if (allocated(p_diag_heat_pv_psi)) then
        deallocate(p_diag_heat_pv_psi,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(p_diag_heat_pv_psi),ierr=',ierr)
     endif

   endif ! diag_eflux_on


   return
   end subroutine update_host_diag
  
end module diag_module_gpu
module bnc_module_gpu
  use bnc_module, only : bnc_nr_host => bnc_nr
  use precision_mod_gpu
  integer, parameter :: bnc_nr=bnc_nr_host
  real (kind=work_p) :: bnc_min_r, bnc_max_r 
  real (kind=work_p) :: bnc_dr
  real (kind=work_p) :: bnc_z_psi_min(bnc_nr)
  ! real (8) :: bnc_arg_pass_r

  attributes(device) :: bnc_min_r, bnc_max_r
  attributes(device) :: bnc_dr, bnc_z_psi_min

  contains

  attributes(host) &
  subroutine update_device_bnc()
  use bnc_module, only : &
    bnc_z_psi_min_host => bnc_z_psi_min, &
    bnc_dr_host => bnc_dr, &
    bnc_min_r_host => bnc_min_r, &
    bnc_max_r_host => bnc_max_r

    bnc_z_psi_min = bnc_z_psi_min_host

    bnc_dr = bnc_dr_host
    bnc_min_r = bnc_min_r_host
    bnc_max_r = bnc_max_r_host

  end subroutine update_device_bnc

end module bnc_module_gpu

module push_mod_gpu
use cudafor
use util_mod_gpu

use dimensions_mod_gpu
use sml_module_gpu
use neu_module_gpu
use fld_module
use eq_module_gpu
use itp_module_gpu
use one_d_cub_mod_gpu
use bicub_mod_gpu
use ptl_module_gpu
use grid_class_gpu
use boundary_class_gpu
use psn_class_gpu
use diag_module_gpu
use bnc_module_gpu
use precision_mod_gpu
implicit none



contains



! -------------------------------------------------
! all code must be in the same module, a limitation
! of PGI cuda compiler
! -------------------------------------------------

!cdir$r unroll 
  attributes(device) &
  subroutine eval_bicub_2_loop(x,y,xc,yc,acoef,f00,f10,f01,f11,f20,f02)
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fx, fy
    real (kind=work_p), dimension(0:ndeg) :: dfx, dfy
    real (kind=work_p), dimension(0:ndeg) :: dfx2, dfy2
    logical, parameter :: use_fy = .false.
    f00 = 0.0d0
    f01 = 0.0d0
    f10 = 0.0d0
    f11 = 0.0d0
    f20 = 0.0d0
    f02 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy
    if (use_fy) then
       fy(0:ndeg) = 0.0d0
       dfy(0:ndeg) = 0.0d0
       dfy2(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j)*yv(j)
          enddo
          do j=1,ndeg
             dfy(i) = dfy(i) + acoef(i,j)*dble(j)*yv(j-1)
          enddo
          do j=2,ndeg
             dfy2(i) = dfy2(i) + acoef(i,j)*dble(j*(j-1))*yv(j-2)
          enddo
       enddo
       ! ----------------
       ! f00 = f(x,y)
       ! f10 = df/dx
       ! f01 = df/dy
       ! f11 = df^2/dx/dy
       ! f20 = df^2/dx/dx
       ! f02 = df^2/dy/dy
       ! ----------------
       do i=0,ndeg
          f00 = f00 + xv(i)*fy(i)
          f01 = f01 + xv(i)*dfy(i)
          f02 = f02 + xv(i)*dfy2(i)
       enddo
       do i=1,ndeg
          dfx(i) = dble(i)*xv(i-1)

          f10 = f10 + dfx(i)*fy(i)
          f11 = f11 + dfx(i)*dfy(i)
       enddo

       do i=2,ndeg
          dfx2(i) = dble(i*(i-1))*xv(i-2)

          f20 = f20 + dfx2(i)*fy(i)
       enddo


    else

       fx(0:ndeg) = 0.0d0
       dfx(0:ndeg) = 0.0d0
       dfx2(0:ndeg) = 0.0d0


       do j=0,ndeg
          do i=0,ndeg
             fx(j) = fx(j) + xv(i)*acoef(i,j)
          enddo
          do i=1,ndeg
             dfx(j) = dfx(j) + dble(i)*xv(i-1)*acoef(i,j)
          enddo
          do i=2,ndeg
             dfx2(j) = dfx2(j) + dble(i*(i-1))*xv(i-2)*acoef(i,j)
          enddo
       enddo

       do j=0,ndeg
          f00 = f00 + fx(j)*yv(j)
          f10 = f10 + dfx(j)*yv(j)
          f20 = f20 + dfx2(j)*yv(j)
       enddo



       do j=1,ndeg
          dfy(j) = dble(j)*yv(j-1)

          f01 = f01 + fx(j)*dfy(j)
          f11 = f11 + dfx(j)*dfy(j)
       enddo

       do j=2,ndeg
          dfy2(j) = dble(j*(j-1))*yv(j-2)

          f02 = f02 + fx(j)*dfy2(j)
       enddo




    endif




    return
  end subroutine eval_bicub_2_loop

!cdir$r unroll 
  attributes(device) &
  subroutine eval_bicub_1_loop(x,y,xc,yc,acoef,f00,f10,f01)
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fx, fy
    real (kind=work_p), dimension(0:ndeg) :: dfx, dfy
    logical, parameter :: use_fy = .false.
    f00 = 0.0d0
    f01 = 0.0d0
    f10 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy
    if (use_fy) then
       fy(0:ndeg) = 0.0d0
       dfy(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j)*yv(j)
          enddo
          do j=1,ndeg
             dfy(i) = dfy(i) + acoef(i,j)*dble(j)*yv(j-1)
          enddo
       enddo
       ! ----------------
       ! f00 = f(x,y)
       ! f10 = df/dx
       ! f01 = df/dy
       ! f11 = df^2/dx/dy
       ! f20 = df^2/dx/dx
       ! f02 = df^2/dy/dy
       ! ----------------
       do i=0,ndeg
          f00 = f00 + xv(i)*fy(i)
          f01 = f01 + xv(i)*dfy(i)
       enddo
       do i=1,ndeg
          dfx(i) = dble(i)*xv(i-1)
          f10 = f10 + dfx(i)*fy(i)
       enddo
    else
       fx(0:ndeg) = 0.0d0
       dfx(0:ndeg) = 0.0d0
       do j=0,ndeg
          do i=0,ndeg
             fx(j) = fx(j) + xv(i)*acoef(i,j)
          enddo
          do i=1,ndeg
             dfx(j) = dfx(j) + dble(i)*xv(i-1)*acoef(i,j)
          enddo
       enddo
       do j=0,ndeg
          f00 = f00 + fx(j)*yv(j)
          f10 = f10 + dfx(j)*yv(j)
       enddo
       do j=1,ndeg
          dfy(j) = dble(j)*yv(j-1)
          f01 = f01 + fx(j)*dfy(j)
       enddo
    endif



    return
  end subroutine eval_bicub_1_loop

!cdir$r unroll 
  attributes(device) &
  subroutine eval_bicub_00_loop(x,y,xc,yc,acoef,f00 )
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv, fy
    real (kind=work_p) :: dx, dy
    dx = (x-xc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    dy = (y-yc)
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy

       fy(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j) * yv(j)
          enddo
       enddo
       f00 = 0.0d0
       do i=0,ndeg
          f00 = f00 + xv(i) * fy(i)
       enddo
    return
  end subroutine eval_bicub_00_loop



!cdir$r unroll 
  attributes(device) &
  subroutine bicub_interpol2(r,z, f00,f10,f01,f11,f20,f02)
    use bicub_mod_gpu
    use precision_mod_gpu

    implicit none
    real (kind=work_p), intent(in) :: r,z
    real (kind=work_p), intent(inout) :: f00,f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: abserr, g00,g10,g01,g11,g20,g02
    integer :: ierr
    logical :: isok
    integer :: nr,nz, i,j
    !real (8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = ubound(rc_cub_gpu,1)
    nz = ubound(zc_cub_gpu,1)
    i = max(1,min(nr, 1 + int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1 + int( (z-zmin_gpu)*dz_inv_gpu ) ))
    !psi_acoef_local(0:ndeg,0:ndeg) = psi_acoef(0:ndeg,0:ndeg,i,j)
    call eval_bicub_2(r,z, rc_cub_gpu(i), zc_cub_gpu(j), acoef_cub_gpu(0,0,i,j), &
         & f00,f10,f01,f11,f20,f02 )
    return
  end subroutine bicub_interpol2


!cdir$r unroll 
   attributes(device) &
  subroutine bicub_interpol1(r,z, f00,f10,f01)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    real (kind=work_p), intent(in) :: r,z
    real (kind=work_p), intent(inout) :: f00,f10,f01
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: abserr, g00,g10,g01
    integer :: ierr
    logical :: isok
    integer :: nr,nz, i,j
    !real (8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = ubound(rc_cub_gpu,1)
    nz = ubound(zc_cub_gpu,1)
    i = max(1,min(nr, 1 + int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1 + int( (z-zmin_gpu)*dz_inv_gpu ) ))
    !psi_acoef_local(0:ndeg,0:ndeg) = psi_acoef(0:ndeg,0:ndeg,i,j)
    call eval_bicub_1(r,z, rc_cub_gpu(i), zc_cub_gpu(j), acoef_cub_gpu(0,0,i,j), &
         & f00,f10,f01)
    return
  end subroutine bicub_interpol1

!cdir$r unroll 
  attributes(device) &
  real (8) function psi_interpol_00(r,z )
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    real (kind=work_p), intent(in) :: r,z
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: g00, abserr
    logical :: isok
    real (kind=work_p) :: f00
    integer :: nr,nz, i,j

    !real(8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = ubound(rc_cub_gpu,1)
    nz = ubound(zc_cub_gpu,1)
    i = max(1,min(nr, 1+int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1+int( (z-zmin_gpu)*dz_inv_gpu ) ))

    !psi_acoef_local(0:ndeg,0:ndeg) = psi_acoef(0:ndeg,0:ndeg,i,j)
    call eval_bicub_00(r,z, rc_cub_gpu(i), zc_cub_gpu(j), &
         & acoef_cub_gpu(0,0,i,j), f00)
    psi_interpol_00 = f00
    return
  end function psi_interpol_00




!cdir$r unroll 
  attributes(device) &
real (kind=8) function psi_interpol_gpu(r_in,z_in,r_der,z_der)
  use itp_module_gpu
  use bicub_mod_gpu
  use precision_mod_gpu

  implicit none
  real (kind=work_p), intent(in) :: r_in, z_in
  integer , intent(in) :: r_der, z_der
! real (kind=8) :: r,z
! real (kind=8) ,external:: psi_interpol_pspline
  real (kind=work_p) :: dpsi_dr, dpsi_dz, psi


  call bicub_interpol1(r_in,z_in,psi,dpsi_dr,dpsi_dz)
  if(r_der==0 .and. z_der==0) then
     psi_interpol_gpu=psi
  elseif(r_der==1 .and. z_der==0) then
     psi_interpol_gpu=dpsi_dr
  elseif(r_der==0 .and. z_der==1) then
     psi_interpol_gpu=dpsi_dz
  else
     !print *, 'Do not use second derivative of psi_interpol with USE_BICUB_MOD.'
     !print *, 'Replace it with psi_interpol_bicub(r,z,psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_dr2,d2psi_dz2).'
     !stop 
  endif

  if (r_der==0 .and. z_der==0) then
     psi_interpol_gpu = max(1D-99, psi_interpol_gpu)
  endif
  
end function psi_interpol_gpu






!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_2_unroll(x,y,xc,yc,A,f00,f10,f01,f11,f20,f02)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
! real (kind=work_p), dimension(0:ndeg) :: fx 
! real (kind=work_p), dimension(0:ndeg) :: dfx
! real (kind=work_p), dimension(0:ndeg) :: dfx2

   real(kind=work_p), dimension(0:ndeg) :: fx,dfx,dfx2



! f00 = 0.0d0
! f01 = 0.0d0
! f10 = 0.0d0
! f11 = 0.0d0
! f20 = 0.0d0
! f02 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy


! fx(0:ndeg) = 0.0d0
! dfx(0:ndeg) = 0.0d0
! dfx2(0:ndeg) = 0.0d0


! do j=0,ndeg
! do i=0,ndeg
! fx(j) = fx(j) + xv(i)*A(i,j)
! enddo
! enddo

        j = 0
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 1
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 2
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 3
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)

!! dfy(j) = dble(j)*yv(j-1)
!! dfy2(j) = dble(j*(j-1))*yv(j-2)
! f01 = f01 + fx(j)*dfy(j), j=1:ndeg
! f02 = f02 + fx(j)*dfy2(j), j=2:ndeg

       f00 = fx(0) + dy*fx(1) + (dy*dy)*fx(2) + ((dy*dy)*dy)*fx(3)
       f01 = fx(1) + 2*dy*fx(2) + 3*(dy*dy)*fx(3)
       f02 = 2*fx(2) + 3*2*dy*fx(3)


! do j=0,ndeg
! do i=1,ndeg
! dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
! enddo
! enddo

        j = 0
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 1
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 2
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 3
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)



! do j=1,ndeg
!! dfy(j) = dble(j)*yv(j-1)
!
! f11 = f11 + dfx(j)*dfy(j)
! enddo


! f11 = f11 + dfx(j)*dfy(j), j=1:ndeg
! f10 = f10 + dfx(j)*yv(j), j=0:ndeg
       f10 = dfx(0) + dy*dfx(1) + (dy*dy)*dfx(2) + ((dy*dy)*dy)*dfx(3)
       f11 = dfx(1) + 2*dy*dfx(2) + 3*dy*dy*dfx(3)


! do j=0,ndeg
! do i=2,ndeg
! dfx2(j) = dfx2(j) + dble(i*(i-1))*xv(i-2)*A(i,j)
! enddo
! enddo

       j = 0
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 1
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 2
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 3
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))

! do j=0,ndeg
! f20 = f20 + dfx2(j)*yv(j)
! enddo

       f20 = dfx2(0) + dy*dfx2(1) + (dy*dy)*dfx2(2) + ((dy*dy)*dy)*dfx2(3)


! do j=1,ndeg
!! dfy(j) = dble(j)*yv(j-1)
!
! f01 = f01 + fx(j)*dfy(j)
! f11 = f11 + dfx(j)*dfy(j)
! enddo

! do j=2,ndeg
!! dfy2(j) = dble(j*(j-1))*yv(j-2)
!
! f02 = f02 + fx(j)*dfy2(j)
! enddo



    return
  end subroutine eval_bicub_2_unroll






!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_1_unroll(x,y,xc,yc,A,f00,f10,f01)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
! real (kind=work_p), dimension(0:ndeg) :: fx, fy
! real (kind=work_p), dimension(0:ndeg) :: dfx, dfy

    real (kind=work_p), dimension(0:ndeg) :: fx,dfx
 
! f00 = 0.0d0
! f01 = 0.0d0
! f10 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy

! fx(0:ndeg) = 0.0d0
! dfx(0:ndeg) = 0.0d0
! do j=0,ndeg
! do i=0,ndeg
! fx(j) = fx(j) + xv(i)*A(i,j)
! enddo
! do i=1,ndeg
! dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
! enddo
! enddo


       j = 0
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 1
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 2
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 3
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)

       f00 = fx(0) + dy*fx(1) + (dy*dy)*fx(2) + ((dy*dy)*dy)*fx(3)
       f01 = fx(1) + 2*dy*fx(2) + 3*(dy*dy)*fx(3)

       j = 0
       dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
       j = 1
       dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
       j = 2
       dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
       j = 3
       dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)

       f10 = dfx(0) + dy*dfx(1) + (dy*dy)*dfx(2) + ((dy*dy)*dy)*dfx(3)

! do j=0,ndeg
! f00 = f00 + fx(j)*yv(j)
! f10 = f10 + dfx(j)*yv(j)
! enddo
!
! do j=1,ndeg
! dfy(j) = dble(j)*yv(j-1)
! f01 = f01 + fx(j)*dfy(j)
! enddo

    return
  end subroutine eval_bicub_1_unroll




!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_00_unroll(x,y,xc,yc,A,f00 )
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00
    integer :: i
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fy
    real (kind=work_p) :: dx, dy
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy

! fy(0:ndeg) = 0.0d0
! do i=0,ndeg
! do j=0,ndeg
! fy(i) = fy(i) + A(i,j) * yv(j)
! enddo
! enddo

        i = 0
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 1
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 2
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 3
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)

! f00 = 0.0d0
! do i=0,ndeg
! f00 = f00 + xv(i) * fy(i)
! enddo

    f00 = fy(0) + dx*fy(1) + (dx*dx)*fy(2) + ((dx*dx)*dx)*fy(3)
    return
  end subroutine eval_bicub_00_unroll



!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_00(x,y,xc,yc,A,f00_ )
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(out) :: f00_

    real (kind=work_p) :: f00

    integer :: i
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p) :: dx, dy
    real (kind=work_p) :: fy_i
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy

! fy(0:ndeg) = 0.0d0
! do i=0,ndeg
! do j=0,ndeg
! fy(i) = fy(i) + A(i,j) * yv(j)
! enddo
! enddo

        f00 = 0

! fy_0 = A(0,0) + A(0,1)*dy + A(0,2)*(dy*dy) + A(0,3)*((dy*dy)*dy)

        fy_i = ((A(0,3)*dy + A(0,2))*dy + A(0,1))*dy + A(0,0)
        f00 = f00 + fy_i

! fy_1 = A(1,0) + A(1,1)*dy + A(1,2)*(dy*dy) + A(1,3)*((dy*dy)*dy)
        fy_i = ((A(1,3)*dy + A(1,2))*dy + A(1,1))*dy + A(1,0)
        f00 = f00 + dx*fy_i
          

! fy_2 = A(2,0) + A(2,1)*dy + A(2,2)*(dy*dy) + A(2,3)*((dy*dy)*dy)
        fy_i = ((A(2,3)*dy + A(2,2))*dy + A(2,1))*dy + A(2,0)
        f00 = f00 + (dx*dx)*fy_i

! fy_i = A(3,0) + A(3,1)*dy + A(3,2)*(dy*dy) + A(3,3)*((dy*dy)*dy)
        fy_i = ((A(3,3)*dy + A(3,2))*dy + A(3,1))*dy + A(3,0)
        f00 = f00 + (dx*dx)*dx*fy_i

! f00 = 0.0d0
! do i=0,ndeg
! f00 = f00 + xv(i) * fy(i)
! enddo
! f00 = fy_0 + dx*fy_1 + (dx*dx)*fy_2 + ((dx*dx)*dx)*fy_3

    f00_ = f00
    return
  end subroutine eval_bicub_00



!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_1(x,y,xc,yc,A,f00_,f10_,f01_)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(out) :: f00_, f10_,f01_
    real (kind=work_p) :: f00, f10,f01
    real (kind=work_p) :: dx,dy
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
! real (kind=work_p), dimension(0:ndeg) :: fx, fy
! real (kind=work_p), dimension(0:ndeg) :: dfx, dfy

    real (kind=work_p) :: dfx_i
    real (kind=work_p) :: fx_i
 
! f00 = 0.0d0
! f01 = 0.0d0
! f10 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy

! fx(0:ndeg) = 0.0d0
! dfx(0:ndeg) = 0.0d0
! do j=0,ndeg
! do i=0,ndeg
! fx(j) = fx(j) + xv(i)*A(i,j)
! enddo
! do i=1,ndeg
! dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
! enddo
! enddo


       f00 = 0
       f01 = 0
! fx_0 = A(0,0) + dx*A(1,0) + (dx*dx)*A(2,0) + ((dx*dx)*dx)*A(3,0)
       fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0)
       f00 = f00 + fx_i

! fx_1 = A(0,1) + dx*A(1,1) + (dx*dx)*A(2,1) + ((dx*dx)*dx)*A(3,1)
       fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1)
       f00 = f00 + dy*fx_i
       f01 = f01 + fx_i


! fx_2 = A(0,2) + dx*A(1,2) + (dx*dx)*A(2,2) + ((dx*dx)*dx)*A(3,2)
       fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2)
       f00 = f00 + (dy*dy)*fx_i
       f01 = f01 + 2.0d0*dy*fx_i


! fx_3 = A(0,3) + dx*A(1,3) + (dx*dx)*A(2,3) + ((dx*dx)*dx)*A(3,3)
       fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3)
       f00 = f00 + dy*((dy*dy)*fx_i)
       f01 = f01 + 3.0d0*((dy*dy)*fx_i)

! f00 = fx_0 + dy*fx_1 + (dy*dy)*fx_2 + ((dy*dy)*dy)*fx_3
! f01 = fx_1 + 2*dy*fx_2 + 3*(dy*dy)*fx_3

       f10 = 0

! dfx_0 = A(1,0) + 2*dx*A(2,0) + 3*(dx*dx)*A(3,0)

       dfx_i = (A(3,0)*3.0d0*dx + A(2,0)*2.0d0)*dx + A(1,0)
       f10 = f10 + dfx_i

! dfx_1 = A(1,1) + 2*dx*A(2,1) + 3*(dx*dx)*A(3,1)
       dfx_i = (A(3,1)*3.0d0*dx + A(2,1)*2.0d0)*dx + A(1,1)
       f10 = f10 + dy*dfx_i

! dfx_2 = A(1,2) + 2*dx*A(2,2) + 3*(dx*dx)*A(3,2)
       dfx_i = (A(3,2)*3.0d0*dx + A(2,2)*2.0d0)*dx + A(1,2)
       f10 = f10 + (dy*dy)*dfx_i

! dfx_3 = A(1,3) + 2*dx*A(2,3) + 3*(dx*dx)*A(3,3)
        dfx_i = (A(3,3)*3.0d0*dx + A(2,3)*2.0d0)*dx + A(1,3)
        f10 = f10 + (dy*dy)*dy*dfx_i

! f10 = dfx_0 + dy*dfx_1 + (dy*dy)*dfx_2 + ((dy*dy)*dy)*dfx_3

! do j=0,ndeg
! f00 = f00 + fx(j)*yv(j)
! f10 = f10 + dfx(j)*yv(j)
! enddo
!
! do j=1,ndeg
! dfy(j) = dble(j)*yv(j-1)
! f01 = f01 + fx(j)*dfy(j)
! enddo

    f00_ = f00
    f10_ = f10
    f01_ = f01

    return
  end subroutine eval_bicub_1



!cdir$r unroll 
   attributes(device) &
  subroutine eval_bicub_2(x,y,xc,yc,A,f00_,f10_,f01_,f11_,f20_,f02_)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(out) :: f00_,f10_,f01_,f11_,f20_,f02_

    real (kind=work_p) :: f00,f10,f01,f11,f20,f02
    real (kind=work_p) :: dx,dy
! real (kind=work_p), dimension(0:ndeg) :: xv,yv
! real (kind=work_p), dimension(0:ndeg) :: fx 
! real (kind=work_p), dimension(0:ndeg) :: dfx
! real (kind=work_p), dimension(0:ndeg) :: dfx2



   real(kind=work_p) :: fx_i,dfx_i,dfx2_i



! f00 = 0.0d0
! f01 = 0.0d0
! f10 = 0.0d0
! f11 = 0.0d0
! f20 = 0.0d0
! f02 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

! xv(0) = 1.0d0
! xv(1) = dx
! xv(2) = dx*dx
! xv(3) = (dx*dx)*dx
! yv(0) = 1.0d0
! yv(1) = dy
! yv(2) = dy*dy
! yv(3) = (dy*dy)*dy


! fx(0:ndeg) = 0.0d0
! dfx(0:ndeg) = 0.0d0
! dfx2(0:ndeg) = 0.0d0


! do j=0,ndeg
! do i=0,ndeg
! fx(j) = fx(j) + xv(i)*A(i,j)
! enddo
! enddo

        f00 = 0
        f01 = 0
        f02 = 0

! fx_0 = A(0,0) + dx*A(1,0) + (dx*dx)*A(2,0) + ((dx*dx)*dx)*A(3,0)
 
        fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0)
        f00 = f00 + fx_i


! fx_1 = A(0,1) + dx*A(1,1) + (dx*dx)*A(2,1) + ((dx*dx)*dx)*A(3,1)

        fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1)
        f00 = f00 + dy*fx_i
        f01 = f01 + fx_i

! fx_2 = A(0,2) + dx*A(1,2) + (dx*dx)*A(2,2) + ((dx*dx)*dx)*A(3,2)

        fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2)
        f00 = f00 + dy*(dy*fx_i)
        f01 = f01 + 2.0d0*(dy*fx_i)
        f02 = f02 + 2.0d0*fx_i
        

! fx_3 = A(0,3) + dx*A(1,3) + (dx*dx)*A(2,3) + ((dx*dx)*dx)*A(3,3)
        fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3)
        f00 = f00 + dy*(dy*(dy*fx_i))
        f01 = f01 + 3.0d0*(dy*(dy*fx_i))
        f02 = f02 + 6.0d0*(dy*fx_i)

!! dfy(j) = dble(j)*yv(j-1)
!! dfy2(j) = dble(j*(j-1))*yv(j-2)
! f01 = f01 + fx(j)*dfy(j), j=1:ndeg
! f02 = f02 + fx(j)*dfy2(j), j=2:ndeg

! f00 = fx_0 + dy*fx_1 + (dy*dy)*fx_2 + ((dy*dy)*dy)*fx_3
! f01 = fx_1 + 2*dy*fx_2 + 3*(dy*dy)*fx_3
! f02 = 2*fx_2 + 3*2*dy*fx_3


! do j=0,ndeg
! do i=1,ndeg
! dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
! enddo
! enddo

        f10 = 0
        f11 = 0


! dfx_0 = A(1,0) + 2*dx*A(2,0) + 3*(dx*dx)*A(3,0)

        dfx_i = (A(3,0)*3.0d0*dx + A(2,0)*2.0d0)*dx + A(1,0)
        f10 = f10 + dfx_i

! dfx_1 = A(1,1) + 2*dx*A(2,1) + 3*(dx*dx)*A(3,1)

        dfx_i = (A(3,1)*3.0d0*dx + A(2,1)*2.0d0)*dx + A(1,1)
        f10 = f10 + dy*dfx_i
        f11 = f11 + dfx_i


! dfx_2 = A(1,2) + 2*dx*A(2,2) + 3*(dx*dx)*A(3,2)

        dfx_i = (A(3,2)*3.0d0*dx + A(2,2)*2.0d0)*dx + A(1,2)
        f10 = f10 + dy*(dy*dfx_i)
        f11 = f11 + 2.0d0*(dy*dfx_i)

! dfx_3 = A(1,3) + 2*dx*A(2,3) + 3*(dx*dx)*A(3,3)

        dfx_i = (A(3,3)*3.0d0*dx + A(2,3)*2.0d0)*dx + A(1,3)
        f10 = f10 + dy*(dy*dy*dfx_i)
        f11 = f11 + 3.0d0*(dy*dy*dfx_i)



! do j=1,ndeg
!! dfy(j) = dble(j)*yv(j-1)
!
! f11 = f11 + dfx(j)*dfy(j)
! enddo


! f11 = f11 + dfx(j)*dfy(j), j=1:ndeg
! f10 = f10 + dfx(j)*yv(j), j=0:ndeg
! f10 = dfx_0 + dy*dfx_1 + (dy*dy)*dfx_2 + ((dy*dy)*dy)*dfx_3
! f11 = dfx_1 + 2*dy*dfx_2 + 3*dy*dy*dfx_3


! do j=0,ndeg
! do i=2,ndeg
! dfx2(j) = dfx2(j) + dble(i*(i-1))*xv(i-2)*A(i,j)
! enddo
! enddo

       f20 = 0

       dfx2_i = (3*2*A(3,0)*dx + 2*1*A(2,0))
       f20 = f20 + dfx2_i

       dfx2_i = (3*2*A(3,1)*dx + 2*1*A(2,1))
       f20 = f20 + dy*dfx2_i

       dfx2_i = (3*2*A(3,2)*dx + 2*1*A(2,2))
       f20 = f20 + (dy*dy)*dfx2_i

       dfx2_i = (3*2*A(3,3)*dx + 2*1*A(2,3))
       f20 = f20 + (dy*dy)*dy*dfx2_i

! do j=0,ndeg
! f20 = f20 + dfx2(j)*yv(j)
! enddo

! f20 = dfx2_0 + dy*dfx2_1 + (dy*dy)*dfx2_2 + ((dy*dy)*dy)*dfx2_3


! do j=1,ndeg
!! dfy(j) = dble(j)*yv(j-1)
!
! f01 = f01 + fx(j)*dfy(j)
! f11 = f11 + dfx(j)*dfy(j)
! enddo

! do j=2,ndeg
!! dfy2(j) = dble(j*(j-1))*yv(j-2)
!
! f02 = f02 + fx(j)*dfy2(j)
! enddo


    f00_ = f00
    f10_ = f10
    f01_ = f01
    f11_ = f11
    f20_ = f20
    f02_ = f02
    

    return
  end subroutine eval_bicub_2
   attributes(device) &
   subroutine I_interpol_wo_pspline(psi, ideriv, ivalue)
   use one_d_cub_mod_gpu
   use precision_mod_gpu
      implicit none

      real (kind=work_p), intent(in) :: psi
      integer, intent(in) :: ideriv
      real (kind=work_p), intent(inout) :: ivalue

      real (kind=work_p) :: pn, wp, acoef(0:ndeg)
      integer :: ip

      pn=psi*one_d_cub_dpsi_inv_gpu
      ip=floor(pn)+1
      ip=min(max(ip,1),ubound(one_d_cub_psi_acoef_gpu,2))
      wp=pn-real(ip-1,kind=work_p)

      acoef=one_d_cub_psi_acoef_gpu(:,ip)

      if (ideriv==0) then
         ivalue=acoef(0)+(acoef(1)+(acoef(2)+acoef(3)*wp)*wp)*wp
      elseif (ideriv==1) then
         ivalue=(acoef(1)+(2.0_work_p*acoef(2)+3.0_work_p*acoef(3)*wp)*wp)*one_d_cub_dpsi_inv_gpu
      else
! print *, 'ideriv in I_interpol_wo_pspline should be 0 or 1'
! stop
      endif
   end subroutine I_interpol_wo_pspline
attributes(device) &
real (kind=8) function I_interpol_gpu(in_psi, ideriv,region)

    use one_d_cub_mod_gpu
    use itp_module_gpu, only : itp_min_psi,itp_max_psi
    use eq_module_gpu, only : eq_x_psi
    use sml_module_gpu, only : sml_bt_sign
    use precision_mod_gpu
    implicit none
    integer , intent(in) :: ideriv,region
    real (kind=work_p) , intent(in) :: in_psi
    real (kind=work_p) :: psi
    integer :: ier
    real (kind=work_p) :: r8value


! sign=-1D0 !-1D0 : cocurrent, 1D0 :counter current 
! sml_bt_sign can be changed in setup.f90 2002/02/08
    
    if(region == 3 ) then
       
       psi=min(eq_x_psi,itp_max_psi) ! for itp_max_psi < eq_x_psi case 2002/01/22
       if(ideriv == 0) then
          call I_interpol_wo_pspline(psi, ideriv, r8value)
          I_interpol_gpu=sml_bt_sign*r8value
       else
          I_interpol_gpu=0D0
       endif
       
    else
       
       psi = in_psi
       if(psi < itp_min_psi) then
          if(psi < itp_min_psi - 1D-4) then
             !print *, 'psi range exceeded',psi
             !call err_count
          endif
          psi=itp_min_psi
       elseif(psi > itp_max_psi) then
          psi=itp_max_psi ! I is constant outside of edge
       endif
       call I_interpol_wo_pspline(psi, ideriv, r8value)
       I_interpol_gpu=sml_bt_sign*r8value
       
    endif
end function I_interpol_gpu
attributes(device) &
subroutine psi_der_all_gpu(r,z,ret)
use bicub_mod_gpu
use precision_mod_gpu
implicit none
  real (kind=work_p) :: r,z,ret(6)
integer, parameter :: i00 = 1
integer, parameter :: i10 = 2
integer, parameter :: i01 = 3
integer, parameter :: i20 = 4
integer, parameter :: i11 = 5
integer, parameter :: i02 = 6

call bicub_interpol2(r,z, &
  ret(i00),ret(i10),ret(i01), &
  ret(i11),ret(i20),ret(i02))

end subroutine psi_der_all_gpu
!!****************************************************************************
!!> calculate field of a given point using interpolation funtions
!! adopted from xorbit
!!
!! first created : 2000/10/19
!! last modified : 2006/02/01
!! B->-B routine added (commented out)
!! time dependance of field is added (2002/6/19)
!! 2006/02/01 fld module update
!! 2002/09/10 code optimization for speed
!! 2002/11/18 code modification for gxc
!!****************************************************************************
attributes(device) &
subroutine field_gpu(fld,t,rz_outside)
    use sml_module_gpu, only :sml_bp_sign, sml_time
    use eq_module_gpu, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z, eq_x_psi, eq_x_z
    use fld_module
    use precision_mod_gpu
    implicit none
    type(fld_type),intent(inout) :: fld !! Field information
    real (kind=work_p), intent(in) :: t !! time
    logical , intent(out) :: rz_outside

    real (kind=work_p) :: psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z
    real (kind=work_p) :: r,z, ret(6)
    !real (kind=8) , external :: I_interpol
    real (kind=work_p) :: r2, over_r,over_r2 !! variables for opimization
    real (kind=work_p) :: cos_rip, sin_rip, ripp, dripp_dr, dripp_dz
    real (kind=work_p) :: rippbphi, drippbphi_dr, drippbphi_dz, drippbphi_dphi

    r=fld%r
    z=fld%z
    r2=r**2
    over_r=1/r
    over_r2=over_r**2
    if(r<eq_min_r) then
        r=eq_min_r
        rz_outside=.true.
    else if (r>eq_max_r)then
        r=eq_max_r
        rz_outside=.true.
    endif
    if(z<eq_min_z) then
        z=eq_min_z
        rz_outside=.true.
    else if (z>eq_max_z)then
       z=eq_max_z
       rz_outside=.true.
    else
       rz_outside=.false.
    endif

    call psi_der_all_gpu(r,z,ret)
    psi =ret(1)
    dpsi_dr =ret(2)
    dpsi_dz =ret(3)
    d2psi_d2r =ret(4)
    d2psi_drdz =ret(5)
    d2psi_d2z =ret(6)

    fld%psi=psi
    fld%dpsidr=dpsi_dr
    fld%dpsidz=dpsi_dz
    ! added 2001/06/01 - lower bound of psi
    if(psi<0D0) then
        psi=0D0
    endif

    !fld_q=q_interpol(psi,0) --> no need
    
! if(psi<eq_x_psi .AND. z<eq_x_z) then
    if(.not. is_rgn12_gpu(r,z,psi) ) then
        fld%I=I_interpol_gpu(psi,0,3)
        fld%dIdpsi = I_interpol_gpu(psi,1,3)
    else
        fld%I=I_interpol_gpu(psi,0,1)
        fld%dIdpsi = I_interpol_gpu(psi,1,1)
    endif

    
    fld%br=- dpsi_dz *over_r * sml_bp_sign
    fld%bz= dpsi_dr *over_r * sml_bp_sign
    fld%bphi=fld%I *over_r 
        
    !derivativs
    fld%dbrdr=dpsi_dz *over_r2 - d2psi_drdz *over_r * sml_bp_sign
    fld%dbrdz=- d2psi_d2z *over_r * sml_bp_sign
    fld%dbrdp=0D0 * sml_bp_sign
    
    fld%dbzdr= - dpsi_dr * over_r2 + d2psi_d2r *over_r * sml_bp_sign
    fld%dbzdz= d2psi_drdz *over_r * sml_bp_sign
    fld%dbzdp=0D0 * sml_bp_sign
           
    fld%dbpdr= dpsi_dr * fld%dIdpsi *over_r - fld%I *over_r2
    fld%dbpdz= fld%dIdpsi * dpsi_dz *over_r
    fld%dbpdp=0D0
! call efield(r,z,fld_psi,t) ! set fld_er, ez, ephi value for a give r,z value -> move to derivs

end subroutine field_gpu
attributes(device) &
  function b_interpol_gpu(r,z,phi)
  use eq_module_gpu, only : eq_x_psi, eq_x_z
  use precision_mod_gpu

  implicit none
  real (kind=work_p) :: b_interpol_gpu
  real (kind=work_p) , intent(in) :: r,z,phi
  real (kind=work_p) :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=work_p) :: ripp,dripp_dr,dripp_dz


  
  call bicub_interpol1(r,z,psi,dpsi_dr,dpsi_dz)
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol_gpu(psi,0,3)
  else
     fi=I_interpol_gpu(psi,0,1)
  endif
  
  br=- dpsi_dz / r ! sign ignored sml_bp_sign
  bz= dpsi_dr / r 
  bphi=fi / r
  
  b_interpol_gpu= sqrt(br**2+bz**2+bphi**2)

end function b_interpol_gpu
attributes(device) &
subroutine bvec_interpol_gpu(r,z,phi,br,bz,bphi)
  use eq_module_gpu
  use sml_module_gpu, only : sml_bp_sign
  use precision_mod_gpu
  implicit none
  real (kind=work_p) , intent(in) :: r,z,phi
  real (kind=work_p) :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=work_p) :: ripp, dripp_dr, dripp_dz

 call bicub_interpol1(r,z,psi,dpsi_dr,dpsi_dz)

 if(psi<eq_x_psi .AND. z<eq_x_z) then
    fi=I_interpol_gpu(psi,0,3)
 else
    fi=I_interpol_gpu(psi,0,1)
 endif

 br=- dpsi_dz / r * sml_bp_sign
 bz= dpsi_dr / r * sml_bp_sign
 bphi=fi / r

end subroutine bvec_interpol_gpu


attributes(device) &
subroutine rho_mu_to_ev_pitch2_gpu(rho,mu,b,ev,pitch,sp_type)
  use sml_module_gpu, only : sml_j2ev
  use ptl_module_gpu, only : ptl_c2_2m, ptl_c_m, ptl_mass
  use precision_mod_gpu
  implicit none
  real (kind=work_p), intent(inout) :: rho,mu,b
  real (kind=work_p), intent(inout) :: ev,pitch
  integer, intent(in) :: sp_type

  real (kind=work_p) :: enj,v_pal,v_pep
  integer, parameter :: idebug = 0

  if(mu<0) then
! if (idebug >= 1) print *, 'minus mu found :',rho,mu,b
    mu=0.0_work_p
  endif
  
  enj=(mu*b+ptl_c2_2m(sp_type)*(rho*b)**2)
  ev=enj*sml_j2ev
  v_pal=ptl_c_m(sp_type)*rho*b
  v_pep=SQRT(2.0_work_p*mu*b/ptl_mass(sp_type))
  pitch=v_pal/SQRT(v_pal**2+v_pep**2)
  
end subroutine rho_mu_to_ev_pitch2_gpu
attributes(device) &
subroutine remove_particle_gpu( i,flag)
  use ptl_module_gpu, only : ptl_ph_gpu,ptl_ct_gpu,ptl_gid_gpu, &
                             ptl_nphase,type_gpu,ptl_nconst,pim,piw0
  use sml_module_gpu, only : sml_mype
  use neu_module_gpu, only : neu_weight_sum_lost_gpu
  use precision_mod_gpu
  implicit none

  integer, intent(in) :: i,flag

! integer :: sp%type
! integer(8) :: sp%ptl(i)%gid
! real(8), dimension(ptl_nphase) :: sp%ptl(i)%ph
! real(8), dimension(ptl_nconst) :: sp%ptl(i)%ct
  !error diagnosis
  real (kind=work_p) :: b_interpol_gpu, b,ekin,pitch
  integer :: ip
  real (kind=work_p) :: dummy
  integer, parameter :: idebug = 0


! sp%ptl(i)%gid = sp%ptl(i)%gid
! sp%ptl(i)%ph = sp%ptl(i)%ph
! sp%ptl(i)%ct = sp%ptl(i)%ct
! sp%type = sp%type
  
  if(ptl_gid_gpu(i) <= 0 ) then
! if (idebug >= 1) print *, 'minus gid particle in remove_particle_gpu'
    return
  endif
  ptl_gid_gpu(i)=-ptl_gid_gpu(i)

  !### debug
   if(ptl_gid_gpu(i) > 0 ) then
     if (idebug >= 1) print *, 'something wrong in remove_particle_gpu'
  endif

!if (idebug >= 1) print *, 'sml_neutral', sml_neutral

  if(flag==-2 .and. sml_mype==0) then
     b=b_interpol_gpu(ptl_ph_gpu(i,1),ptl_ph_gpu(i,2),0.0_work_p)
     call rho_mu_to_ev_pitch2_gpu(ptl_ph_gpu(i,4), &
                 ptl_ct_gpu(i,pim), b,ekin,pitch,type_gpu)
! write(400+sp%type,*) sp%ptl(i)%ph(1:2),&
! ekin, pitch,& 
! sp%phase0(1,i),sp%phase0(2,i) 
  endif

  if(sml_neutral .and. flag==-1.and.(ptl_gid_gpu(i)<0).and.(type_gpu==1))then
! neu_weight_sum_lost_gpu(ith)=neu_weight_sum_lost_gpu(ith) + ptl_ct_gpu(i,piw0)
     ip = mod( i, size(neu_weight_sum_lost_gpu)) + 1
     dummy = atomicAdd( neu_weight_sum_lost_gpu(ip), ptl_ct_gpu(i,piw0))
  endif

  if(flag==-1 .and. sml_mype==0) then
! write(450+sp%type,*) sp%phase(1,i), sp%phase(2,i),&
! ekin, pitch,& 
! sp%phase0(1,i),sp%phase0(2,i) 
  endif

! if(flag==-3) then
! print *, 'unexpected search fail in search_pid', &
! sp%phase(1,i), sp%phase(2,i)
! endif
end subroutine remove_particle_gpu

attributes(device) &
real (kind=8) function z_psi_min_gpu(r)
  use sml_module_gpu, only : sml_mype
  use bnc_module_gpu, only : bnc_min_r, bnc_nr, bnc_dr, bnc_z_psi_min
  use precision_mod_gpu
  implicit none
  real (kind=work_p),intent(in) :: r
  integer :: i
  real (kind=work_p) :: aa,bb
  
  i=int((r-bnc_min_r)/bnc_dr) +1
  i=min(bnc_nr-1,max(1,i))
  bb=(r-bnc_min_r)/bnc_dr + 1D0 -i
  aa=1D0-bb
! print *, 'zaa', r, i, aa,bb,bnc_z_psi_min(i)
  z_psi_min_gpu=bnc_z_psi_min(i)*aa + bnc_z_psi_min(i+1)*bb

end function z_psi_min_gpu
attributes(device) &
real (kind=8) function b_interpol_sym_gpu(r,z) ! axi-symmetric b interpolation
  use eq_module_gpu
  use bicub_mod_gpu
  use precision_mod_gpu

  implicit none
  real (kind=work_p) , intent(in) :: r,z
  real (kind=work_p) :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=work_p) :: ripp,dripp_dr,dripp_dz

 call bicub_interpol1(r,z,psi,dpsi_dr,dpsi_dz)

  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol_gpu(psi,0,3)
  else
     fi=I_interpol_gpu(psi,0,1)
  endif

  br=- dpsi_dz / r
  bz= dpsi_dr / r
  bphi=fi / r

  b_interpol_sym_gpu= sqrt(br**2+bz**2+bphi**2)

end function b_interpol_sym_gpu


  ! function evaluation
  attributes(device) &
    function eq_ftn_gpu2(p,r,z,ftn)
    use eq_module_gpu, only : eq_x_psi, eq_x_z,eq_ftn_type
    use precision_mod_gpu
    !use sml_module
    implicit none
    real(kind=work_p) :: eq_ftn_gpu2
    type(eq_ftn_type) :: ftn
    real(kind=work_p) :: p,r, z 


    real(kind=work_p) :: tmp, tmp2, tmp3, tmp4

   if(is_rgn12_gpu(r,z,p)) then
       tmp=p
    else
       tmp=eq_x_psi
    endif

    !tmp=min(max(p,sml_pmin), sml_pmax)

    select case(ftn%shape)
    case(0) ! constant
       eq_ftn_gpu2 = ftn%iny(1)
    case(1) ! hyperbolic tanh
       eq_ftn_gpu2 = ftn%sv(2)*tanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
    case(2) ! linear
       eq_ftn_gpu2 = ftn%sv(1)*tmp + ftn%sv(2)
    case(3) ! a exp(-B w tanh( (r-r0)/w ))
       eq_ftn_gpu2 = ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2)))
    case(4) 
       eq_ftn_gpu2 = ftn%sv(2)*tanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
       if(tmp < ftn%inx(1)-0.5_work_p*ftn%inx(2) ) then
          eq_ftn_gpu2 = eq_ftn_gpu2 + ftn%sv(4)*sqrt(tmp) + ftn%sv(5)
       endif
    case(5)
       tmp2=ftn%sv(3)*( ftn%inx(1)-sqrt(ftn%inx(1)*tmp) ) ! z
       tmp3=exp(tmp2) ! expz
       tmp4=exp(-tmp2) ! expmz
       ! A * ( (1+z*slope)*expz - expmz )/(expz + expmz) + B
       eq_ftn_gpu2 = ftn%sv(2)*( (1+tmp2*ftn%sv(4))*tmp3 - tmp4 )/(tmp3+tmp4) + ftn%sv(1)
    case(-1)
       ! print*,'shape number -1 not implemented in eq_ftn_gpu2'
       ! tmp=min(max(tmp,ftn_min),ftn_max)
       ! call ftn_evaluation(ftn_spl,tmp,eq_ftn_gpu2)
       tmp=min(max(tmp,ftn%min),ftn%max)
       eq_ftn_gpu2 = huge(0.0_work_p)
! print *, "GPU code is not properly implemented for this case"
    case default
       ! print *, 'Invalid shape number in eq_ftn_gpu2', ftn%shape
       ! stop
    end select
    
  end function eq_ftn_gpu2
  ! function evaluation
  attributes(device) &
    function eq_dftn_gpu2(p,r,z, ftn)
         
    use eq_module_gpu, only : eq_x_z, eq_x_psi,eq_ftn_type
    use precision_mod_gpu
    ! use sml_module
    implicit none
    real(kind=work_p) :: eq_dftn_gpu2
    type(eq_ftn_type) :: ftn
    real(kind=work_p) :: p,r, z 
    real(kind=work_p) :: tmp


    ! region searching and enforcing region 3 value to x-point value
    if(is_rgn12_gpu(r,z,p)) then
       tmp=p
    else
       tmp=eq_x_psi
       eq_dftn_gpu2=0.0_work_p
       return
    endif


    select case(ftn%shape)
    case(0) ! constant
       eq_dftn_gpu2=0
    case(1) ! hyperbolic tanh
       eq_dftn_gpu2 = -ftn%sv(3)*ftn%sv(2)*(1.0_work_p/cosh((ftn%inx(1)-tmp)*ftn%sv(3)))**2
    case(2) ! linear
       eq_dftn_gpu2 = ftn%sv(1)
    case(3)
       eq_dftn_gpu2 = - ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2))) &
            /((cosh((ftn%inx(1)-tmp)*ftn%sv(2)))**2*ftn%inx(3))
       
       
    case(-1)
       if(ftn%min < tmp .and. tmp < ftn%max) then
          eq_dftn_gpu2=-huge(0.0_work_p)
! print *, "GPU code is not properly implemented for this case"
       else
          ! EFD not sure how to work around this
          ! call ftn_derivative(ftn_spl,tmp,1,eq_dftn_gpu2)
          eq_dftn_gpu2 = -0.0_work_p
          ! print*,'eqn_dftn_gpu2: ftn_derivative not implemented'
       endif
    case default
       ! print *, 'Invalid shape number in eq_dftn_gpu2', ftn%shape
       ! stop
    end select
    
    ! case 4 and case 5 does not have derivative function, yet. --> case 5 is required

  end function eq_dftn_gpu2
  ! function evaluation
  attributes(device) &
  logical function is_rgn1_gpu(r,z,psi)
    use precision_mod_gpu
    use eq_module_gpu, only : eq_x_psi, eq_x_r, eq_x_z, eq_x_slope, &
                              eq_x2_r, eq_x2_z, eq_x2_slope
    implicit none
    real (kind=work_p), intent(in) :: r,z,psi
    real (kind=work_p) , parameter :: epsil_psi = 1D-5

    if(psi <= eq_x_psi-epsil_psi .and. &
         -(r-eq_x_r)*eq_x_slope + (z-eq_x_z) > 0.0_work_p .and. -(r-eq_x2_r)*eq_x2_slope &
         + (z-eq_x2_z) < 0.0_work_p ) then
       is_rgn1_gpu=.true.
    else
       is_rgn1_gpu=.false.
    endif
  end function is_rgn1_gpu


  attributes(device) &
  logical function is_rgn12_gpu(r,z,psi)
    use eq_module_gpu, only : eq_x_r, eq_x_z, eq_x_psi, eq_x2_r, eq_x2_z, eq_x2_slope
    use precision_mod_gpu
    implicit none
    real (kind=work_p) :: r,z,psi
    real (kind=work_p) , parameter :: epsil_psi = 1D-5

    if(psi > eq_x_psi -epsil_psi .or. &
         -(r-eq_x_r)*eq_x_slope + (z-eq_x_z) > 0.0_work_p .and. -(r-eq_x2_r)*eq_x2_slope + (z-eq_x2_z) < 0.0_work_p ) then
       is_rgn12_gpu=.true.
    else
       is_rgn12_gpu=.false.
    endif
  end function is_rgn12_gpu

attributes(device) &
subroutine efield_gk_elec_gpu(i,fld,itr,p) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module_gpu, only : sml_turb_efield, sml_mype, sml_turb_efield
  use grid_class_gpu, only : grid_delta_phi, grid_nd, grid_rhomax, grid_drho, grid_nrho 

  use psn_class_gpu, only : psn_pot_rho_ff, &
                            psn_E_rho_ff, &
                            psn_E_phi_ff, &
                            psn_ddpotdt_phi !, psn_ddpotdt, psn_E00_ff,psn_pot_phi_ff


  use fld_module, only : fld_type
! use ptl_module_gpu, only : type_gpu, rhoi_gpu
  use precision_mod_gpu
  implicit none
  integer, intent(in) :: i,itr
  type(fld_type), intent(inout) :: fld

  integer :: ip,node,irho, iphi
  real (kind=work_p) :: wphi(0:1), wp, rho, rhon, wrho(2), p(3)
  real (kind=work_p) :: E(3), B !, D(3), pot, E00(2), ddpotdt


  iphi=floor(fld%phi/grid_delta_phi)
! wphi(1)=(fld%phi/grid%delta_phi) - grid%iphi_offset
  wphi(1)=(fld%phi/grid_delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)

! pot=0D0
  E=0D0
! E00=0D0
  !D=0D0 !debug
! ddpotdt=0D0

  !get E
 
  if(itr>0) then 
     do ip=1, 3
        node=grid_nd(ip,itr)
        wp=p(ip)
        
           ! for electron -- rho=0 case, optimized. 
! pot = pot + wp*wphi(0)*psn_pot_phi_ff(0,node,iphi) 
! pot = pot + wp*wphi(1)*psn_pot_phi_ff(1,node,iphi)



           E = E + wp*wphi(0)*psn_E_phi_ff(:,0,node,iphi)
           E = E + wp*wphi(1)*psn_E_phi_ff(:,1,node,iphi)


! E00 = E00 + wp*wphi(0)*psn_E00_ff(:,0,node)
! E00 = E00 + wp*wphi(1)*psn_E00_ff(:,1,node)

! ddpotdt = ddpotdt + wp*wphi(0)*psn_ddpotdt_phi(node,0,iphi)
! ddpotdt = ddpotdt + wp*wphi(1)*psn_ddpotdt_phi(node,1,iphi)



     enddo
  else
! print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
  endif
    
  !E(3) was para E-field and becomes Ephi
  if(sml_turb_efield) then
     B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     E(3)=(E(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi

     !debug
     !D(3)=(D(3)*B- D(1)*fld%Br - D(2)*fld%Bz)/fld%Bphi
  else
     E(3)=0D0
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
! fld%Epot=pot
! fld%Er00=E00(1)
! fld%Ez00=E00(2) 
! fld%ddpotdt=ddpotdt


end subroutine efield_gk_elec_gpu
attributes(device) &
subroutine derivs_sp_elec_gpu(fld,t, ptli, &
               yprime, diag_on, vf_diag)
  use fld_module, only : fld_type
  use sml_module_gpu, only : &
    sml_n_vf_diag, sml_inpsi, sml_outpsi, sml_deltaf, &
    sml_e_charge, sml_deltaf_f0_mode, sml_f0_1_lt_e, &
    sml_f0_1_lt, sml_f0_1_ln, sml_dwdt_fix_bg,sml_dwdt_exb_only, sml_ignore_drift_near_wall
  use ptl_module_gpu, only : pim, pirho, piw1,piw2,&
      ptl_nphase, ptl_charge, ptl_mass, &
      ptl_nconst, ptl_ph_gpu, ptl_ct_gpu, ptl_deltaf_sp
  use eq_module_gpu, only : eq_den, eq_tempe, eq_tempi
  use precision_mod_gpu
  implicit none
  type(fld_type), intent(in) :: fld ! field variable
  real (kind=work_p), intent(in) :: t ! time

  integer, intent(in) :: ptli ! particle info
! integer, intent(in) :: sp_type ! particle species type (ion/elec)
  real (kind=work_p), intent(inout) :: yprime(ptl_nphase) !
  logical, intent(in) :: diag_on
  real (kind=work_p), intent(inout) :: vf_diag(sml_n_vf_diag) ! variables for diagnosis 
  !
  real (kind=work_p) :: mass, charge,c_m ! charge and mass 
  real (kind=work_p) :: r, z, phi, rho, mu,inv_r
  real (kind=work_p) :: B, B2, over_B, over_B2
  real (kind=work_p) :: D, nb_curl_nb
  real (kind=work_p) :: dbdr, dbdz, dbdphi
  real (kind=work_p) :: cmrho2, cmrho, murho2b, murho2b_c, vp
  real (kind=work_p) :: fr, fp, fz
  real (kind=work_p) :: fr_exb, fp_exb, fz_exb, yp_exb(3)

  ! for weight calculation
  real (kind=work_p) :: energy, vmag, pitch, dvpdt, denergy_dt, dpitch_dt
  real (kind=work_p) :: psi, den, dden, temp, dtemp, dfdp, tmp, envelop, df0(5)
  real (kind=work_p) :: one_m_w
  real (kind=work_p) :: total_ddpotdt ! for electron -- from adiabatic reponse
  !
  real (kind=work_p) :: bp, f0_1_lt
  real (kind=work_p) :: i_factor
  ! some global parameter to be
  real (kind=work_p) :: sml_wdot_energy_max, sml_f0_psi_c, sml_f0_1_psi_w
  
  ! these are constant for deltaf method -- maybe should be moved to sml_module

  if(sml_deltaf_f0_mode==-1) then
     ! these are constant for deltaf method -- maybe should be moved to sml_module
     sml_wdot_energy_max=10D0
     sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
     sml_f0_1_psi_w=1D0/( 0.4*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )
  endif

  ! prepare constant
  mass=ptl_mass(0) !
  charge=ptl_charge(0) !-1.D0/ptl_charge
  c_m=charge/mass
    
  r= ptl_ph_gpu(ptli,1)
  z= ptl_ph_gpu(ptli,2)
  phi= ptl_ph_gpu(ptli,3)
  rho= ptl_ph_gpu(ptli,4)
  mu= ptl_ct_gpu(ptli,pim) 
  inv_r=1D0/r

  B = sqrt( fld%br**2 + fld%bz**2 + fld%bphi **2 )
  b2=b*b
  over_B=1/B
  over_B2=over_B*over_B
  
  ! normalized b dot curl of normalized b
  nb_curl_nb= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/r + fld%dbpdr)*fld%bz )
  D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
      
  dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  dbdphi=0D0 ! no B perturbation


  ! bug fix -- rho**2 B grad B term added 2002/1/15
  !! optimization variables
  
  cmrho2=c_m*rho**2
  cmrho =c_m*rho
  murho2b=(mu+charge*cmrho2 *B)
  murho2b_c=murho2b/charge
  vp=cmrho*B

  ! F is grad H/q
  fr = fld%Er - (murho2B_c) * dbdr ! fr is force over q . -gradient of Hamiltonian
  fp = fld%Ephi - (murho2B_c) * dbdphi/r ! modified by shlee 5/30/2001 -- /r added
  fz = fld%Ez - (murho2B_c) * dbdz

  ! ignore all drifts when it close to divertor
  ! ion & electron --> pushe for electron
  ! now drift and parallel E-field are ignored.
  if(sml_ignore_drift_near_wall ) then
     call ignore_factor_gpu(r,z,i_factor)
     !ep_r = (fld%Er *fld%br )*fld%br /b2
     !ep_p = (fld%Ephi*fld%bphi)*fld%bphi/b2 
     !ep_z = (fld%Ez *fld%bz )*fld%bz /b2 

     fr = fr*i_factor !+ ep_r*(1D0-i_factor)
     fp = fp*i_factor !+ ep_p*(1D0-i_factor)
     fz = fz*i_factor !+ ep_z*(1D0-i_factor) 
  endif
  
  yprime(1)= D*( (fld%bz*Fp - fld%Bphi * Fz) * over_B2 &
       + cmrho * fld%br &
       + cmrho2 * (fld%dbzdp*inv_r - fld%dbpdz ) )
  yprime(2)= D*( (fld%bphi * fr - fld%br * fp ) * over_B2 &
       + cmrho * fld%bz &
       + cmrho2 * (fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r) )
  yprime(3)= D*( (fld%br * fz - fld%bz * fr) * over_B2 &
       + cmrho * fld%bphi &
       + cmrho2 * ( fld%dbrdz - fld%dbzdr) ) * inv_r


  fr_exb = fld%Er ! fr is force over q . -gradient of Hamiltonian
  fp_exb = fld%Ephi! modified by shlee 5/30/2001 -- /r added
  fz_exb = fld%Ez 
  
  yp_exb(1)= D*(fld%bz *fp_exb - fld%Bphi * fz_exb) * over_B2 
  yp_exb(2)= D*(fld%bphi *fr_exb - fld%br * fp_exb) * over_B2 
  yp_exb(3)= D*(fld%br *fz_exb - fld%bz * fr_exb) * over_B2 *inv_r

  yprime(pirho)=D*over_B2 *( &
       fld%br*fr + fld%bz*fz + fld%bphi*fp &
       + rho*( fr*(fld%dbzdp*inv_r-fld%dbpdz) + fz*(fld%bphi*inv_r+fld%dbpdr-fld%dbrdp*inv_r) + fp*(fld%dbrdz-fld%dbzdr)) &
       )
  
  if( ptl_deltaf_sp(0) .and. sml_extra_dwdt ) then
     energy=(0.5D0*mass*vp*vp + mu*b)
     vmag=sqrt(2D0*energy/mass)
     pitch = vp/vmag
     dvpdt=c_m*yprime(pirho)*B + (dbdr*yprime(1)+dbdz*yprime(2))*vp/B !
     denergy_dt=mu*(dbdr*yprime(1)+dbdz*yprime(2)+dbdphi*yprime(3))+(mass*vp*dvpdt)
     dpitch_dt=dvpdt/vmag - 0.5D0*pitch*denergy_dt/energy


          ! dlog(f0) ------ need to be separate routine later 
     ! temp routine -- finish it later
     ! no space gradient
     psi=fld%psi

     !den=eq_ftn(psi,z,eq_den)
     den = eq_ftn_gpu2(psi,r,z,eq_den)


     !dden=eq_dftn(psi,z,eq_den)
     dden = eq_dftn_gpu2(psi,r,z,eq_den)

! if(sp_type==0) then !electron
        !temp=eq_ftn(psi,z,eq_tempe)*sml_e_charge
        temp=eq_ftn_gpu2(psi,r,z,eq_tempe)*sml_e_charge
        ! dtemp=eq_dftn(psi,z,eq_tempe)*sml_e_charge
        dtemp=eq_dftn_gpu2(psi,r,z,eq_tempe)*sml_e_charge
! else
        ! temp=eq_ftn(psi,z, eq_tempi)*sml_e_charge
! temp=eq_ftn_gpu2(psi,r,z, eq_tempi)*sml_e_charge
        ! dtemp=eq_dftn(psi,z,eq_tempi)*sml_e_charge
! dtemp=eq_dftn_gpu2(psi,r,z,eq_tempi)*sml_e_charge
! endif
     
     if(sml_deltaf_f0_mode==-2) then !consistent grad f with real profile
        dfdp=dden/den + dtemp/temp*(energy/temp - 1.5D0)
     elseif(sml_deltaf_f0_mode==-1) then
        tmp=1D0/sqrt(fld%dpsidr**2+fld%dpsidz**2) ! drdpsi 
        envelop= exp( - ((sqrt(psi) - sml_f0_psi_c)*sml_f0_1_psi_w )**8 )
        !envelop=1D0
! if(sp_type==0) then
           f0_1_Lt=sml_f0_1_Lt_e
! else
! f0_1_Lt=sml_f0_1_Lt
! endif
        dfdp=-envelop*(sml_f0_1_Ln + f0_1_Lt*(energy/temp - 1.5D0))*tmp
     endif
     df0(1)=dfdp *fld%dpsidr !/(2D0*q+1D-99)
     df0(2)=dfdp *fld%dpsidz !/(2D0*q+1D-99)
     df0(3)=0D0 * r ! r : yprime(3) is phi dot --> dimension matching
     df0(4)=0D0 ! zero for extra dwdt
     df0(5)=0D0 ! pitch deriv

     if(sml_dwdt_fix_bg)then 
        one_m_w=1D0
     else
        one_m_w=1D0 - ptl_ph_gpu(ptli,piw2)
     end if

     if( .not. sml_dwdt_exb_only ) then
        yprime(piw1)= -one_m_w* (&
             df0(1)*yprime(1) + df0(2)*yprime(2) + df0(3)*yprime(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )
     else if(sml_deltaf_f0_mode == -1) then
        yprime(piw1)= -one_m_w* (&
             df0(1)*yp_exb(1) + df0(2)*yp_exb(2) + df0(3)*yp_exb(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )
     else if(sml_deltaf_f0_mode == -2) then ! remove weight evolution from v_grad B - total will be updated later
        yprime(piw1)= -one_m_w* (&
             df0(1)*(yp_exb(1)-yprime(1)) + df0(2)*(yp_exb(2)-yprime(2)) )

     endif
     ! electron -- minus adibatic response
! if(sp_type==0) then
        total_ddpotdt=fld%ddpotdt - yprime(1)*(fld%Er-fld%Er00) - yprime(2)*(fld%Ez-fld%Ez00) -r*yprime(3)*fld%Ephi

        yprime(piw1) = yprime(piw1) - one_m_w*total_ddpotdt/temp*sml_e_charge
! endif
     yprime(piw2)=yprime(piw1)
  else
     yprime(piw1:piw2)=0D0
  endif
  if(diag_on) then
     vf_diag(1)=fld%psi
     vf_diag(2)=B
     vf_diag(3)=fld%Bphi
     vf_diag(4)=yp_exb(1)*fld%dpsidr + yp_exb(2)*fld%dpsidz ! V_ExB dot grad psi
     vf_diag(5)=yprime(1)*fld%dpsidr + yprime(2)*fld%dpsidz ! V_d dot grad psi
     vf_diag(6)=fld%dpsidr*fld%dpsidr + fld%dpsidz*fld%dpsidz ! |grad psi|^2
     bp=sqrt( fld%br*fld%br + fld%bz*fld%bz ) 
     vf_diag(7)=(yprime(1)*fld%br + yprime(2)*fld%bz)/bp
     vf_diag(8)=(yp_exb(1)*fld%br + yp_exb(2)*fld%bz)/bp
  endif

end subroutine derivs_sp_elec_gpu
attributes(device) &
subroutine efield_gpu(i,fld,time,itr,p)
  use sml_module_gpu, only : &
     sml_00_efield, sml_turb_efield
  use fld_module
  !use eq_module
  use precision_mod_gpu
  implicit none
  type(fld_type), intent(inout) :: fld
  integer, intent(in) :: i,itr
  real (kind=work_p), intent(in) :: time, p(3)

  
  if(sml_00_efield .or. sml_turb_efield) then
     !Gyrokinetic E
     call efield_gk_elec_gpu(i,fld,itr,p)
  else
       fld%Er=0D0
       fld%Ez=0D0
       fld%Ephi=0D0
       fld%Epot=0D0
       fld%Er00=0D0
       fld%Ez00=0D0
       fld%ddpotdt=0D0
  endif

end subroutine efield_gpu
! port1 called in push 
attributes(device) &
subroutine diag_1d_port1_gpu(ptli,derivs,sp_type,vd,ith)



  use ptl_module_gpu, only : ptl_ph_gpu,ptl_ct_gpu, ptl_gid_gpu
  use diag_module, only : diag_1d_npv1
  use sml_module_gpu, only : sml_n_vf_diag, sml_deltaf, sml_ev2j, sml_mype, sml_inpsi
  use eq_module_gpu, only : eq_x_z, eq_x2_z, eq_axis_r, eq_axis_z, eq_tempi
  use diag_module_gpu !, only : diag_1d_npv1, diag_1d_f_pv1, diag_1d_df_pv1, diag_1d_pin, &
! diag_1d_dp_inv, diag_1d_npsi
  use precision_mod_gpu
  implicit none
  integer, intent(in) :: ptli
  real (kind=work_p), intent(in) :: derivs(ptl_nphase)
  integer, intent(in) :: sp_type
  real (kind=work_p) :: vd(sml_n_vf_diag) ! 1: R_major, 2: B_toroidal, 3: B_total 4: radial ExB 5: v_para

  integer, intent(in) :: ith
  !
  real (kind=work_p) :: psi,z, pn, wp, b, r, rho, mu, w, dw, vp, en, den, pe, we
  real (kind=work_p) :: diag_1d_de
  integer :: ip, j, ie
  real (kind=work_p) :: v(diag_1d_npv1)

  real (kind=work_p) :: dval0, dummy0,dval1,dummy1
  real (kind=work_p) :: dval2,dummy2,dval3,dummy3,dval4,dummy4,dval5,dummy5,dval6,dummy6,dval7,dummy7
  integer :: lb,ub,i,i0,ith0


  if(ptl_gid_gpu(ptli) > 0) then
     psi = vd(1)
     r=ptl_ph_gpu(ptli,1)
     z=ptl_ph_gpu(ptli,2)
     if(.not. is_rgn12_gpu(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) return ! EXIT for private region of lower divertor
     pn=(psi-diag_1d_pin)*diag_1d_dp_inv
     ip=floor(pn)+1
     if(ip <1 .or. diag_1d_npsi <= ip) return ! EXIT for out of diag_1d_pin/pout range
     wp=1.0_work_p - pn + real(ip-1,kind=work_p) 
     
     ! local variables for readability 
     b=vd(2)
! r=ptl_ph_gpu(ptli,1)
     rho=ptl_ph_gpu(ptli,4)
     mu=ptl_ct_gpu(ptli,1)
     w=ptl_ct_gpu(ptli,2)
     dw=ptl_ph_gpu(ptli,5)*w
     vp=ptl_c_m(sp_type)*rho*B

     ! obtain variables
     v(1) = 1D0 ! full weight
     v(2) = derivs(3)*r ! g.c. toroidal velocity -- dzeta/dt * R_major
     v(3) = vd(7) ! g.c. poloidal velocity -- v dot Bp/Bp
     v(4) = vp ! v_|| parallel velocity 
     v(5) = vp*r*vd(3)/vd(2) ! v_zeta_|| / R --toroidal angular momentum of v_|| : v_|| * B_zeta /(R*B)
     v(6) = vd(5) ! radial drift - psi_dot
     v(7) = ptl_c2_2m(sp_type)*(rho*B)**2 ! parallel mean energy
     v(8) = mu*B ! perp temperature
     v(9) = (v(7)+v(8))*v(6) ! Energy radial flux
     v(10)= v(4)*vd(4)*vd(3)/vd(2) ! V_exb * V_|| * B_phi/B
     v(11)= v(4)*vd(5)*vd(3)/vd(2) ! V_r * V_|| * B_phi/B
     v(12)= vd(4) ! radial drift by exb - V_exb dot grad_psi
     v(13)= (v(7)+v(8))*vd(4) ! heat flux by radial exb
     v(14)= vd(8) ! poloidal comp. of V_ExB
     v(15)= vd(6) ! grad_psi ^2


     ith0 = 1 + mod(ith, size(diag_1d_f_pv1,4))
     lb = lbound(diag_1d_f_pv1,1)
     ub = ubound(diag_1d_f_pv1,1)

     do i=lb,ub
        i0 = (i-lb) + lbound(v,1)
        dval0 = (w * wp)*v(i0)
        dummy0 = atomicAdd( diag_1d_f_pv1(i,ip,sp_type,ith0), dval0 )

        dval1 = (w * (1.0_work_p-wp)) * v(i0)
        dummy1 = atomicAdd( diag_1d_f_pv1(i,ip+1,sp_type,ith0), dval1)
      enddo


     if(ptl_deltaf_sp(sp_type)) then

         ith0 = 1 + mod(ith,size(diag_1d_df_pv1,4))
         lb = lbound(diag_1d_df_pv1,1)
         ub = ubound(diag_1d_df_pv1,1)

         do i=lb,ub
           i0 = (i-lb) + lbound(v,1)
           dval0 = (dw*wp)*v(i0)
           dummy0 = atomicAdd( diag_1d_df_pv1(i,ip ,sp_type,ith0), dval0)

           dval1 = dw*(1.0_work_p-wp) * v(i0)
           dummy1 = atomicAdd(diag_1d_df_pv1(i,ip+1,sp_type,ith0), dval1)
         enddo



     endif

! if(diag_tavg_on .and. diag_omid_on) then
! if( r-eq_axis_r > abs(z-eq_axis_z) ) then
! diag_1d_omid_f_pv1(:,ip ,sp_type,ith)=v(:)*w*wp +diag_1d_omid_f_pv1(:,ip ,sp_type,ith)
! diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1.0_work_p-wp) +diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)
! endif
! endif
! if(diag_eflux_on)then
     if(.true.)then
        diag_1d_emin = 0.0_work_p
        diag_1d_emax = 3.0_work_p*eq_ftn_gpu2(sml_inpsi,eq_axis_r,eq_axis_z,eq_tempi)*sml_ev2j
        diag_1d_de=(diag_1d_emax-diag_1d_emin)/diag_1d_ne

        en=v(7)+v(8)
        pe=(en-diag_1d_emin)/diag_1d_de
        ie=floor(pe)+1
        if(ie >=1 .and. diag_1d_ne > ie) then ! EXIT for out of diag_1d_emin/emax range
          we=1.0_work_p - pe + real(ie-1,kind=work_p)
          ith0 = 1 + mod(ith, size(diag_1d_eflux_pv,5))

          dval0 = w * wp * we
          dummy0 = atomicAdd( diag_1d_eflux_pv(1,ip,ie,sp_type,ith0), dval0)

          dval1 = w * (1.0_work_p-wp) * we
          dummy1 = atomicAdd( diag_1d_eflux_pv(1,ip+1,ie,sp_type,ith0), dval1)

          dval2 = w * wp * (1.0_work_p-we)
          dummy2 = atomicAdd( diag_1d_eflux_pv(1,ip,ie+1,sp_type,ith0), dval2 )

          dval3 = w * (1.0_work_p-wp) * (1.0_work_p-we)
          dummy3 = atomicAdd( diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith0), dval3)

          dval4 = w * wp * we * v(9)
          dummy4 = atomicAdd( diag_1d_eflux_pv(2,ip,ie,sp_type,ith0), dval4 )

          dval5 = w * (1.0_work_p-wp) * we * v(9)
          dummy5 = atomicAdd( diag_1d_eflux_pv(2,ip+1,ie,sp_type,ith0), dval5)

          dval6 = w * wp * (1.0_work_p-we) * v(9)
          dummy6 = atomicAdd( diag_1d_eflux_pv(2,ip,ie+1,sp_type,ith0), dval6 )

          dval7 = w * (1.0_work_p-wp) * (1.0_work_p-we) * v(9)
          dummy7 = atomicAdd( diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith0), dval7)

          if(ptl_deltaf_sp(sp_type)) then
            ith0 = 1 + mod(ith, size(diag_2d_dflux_pv,5))
         
            dval0 = dw * wp * we
            dummy0 = atomicAdd( diag_2d_dflux_pv(1,ip,ie,sp_type,ith0), dval0 )

            dval1 = dw * (1.0_work_p-wp) * we
            dummy1 = atomicAdd( diag_2d_dflux_pv(1,ip+1,ie,sp_type,ith0), dval1)

            dval2 = dw * wp * (1.0_work_p-we)
            dummy2 = atomicAdd( diag_2d_dflux_pv(1,ip,ie+1,sp_type,ith0), dval2 )

            dval3 = dw * (1.0_work_p-wp) * (1.0_work_p-we) 
            dummy3 = atomicAdd( diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith0), dval3)

            dval4 = dw * wp * we * v(9)
            dummy4 = atomicAdd( diag_2d_dflux_pv(2,ip,ie,sp_type,ith0), dval4 )

            dval5 = dw * (1.0_work_p-wp) * we * v(9)
            dummy5 = atomicAdd( diag_2d_dflux_pv(2,ip+1,ie,sp_type,ith0), dval5)

            dval6 = dw * wp * (1.0_work_p-we) * v(9)
            dummy6 = atomicAdd( diag_2d_dflux_pv(2,ip,ie+1,sp_type,ith0), dval6 )

            dval7 = dw * (1.0_work_p-wp) * (1.0_work_p-we) * v(9)
            dummy7 = atomicAdd( diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith0), dval7)


          endif ! ptl_deltaf_sp
        endif ! ie
     endif ! diag_eflux_on
  endif
  
end subroutine diag_1d_port1_gpu
attributes(device) &
subroutine diag_heat_port_gpu(w, pot, epara, eperp, ct, old_ph, new_ph, dphi, stype, ith)
    use sml_module_gpu
    use ptl_module_gpu
    use diag_module_gpu !, only : diag_heat_nvar
    use precision_mod_gpu

    implicit none
    real (kind=work_p), intent(in) :: w, pot, epara, eperp
    real (kind=work_p), intent(in) :: ct(ptl_nconst), old_ph(ptl_nphase), new_ph(ptl_nphase), dphi
    integer, intent(in) :: stype, ith
    real (kind=work_p) :: wp, r, z, psi, rn, zn, pn, wr,wz, ws, v(5)
    integer :: ir, iz, ip, itmp

    real (kind=work_p) :: dval0, dummy0,dval1,dummy1, dval2, dummy2, dval3, dummy3
    integer :: lb,ub,i,i0,ith0
    real (kind=work_p) :: x(2), phi, xff(2), phi_mid

    ! chracteristic r and z - use mean value for simplicity
! r=(old_ph(pir)+new_ph(pir))*0.5_work_p
! z=(old_ph(piz)+new_ph(piz))*0.5_work_p

    !get field following
    x=new_ph(1:2)
    phi=new_ph(3)
    phi_mid=(floor(phi/dphi) + 0.5_work_p) * dphi
    call field_following_pos2_gpu(x,phi,phi_mid,xff)
    r=xff(1)
    z=xff(2)

    wp=w*ct(piw0)

    v(1)=1_work_p
    v(2)=epara
    v(3)=eperp
    v(4)=pot
    v(5)=sqrt((old_ph(pir)-new_ph(pir))**2+(old_ph(piz)-old_ph(piz))**2)

    ! for all sections
    do i=1, diag_heat_nsection
        ! check range
        if(diag_heat_rmin(i) < r .and. r < diag_heat_rmax(i) &
            .and. diag_heat_zmin(i) < z .and. z < diag_heat_zmax(i)) then

        !r index
        rn=(r-diag_heat_rmin(i))/diag_heat_dr(i)
        ir=floor(rn)+1
        if(ir<1 .or. diag_heat_nr<= ir) cycle
        wr=1.0_work_p - rn + real(ir-1,kind=work_p)

        !z index
        zn=(z-diag_heat_zmin(i))/diag_heat_dz(i)
        iz=floor(zn)+1
        if(iz<1 .or. diag_heat_nz<= iz) cycle
        wz=1.0_work_p - zn + real(iz-1,kind=work_p)

        ith0 = 1 + mod(ith,size(diag_heat_pv,6))
        lb = lbound(diag_heat_pv,1)
        ub = ubound(diag_heat_pv,1)

        do itmp=lb,ub
           i0 = (itmp-lb) + lbound(v,1)
           dval0 = v(i0)*wp*wr *wz
           dummy0 = atomicAdd( diag_heat_pv(itmp,ir,iz,i,stype,ith0), dval0 )

           dval1 = v(i0)*wp*(1.0_work_p-wr)*wz
           dummy1 = atomicAdd( diag_heat_pv(itmp,ir+1,iz,i,stype,ith0), dval1 )

           dval2 = v(i0)*wp*wr *(1.0_work_p-wz)
           dummy2 = atomicAdd( diag_heat_pv(itmp,ir,iz+1,i,stype,ith0), dval2 )

           dval3 = v(i0)*wp*(1.0_work_p-wr)*(1.0_work_p-wz)
           dummy3 = atomicAdd( diag_heat_pv(itmp,ir+1,iz+1,i,stype,ith0), dval3 )
        enddo
        !psi 
        psi=psi_interpol_gpu(r,z,0,0)
        pn=(psi-diag_heat_pmin(i))/diag_heat_dp(i)
        ip=floor(pn)+1
        if(ip<1 .or. diag_heat_npsi<= ip) cycle
        ws=1.0_work_p - pn + real(ip-1,kind=work_p)

        ith0 = 1 + mod(ith,size(diag_heat_pv_psi,5))
        lb = lbound(diag_heat_pv_psi,1)
        ub = ubound(diag_heat_pv_psi,1)

        do itmp=lb,ub
           i0 = (itmp-lb) + lbound(v,1)
           dval0 = v(i0)*wp*ws
           dummy0 = atomicAdd( diag_heat_pv_psi(itmp,ip,i,stype,ith0), dval0 )

           dval1 = v(i0)*wp*(1.0_work_p-ws)
           dummy1 = atomicAdd( diag_heat_pv_psi(itmp,ip+1,i,stype,ith0), dval1 )
        enddo
        endif
    enddo 
return
end subroutine diag_heat_port_gpu
!obtain derivatives of phase variable : actual calculation is done in derivs_sp. 
! prepare E-field and B-field
attributes(device) &
subroutine derivs_single_gpu(ptli,dy,i,time,fld,ith,diag_on,itr,p)

  use sml_module_gpu
  use fld_module, only : fld_type
  use ptl_module_gpu, only : ptl_ph_gpu,ptl_ct_gpu, type_gpu
  use precision_mod_gpu

  implicit none
  integer,intent(in) :: ptli, itr
  real (kind=work_p), intent(inout) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (kind=work_p), intent(in) :: time, p(3)
  type(fld_type), intent(inout) :: fld
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  !
  logical rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag) ! variables for diagnosis 

  integer, parameter :: idebug = 0

  ! Save space information
  fld%r=ptl_ph_gpu(ptli,1)
  fld%z=ptl_ph_gpu(ptli,2)
  fld%phi=ptl_ph_gpu(ptli,3)

  ! obtain B-field information 
  call field_gpu(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field for each particle : use position information from 'charge'
     call efield_gpu(i,fld,time,itr,p)
! if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp_elec_gpu(fld,time,ptli,dy,diag_on,vf_diag)
     
  else
     call remove_particle_gpu(i,-1)
! call remove_particle_gpu(sp%ptl(i)%gid, &
! sp%ptl(i)%ph, sp%ptl(i)%ct, &
! sp%type, i, -1)

! if (idebug >= 1) then
! print *, 'particle eliminated due to rz_outside', &
! i, sml_mype, sp%type, sp%ptl(i)%gid
! endif
  endif
     

  if(diag_on) call diag_1d_port1_gpu(ptli,dy,type_gpu,vf_diag,ith)

end subroutine derivs_single_gpu
    ! set minimum weight to prevent weight to be smaller than 0
    attributes(device) &
    subroutine restrict_weight_gpu(w)
      use precision_mod_gpu
      implicit none
      real (kind=work_p) :: w(2)
      real (kind=work_p), parameter :: weight_min = -50.0_work_p
      real (kind=work_p), parameter :: weight_max = 0.999_work_p

      w(1)=max(w(1),weight_min)
      w(2)=max(w(2),weight_min)

      w(1)=min(w(1),weight_max)
      w(2)=min(w(2),weight_max)


    end subroutine restrict_weight_gpu

attributes(device) &
subroutine derivs_single_with_e_gpu(ptli,dy,i,time,fld,E_mag)
  use sml_module_gpu, only : sml_mype, sml_n_vf_diag, sml_2pi
  use fld_module
  use ptl_module_gpu, only : ptl_ph_gpu, ptl_ct_gpu, type_gpu
  use precision_mod_gpu
  implicit none

  integer,intent(in) :: ptli
  real (kind=work_p) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (kind=work_p), intent(in) :: time 
  type(fld_type), intent(inout) :: fld
  real (kind=work_p), intent(in) :: E_mag(3)
! integer, intent(in) :: ith
  !
  logical :: rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)
  logical, parameter :: diag_on=.false.
  !
  real (kind=work_p) :: dpsi(2), E(3), bp, dtheta_norm(2), B
  real (kind=work_p) :: p(3), x(2), phi, phi_mid, xff(2)
  integer :: itr

  ! Save space information
  fld%r=ptl_ph_gpu(ptli,1)
  fld%z=ptl_ph_gpu(ptli,2)
  fld%phi=ptl_ph_gpu(ptli,3)

  ! obtain B-field information 
  call field_gpu(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field from previous E-field

     ! same E-field considering B-field curvature
     dpsi(1:2)=(/ fld%dpsidr, fld%dpsidz /)
     E(1:2)=E_mag(1)*dpsi

     bp=sqrt(fld%br**2+fld%bz**2)
     dtheta_norm(1:2)=(/ fld%br, fld%bz /)/bp
     E(1:2)=E(1:2) + E_mag(2)*dtheta_norm

     B=sqrt(fld%br**2+fld%bz**2+fld%bphi**2)
     E(3)=(E_mag(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi

     fld%Er=E(1)
     fld%Ez=E(2)
     fld%Ephi=E(3)
! if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp_elec_gpu(fld,time,ptli,dy,diag_on,vf_diag)
     
  else
     call remove_particle_gpu(i,-1)
! call remove_particle_gpu( sp%ptl(i)%gid, &
! sp%ptl(i)%ph, &
! sp%ptl(i)%ct, &
! sp%type, &
! i, -1 )

! print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  end if
end subroutine derivs_single_with_e_gpu
! single particle push -- get new_phase using rk4 with initial E-field
attributes(device) &
subroutine push_single_gpu(i,y,new_phase,dt,ith,diag_on,itr,p)
  use sml_module_gpu
  use ptl_module_gpu, only : ptl_gid_gpu, ptl_ph_gpu, piw1, piw2
  use fld_module
  use precision_mod_gpu
  !gpu use perf_monitor
  implicit none

  integer, intent(in) :: i, itr
  real (kind=work_p), intent(in) :: y(ptl_nphase), p(3)
  real (kind=work_p), intent(inout) :: new_phase(ptl_nphase)
  real (kind=work_p), intent(in) :: dt
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  !
  integer :: ptli
  type(fld_type) :: fld
  real (kind=work_p) :: dy(ptl_nphase),dyt(ptl_nphase),dym(ptl_nphase)
  integer :: rtn,j,i1
  logical :: rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)
  real (kind=work_p) , external :: psi_interpol
  !gpu character (len=5) :: err_str(0:1)
  real (kind=work_p) :: hdt, time, th ! -- time should be input 
  real (kind=work_p) :: bp,B, E_mag(3)

  !gpu err_str(1)='ion'
  !gpu err_str(0)='elec'


  !pushing particle "i"


  time=sml_time ! need to update -- rk4 time ###

  hdt=dt*0.5D0
  th=time + hdt

  if(ptl_gid_gpu(i)>0) then
     
     !set ptli -- ptli%ct does not change
! ptli%ct=sp%ptl(i)%ct
! ptli%ph=sp%ptl(i)%ph
     ptli = i 
     !get derivs with updating E-field information - assigned on fld
     !diag-ports are called when diag_on is .true.
     call derivs_single_gpu(ptli,dy,i,time,fld,ith,diag_on,itr,p)




     ! RK2 - RK4 hybrid -- time advance with RK4 with time-constant E-field
     
     ! get E-field in magnetic field
     bp=sqrt(fld%br**2 + fld%bz**2) 
     B=sqrt(bp**2 + fld%Bphi**2 )
     E_mag(2)=(fld%Er*fld%Br + fld%Ez*fld%Bz)/bp
     E_mag(3)=(E_mag(2)*bp + fld%Ephi*fld%Bphi)/B ! parallel field
     E_mag(1)=(fld%Er*fld%dpsidr + fld%Ez*fld%dpsidz )/(fld%dpsidr**2 + fld%dpsidz**2)
     
     ! get derivs with existing E-field

     ptl_ph_gpu(ptli,:)=y
!write(*,*)'y = ', y
     call derivs_single_with_e_gpu(ptli,dy ,i,time ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + hdt * dy 
     call derivs_single_with_e_gpu(ptli,dyt,i,th ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + hdt * dyt
     call derivs_single_with_e_gpu(ptli,dym,i,th ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + dt * dym
     dym = dyt + dym
     call derivs_single_with_e_gpu(ptli,dyt,i,time+dt,fld,E_mag)
     
     ! Obtain new_phase
     new_phase = y + dt/6D0 * ( dy + dyt + 2D0*dym )
! call restrict_weight_gpu(new_phase(piw1:piw2))
  endif
    
end subroutine push_single_gpu
! single particle push -- get new_phase using rk4 with initial E-field
attributes(device) &
subroutine pushe_single_gpu(i,y,new_phase,dt,ith,diag_on,itr,p)
  use sml_module_gpu
  use ptl_module_gpu, only : ptl_gid_gpu, ptl_ph_gpu, piw1, piw2
  use fld_module
  use precision_mod_gpu
  !gpu use perf_monitor
  implicit none

  integer, intent(in) :: i, itr
  real (kind=work_p), intent(in) :: y(ptl_nphase), p(3)
  real (kind=work_p), intent(inout) :: new_phase(ptl_nphase)
  real (kind=work_p), intent(in) :: dt
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  !
  integer :: ptli
  type(fld_type) :: fld
  real (kind=work_p) :: dy(ptl_nphase),dyt(ptl_nphase),dym(ptl_nphase)
  integer :: rtn,j,i1
  logical :: rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)
  real (kind=work_p) , external :: psi_interpol
  !gpu character (len=5) :: err_str(0:1)
  real (kind=work_p) :: hdt, time, th ! -- time should be input 
  real (kind=work_p) :: bp,B, E_mag(3)

  !gpu err_str(1)='ion'
  !gpu err_str(0)='elec'


  !pushing particle "i"


  time=sml_time ! need to update -- rk4 time ###

  hdt=dt*0.5D0
  th=time + hdt

  if(ptl_gid_gpu(i)>0) then
     
     !set ptli -- ptli%ct does not change
! ptli%ct=sp%ptl(i)%ct
! ptli%ph=sp%ptl(i)%ph
     ptli = i 
     !get derivs with updating E-field information - assigned on fld
     !diag-ports are called when diag_on is .true.
     call derivs_single_gpu(ptli,dy,i,time,fld,ith,diag_on,itr,p)


     

     ! RK2 - RK4 hybrid -- time advance with RK4 with time-constant E-field
     
     ! get E-field in magnetic field
     bp=sqrt(fld%br**2 + fld%bz**2) 
     B=sqrt(bp**2 + fld%Bphi**2 )
     E_mag(2)=(fld%Er*fld%Br + fld%Ez*fld%Bz)/bp
     E_mag(3)=(E_mag(2)*bp + fld%Ephi*fld%Bphi)/B ! parallel field
     E_mag(1)=(fld%Er*fld%dpsidr + fld%Ez*fld%dpsidz )/(fld%dpsidr**2 + fld%dpsidz**2)
     
     
     ! get derivs with existing E-field
     ptl_ph_gpu(ptli,:)=y
     call derivs_single_with_e_gpu(ptli,dy ,i,time ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + hdt * dy 
     call derivs_single_with_e_gpu(ptli,dyt,i,th ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + hdt * dyt
     call derivs_single_with_e_gpu(ptli,dym,i,th ,fld,E_mag)
     
     ptl_ph_gpu(ptli,:) = y + dt * dym
     dym = dyt + dym
     call derivs_single_with_e_gpu(ptli,dyt,i,time+dt,fld,E_mag)
     
     ! Obtain new_phase
     new_phase = y + dt/6D0 * ( dy + dyt + 2D0*dym )
! call restrict_weight_gpu(new_phase(piw1:piw2))
  endif
    
end subroutine pushe_single_gpu
attributes(device) &
subroutine bounce_gpu(new_phase,old_phase,rtn)
  use eq_module_gpu , only : eq_axis_r,eq_x_psi
  use ptl_module_gpu, only : ptl_nphase, piw1, piw2
  use sml_module_gpu
  use bnc_module_gpu, only : bnc_max_r, bnc_min_r
  use precision_mod_gpu
  implicit none
  real (kind=work_p),intent(inout) :: new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer, intent(out) :: rtn
  integer :: inner
  real (kind=work_p) :: r, z, psi,sign
  real (kind=work_p) :: b , br, bz, bphi,mid,z1,z2,z_pmin
  integer :: i,j,count
  integer, parameter :: JMAX=5
  !real (kind=8), external :: psi_interpol,z_psi_min,I_interpol,B_interpol_sym
  ! real (kind=8), parameter :: ZTOL=1D-6, NPSI_TOL=1D-5
  ! real (kind=8), parameter :: ZTOL=1D-10, NPSI_TOL=1D-5
  real (kind=work_p), parameter :: ZTOL=1E-6_work_p, NPSI_TOL=1E-4_work_p, BTOL=1E-8_work_p, PTOL=1E-8_work_p
  real (kind=work_p) :: psitol
  real (kind=work_p) :: delz,deltab2,deltap2
  real (kind=work_p) :: RBt,oldB,oldP,newB,newP,deltab,deltap,sign2,deltar,deltaz,zh,rh,dz1,dz2,dz3
  real (kind=work_p) :: dpsi_dr,dpsi_dz,d2psi_d2r,d2psi_drdz,d2psi_d2z,db_dr,db_dz,denom
  logical :: use_current, diverge
  psitol=NPSI_TOL*eq_x_psi
  rtn=0

  ! reverse angle -- z -> -z, phase variable -> old phase variable except z
  !if(sml_concentric) then
  if(.true.) then
     psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
! if(psi<=sml_inpsi .or. sml_bounce==2 .and. psi>= sml_outpsi) then
     if( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi ) then
        new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
        new_phase(2)=-old_phase(2)
        if(sml_bounce_zero_weight==1) new_phase(piw1:piw2)=0.0_work_p
     endif
  else 
     ! bounce at inner (& outer )boundary
     ! Finding z value which has same psi and r of old particle position
     ! This position gives same mu, P_phi.
     ! similar energy if B is similar.
     ! for exact calculation, 1st order correction is needed.
     psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
     r=old_phase(1) ! radius of old position
     ! if((psi <= sml_inpsi .or. sml_bounce==2 .and. psi >= sml_outpsi) .and. &
     ! r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then 
     if(( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi) .and. &
          r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then 
        rtn=1 ! a big jump will be happen
        if(is_rgn1_gpu(r,z,psi)) then
           inner=1
        else
           inner=0
        endif
        psi=psi_interpol_gpu(old_phase(1),old_phase(2),0,0) ! psi is old position value
        ! if(psi > (sml_inpsi-psitol) .and. psi < (sml_outpsi+psitol)) then !old psi is inside simulation range
        if(psi >= (sml_inpsi-psitol) .and. psi <=(sml_outpsi+psitol)) then !old psi is inside simulation range
           z=old_phase(2)
           z_pmin=z_psi_min_gpu(r)

           ! bug fixed. sign should be -1 
           ! if(r>eq_axis_r) then
           ! sign=1
           ! else
           sign=-1
           ! endif

           ! set boundary of binary search
           z1=z_pmin 
           if(z>z_pmin) then 
              ! z2=max(z_pmin - (z-z_pmin+0.01)*2. , sml_bd_min_z) ! heuristic formula
              if( z_pmin- (z-z_pmin+0.01)*2. > sml_bd_min_z ) then
                 z2= z_pmin- (z-z_pmin+0.01)*2.
              else
                 z2=sml_bd_min_z
              endif
           else
              ! z2=min(z_pmin - (z-z_pmin+0.01/sml_norm_r)*2. , sml_bd_max_z) ! heuristic formula
              if(z_pmin - (z-z_pmin+0.01)*2. < sml_bd_max_z) then
                 z2=z_pmin - (z-z_pmin+0.01)*2
              else
                 z2=sml_bd_max_z
              endif
           endif

           ! find z value using binary search.-------------------------
           do while (abs(z1-z2) > ZTOL)
              mid=0.5_work_p*(z1+z2)
              if(sign *(psi_interpol_gpu(r,mid,0,0)-psi) < 0 ) then
                 z2=mid
              else
                 z1=mid
              endif
              ! write(400,*) mid,z1,z2,z1-z2,psi,psi_interpol(r,mid,0,0)
           enddo

           ! z1=0.5D0*(z1+z2)
           if(inner==1) then ! z1 gives larger psi for inner bd.
              z1=z2
           endif
           !------------------------------------------------------------

           new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
           new_phase(2)=z1 

           ! 1st order correction 
           ! delz=(psi-new_phase(9))/psi_interpol(new_phase(1),z1,0,1)
           ! new_phase(2)=new_phase(2)+delz
           ! new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)

           ! (dR,dZ) correction for attaining position with same B and psi values as old position 
           oldB=B_interpol_sym_gpu(old_phase(1),old_phase(2))
           oldP=psi_interpol_gpu(old_phase(1),old_phase(2),0,0)
           deltab2=oldB-B_interpol_sym_gpu(new_phase(1),new_phase(2))
           deltap2=oldP-psi_interpol_gpu(new_phase(1),z1,0,0)

           !init
           if( .not. (z > sml_bd_min_z .and. z < sml_bd_max_z) ) then
              rtn=-1
              return
           endif
           use_current=.false.
           diverge=.false.
           do j=1, JMAX
              ! Newton-Raphson procedure is iterated 
              psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol_gpu(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol_gpu(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol_gpu(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              if(((dabs(deltab)/oldB)<BTOL).and.((dabs(deltap)/oldP)<PTOL)) then 
                 use_current=.true.
                 exit
              endif
              
              d2psi_d2r=psi_interpol_gpu(new_phase(1),new_phase(2),2,0)
              d2psi_drdz=psi_interpol_gpu(new_phase(1),new_phase(2),1,1)
              d2psi_d2z=psi_interpol_gpu(new_phase(1),new_phase(2),0,2)
              
              db_dr=-newB/new_phase(1)+1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_d2r+dpsi_dz*d2psi_drdz)
              db_dz=1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_drdz+dpsi_dz*d2psi_d2z)
              
              denom=dpsi_dr*db_dz-dpsi_dz*db_dr
              deltar=(db_dz*deltap-dpsi_dz*deltab)/denom
              deltaz=(-db_dr*deltap+dpsi_dr*deltab)/denom
              
              ! move to new (r,z) 
              new_phase(1)=new_phase(1)+deltar
              new_phase(2)=new_phase(2)+deltaz
              
              ! check diverge
              if(new_phase(1) < sml_bd_min_r .or. new_phase(1) > sml_bd_max_r &
                   .or. new_phase(2) < sml_bd_min_z .or. new_phase(2) > sml_bd_max_z ) then
                 diverge=.true.
                 use_current=.false. ! for safety
                 exit
              endif
           enddo
           
           if(.NOT. use_current .and. .NOT. diverge ) then
              ! loop ended not satisfying the TOL codition, nor diverge
              psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol_gpu(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol_gpu(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol_gpu(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              ! use original (binary search position)
              if(((deltab2/oldB)**2+(deltap2/oldP)**2)<((deltab/oldB)**2+(deltap/oldP)**2)) then
                 use_current=.false.
              endif
           endif
           
           if(.NOT. use_current) then 
              new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
              new_phase(2)=z1
           endif
           
           psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
           ! end of second (dR,dZ) correction for same B

           if(psi < (sml_inpsi-psitol) .or. psi > (sml_outpsi+psitol)) then
              ! if(new_phase(9) <= sml_inpsi-psitol .or. new_phase(9) >= sml_outpsi+psitol) then
              !print *, 'Fail finding proper psi', psi/eq_x_psi, oldP/eq_x_psi, psi/eq_x_psi
              !print *, 'oldB, newB', oldB, newB
              rtn=-1
           endif
           if(sml_bounce_zero_weight==1) new_phase(piw1:piw2)=0D0 
        else ! old position was outside of boundary
           rtn=-1
           ! print *, 'ptl_elimination', i, sml_mype
        endif
     endif
  endif

end subroutine bounce_gpu

attributes(device) &
  subroutine derivs_gpu(x,phi,dx)
    use sml_module_gpu
    use precision_mod_gpu

    implicit none
    real (kind=work_p), intent(in) :: x(2), phi
    real (kind=work_p), intent(out) :: dx(2)
    real (kind=work_p) :: b(2), bphi, r,z

    r=min(max(x(1),sml_bd_min_r),sml_bd_max_r)
    z=min(max(x(2),sml_bd_min_z),sml_bd_max_z)


    call bvec_interpol_gpu(r,z,phi,b(1),b(2),bphi)
    dx = b/bphi*x(1)

  end subroutine derivs_gpu


attributes(device) &
subroutine field_following_pos2_gpu(x_org,phi_org,phi_dest,x_dest)
  use precision_mod_gpu
  implicit none
  real (kind=work_p), intent(in) :: x_org(2),phi_org, phi_dest
  real (kind=work_p), intent(out) :: x_dest(2)
  real (kind=work_p) :: phi, x(2),dphi
  real (kind=work_p) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)
  integer :: sml_ff_step=1
  integer :: i

  dphi=(phi_dest-phi_org)/real(sml_ff_step)

  ! 0 step
  phi=phi_org
  x=x_org

  do i=1, sml_ff_step
     ! get first derivative
     call derivs_gpu(x,phi,dx1)

     if( .false. ) then ! first order calculation
        x = x + dx1*dphi
     else if( .true. ) then ! second order calculation - rk2 
        ! obtain mid point
        hh=dphi*0.5D0

        x_tmp = x + dx1*hh

        ! get new derivative
        call derivs_gpu(x_tmp,phi+hh,dx2)

        ! advance one step using mid-point derivative
        x = x + dx2*dphi
     else
        ! 4th Order Calculation - rk4
        !
        hh=dphi*0.5D0
        h6=dphi/6D0
        ! derivative 1 (from x)- obtained already
        x_tmp=x + hh*dx1 ! yt=y+hh*dydx : yt -> x_tmp
        !
        ! derivative 2 (from x_tmp) 
        call derivs_gpu(x_tmp,phi+hh,dx2)! dyt from yt : dyt -> dx2
        x_tmp=x + hh*dx2 ! yt=y+hh*dyt : yt -> x_tmp
        !
        ! derivative 3 (from x_tmp)
        call derivs_gpu(x_tmp,phi+hh,dx3)! dym from yt : dym -> dx3 
        x_tmp=x + dphi*dx3 ! yt=y + h*dym : yt -> x_tmp
        dx3 = dx2 + dx3 ! dym = dyt + dym : dym -> dx3 , dyt -> dx2
        !
        ! derivative 4 (from x_tmp)
        call derivs_gpu(x_tmp,phi+dphi,dx2) ! dyt from yt : dyt -> dx2, yt -> x_tmp 
        x = x + h6 * (dx1 + dx2 + 2D0*dx3) ! yout = y + h6* (dydx+dyt+2D0*dym) 
     endif
     phi=phi+dphi
  enddo
  x_dest=x

end subroutine field_following_pos2_gpu





  attributes(device) &
  subroutine t_coeff_gpu(itr,x,p)
  use grid_class_gpu, only : grid_mapping
  use precision_mod_gpu
    implicit none
    integer, intent(in) :: itr
    real (kind=work_p), intent(in) :: x(2)
    real (kind=work_p), intent(out) :: p(3)

    integer :: nd
    real (kind=work_p) :: dx(2)

    ! dx=x - grid%x(:,grid%nd(3,itr))
    dx(1:2) = x(1:2) - grid_mapping(1:2,3,itr)
    p(1:2)= grid_mapping(1:2,1,itr)*dx(1) + grid_mapping(1:2,2,itr)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)


  end subroutine t_coeff_gpu





      attributes(host) &
      subroutine init_push_mod_gpu( grid )

      use grid_class, only : grid_type


      use grid_class_gpu, only : &
          update_device_grid_type

      implicit none
      type(grid_type) :: grid

      ! --------------------------------------------------------
      ! perform one-time initialization related to push() on GPU
      ! --------------------------------------------------------
      call push_update_device_gpu()

      call update_device_grid_type(grid )



      return
      end subroutine init_push_mod_gpu
  !!return initial guess value for triangle search
  attributes(device) &
  subroutine guess_gpu(x,init)
    use precision_mod_gpu
    use grid_class_gpu, only :grid_guess_min, grid_inv_guess_d, grid_guess_n, grid_guess_table 
    implicit none
    real (kind=work_p), intent(in) :: x(2)
    integer, intent(out) :: init
    integer :: i(2)

    i= (x-grid_guess_min)*grid_inv_guess_d +1
    !error message for debug only
! if(i(1)<=0 .or. i(1)>grid%guess_n(1)) then
! print *, 'Invaild number for guess table- R',i(1),x(1),grid%guess_min(1),grid%guess_max(1)
! endif
! if(i(2)<=0 .or. i(2)>grid%guess_n(2)) then
! print *, 'Invaild number for guess table- Z',i(2),x(2),grid%guess_min(2),grid%guess_max(2)
! endif

    i=min(max(i,1),grid_guess_n)
    init=grid_guess_table(i(1),i(2))

  end subroutine guess_gpu

attributes(device) &
subroutine search_tr2_gpu( xy, itr, p )
  use grid_class_gpu, only : grid_guess_table, grid_guess_min, grid_inv_guess_d, &
                             grid_guess_xtable, grid_guess_count, grid_guess_list 
  use precision_mod_gpu
  implicit none
  real(kind=work_p) :: xy(2)
  integer :: itr
  real(kind=work_p) :: p(3)

  real(kind=work_p), parameter :: zero = 0.0d0
  real(kind=work_p), parameter :: eps = 10.0d0*epsilon(zero)
  integer :: ij(2), istart,iend, k, itrig
  integer :: i,j, ilo,ihi, jlo,jhi
  real(kind=work_p) :: dx(2), pmin, pmax, dp
  logical :: is_found


  ilo = lbound( grid_guess_table, 1 )
  ihi = ubound( grid_guess_table, 1 )

  jlo = lbound( grid_guess_table, 2 )
  jhi = ubound( grid_guess_table, 2 )

  ij = (xy - grid_guess_min)*grid_inv_guess_d + 1
  i = max(ilo, min(ihi, ij(1)) )
  j = max(jlo, min(jhi, ij(2)) )


  istart = grid_guess_xtable(i,j)
  iend = istart + grid_guess_count(i,j) - 1


  itr = -1
  do k=istart,iend
     itrig = grid_guess_list(k)
     ! call t_coeff( grid, itrig, xy, p )

    dx(1:2) = xy(1:2) - grid_mapping(1:2,3,itrig)
    p(1:2)= grid_mapping(1:2,1,itrig)*dx(1) + &
            grid_mapping(1:2,2,itrig)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)

     if (minval(p) .ge. -eps) then
        itr = itrig
        exit
     endif
  enddo

  return
end subroutine search_tr2_gpu



      subroutine push_update_host_gpu( sp, &
                                       psn, diag_on, &
                                       gpu_ibegin, gpu_iend )
      use sml_module, only : sml_neutral_use_ion_loss
      use psn_class, only : psn_type
      use ptl_module, only : species_type
      use ptl_module_gpu, only : &
         update_host_species_type

      use diag_module_gpu, only : &
          update_host_diag 

      use psn_class_gpu, only : &
          update_host_psn_type

      use neu_module_gpu, only : &
          update_host_neu

      implicit none
      logical, intent(in) :: diag_on
      integer, intent(in) :: gpu_ibegin, gpu_iend

      type(psn_type) :: psn
      type(species_type) :: sp


      integer, parameter :: idebug = 0
      ! need to copy 
      ! ptl(1:sp%num)%ph 
      ! ptl(1:sp%num)%ct
      ! particle-triangle data %itr, %p(3) for E-field calculation
      ! psn%E_rho_ff
      !
      if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_host_sml_type '
! if (sml_mype == 0) print *,'diag_on = ', diag_on
      endif
      ! call update_host_sml()

      if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_host_psn_type '
      endif
      call update_host_psn_type(psn)

      if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_host_neu '
      endif
      if(.not. sml_neutral_use_ion_loss) then
        call update_host_neu()
      endif


      if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_host_species_type '
      endif
      call update_host_species_type(sp, &
                                    gpu_ibegin,gpu_iend)


! call update_host_grid_type()

      !if (diag_on) then -- diag_on is true only when ipc==1
        if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_host_diag '
        endif
        call update_host_diag()
      !endif


      return
      end subroutine push_update_host_gpu
attributes(host) &
subroutine push_update_device_gpu( )
 use sml_module, only : sml_mype
! use bicub_mod, only : psi_bicub

 use sml_module_gpu, only : update_device_sml 
 use eq_module_gpu, only : update_device_eq 
 use itp_module_gpu, only : update_device_itp 
 use one_d_cub_mod_gpu, only : update_device_one_d 
 use bicub_mod_gpu, only : update_device_bicub 
 use ptl_module_gpu, only : update_device_ptl 
 use diag_module_gpu, only : update_device_diag 
 use bnc_module_gpu, only : update_device_bnc 

 integer, parameter :: idebug = 0

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_sml()'
 endif
! call update_device_sml()

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_eq()'
 endif
 call update_device_eq()

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_itp()'
 endif
 call update_device_itp()

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_one_d()'
 endif
 call update_device_one_d()

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_bicub()'
 endif
 call update_device_bicub()

 if (idebug >= 1) then
! if (sml_mype == 0) print*,'before update_device_ptl()'
 endif
! call update_device_ptl( )

 if (idebug >= 1) then
! if (sml_mype) print*,'before update_device_diag()'
 endif
 call update_device_diag()

 if (idebug >= 1) then
! if (sml_mype) print*,'before update_device_bnc()'
 endif
 call update_device_bnc()

 return
end subroutine push_update_device_gpu
attributes(host) &
subroutine pushe_gpu(istep,ihybrid,ncycle,grid,psn,sp,diag_on_input)
  use grid_class
  use psn_class
  use ptl_module
  use perf_monitor
  use sml_module
  use omp_module , only : split_indices
  use eq_module
  use cudafor
  use precision_mod_gpu

  use ptl_module_gpu, only : &
     phase0_gpu, ptl_ph_gpu, &
     ptl_ct_gpu, &
! iperm_gpu, xstart_gpu, &
     update_device_species_type, &
     update_host_species_type

  use psn_class_gpu, only : &
    update_device_psn_type

! use bicub_mod_gpu, only : &
! update_device_bicub

  use gen_perm_gpu_mod, only : &
    gen_perm_gpu

  use reorder_gpu_mod, only : &
    reorder1d_gpu, &
    reorder2d_gpu

  use util_mod_gpu, only : get_gpu_streamid

  implicit none
  integer, intent(in) :: istep, ihybrid
  integer, intent(in) :: ncycle
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type),target :: sp
  logical :: diag_on,need_copyin,need_copyout,diag_on_input
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  integer :: istatus
  real(kind=work_p) :: gpu_ratio
  integer :: perm_gpu_freq
  logical :: perm_gpu
  character(len=255) :: sml_gpu_ratio_str
  integer :: gpu_ibegin, gpu_iend
  integer :: cpu_ibegin, cpu_iend

  integer :: ierr,stream1
  integer :: nblocks
  type(dim3) :: tgrid, tblock
  integer :: lb1,ub1,lb2,ub2
    logical :: rz_outside
! type(ptl_type), dimension(:), pointer :: ptl
 
  real (kind=work_p) :: dt_now,dt,time_now,new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer, parameter :: idebug =0 
  real(kind=work_p), parameter :: dzero = 0.0_work_p

  integer :: icycle
  integer :: epc, i, rtn

  real (kind=work_p) :: phi_mid, x(2), phi, mu, rho, xff(2)
  real (kind=work_p) :: p(3)
  integer :: itr, ip
! real(kind=work_p) :: phase_tmp(ptl_nphase)

  attributes(device) :: d_ibegin,d_iend !, phase_tmp


! -------------------------------------------
! local variables related to particle sorting
! -------------------------------------------
  logical, parameter :: use_sort_particles = .true.
  logical, parameter :: use_reorder_array = .false.
  integer :: ilo,ihi,jlo,jhi
  integer :: nx, ny, n, xydim
  real(kind=work_p) :: xmin,ymin, inv_dx,inv_dy
  integer*8, allocatable, dimension(:) :: gid
  real(kind=work_p), allocatable, dimension(:) :: xcoord,ycoord 
  real(8), allocatable, dimension(:,:) :: tphase0
  type(ptl_type), allocatable, dimension(:) :: tptl
  integer, allocatable, dimension(:) :: iperm
  integer, allocatable, dimension(:) :: iperm_gpu, xstart_gpu
  integer :: itotal
  integer :: mm, nn, lld
  integer :: streamid
  attributes(device) :: d_ibegin,d_iend, iperm_gpu, xstart_gpu

  if(sp%num==0) return ! nothing to push

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

  if (use_sort_particles) then
    call t_startf("pushe_sort_particles")
! --------------
! sort particles
! --------------
    ilo = lbound( grid%guess_table,1)
    ihi = ubound( grid%guess_table,1)
    jlo = lbound( grid%guess_table,2)
    jhi = ubound( grid%guess_table,2)

    allocate( iperm(sp%num))
    allocate( gid(sp%num), xcoord(sp%num), ycoord(sp%num) )

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
         gid(i) = sp%ptl(i)%gid
         xcoord(i) = sp%ptl(i)%ph(1)
         ycoord(i) = sp%ptl(i)%ph(2)

         iperm(i) = i
      enddo
    enddo

    call gen_perm( ilo,ihi,jlo,jhi, &
             grid%guess_min, grid%inv_guess_d, &
             sp%num, gid, xcoord, ycoord, iperm )

! write(20+sml_mype,*) 'iperm(1:5) ', iperm(1:5)
! write(20+sml_mype,*) 'iperm(n-5:n) ', iperm(sp%num-1:sp%num)


    deallocate( gid, xcoord, ycoord )
! ------------------------
! rearrange particle data in sp
! 
! which arrays should be permutated
! ptl
!
! should other arrays such as phase0(:) or rhoi(:) 
! be rearranged????y
! ------------------------


! --------------
! permute sp%ptl
! --------------
    allocate(tptl(sp%num))
!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
        tptl(i) = sp%ptl( iperm(i) )
      enddo
    enddo

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
      sp%ptl(i) = tptl(i)
     enddo
    enddo
    deallocate( tptl )



! -----------------
! permute sp%phase0
! -----------------
    allocate(tphase0(lbound(sp%phase0,1):ubound(sp%phase0,1),sp%num))

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
        tphase0(:,i) = sp%phase0(:, iperm(i) )
      enddo
    enddo

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
      sp%phase0(:,i) = tphase0(:,i)
     enddo
    enddo
    deallocate( tphase0 )



    deallocate( iperm )
    call t_stopf("pushe_sort_particles")
   endif

! ptl => sp%ptl
  ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)

  diag_on=diag_on_input

  if ((sml_gpu_ratio < 0.0) .or. (sml_gpu_ratio > 1.0)) then
     sml_gpu_ratio = 0.7
     call getenv("SML_GPU_RATIO",sml_gpu_ratio_str)
     if (len(trim(sml_gpu_ratio_str)) >= 1) then
        read(sml_gpu_ratio_str,*,iostat=istatus) gpu_ratio
        if (istatus.eq.0) then
           gpu_ratio = max( 0.0_work_p, min(1.0_work_p, dble(gpu_ratio)))
           sml_gpu_ratio = gpu_ratio
           if (idebug >= 1) then
              if (sml_mype.eq.0) write(*,*) 'sml_gpu_ratio ', sml_gpu_ratio
           endif
        endif
     endif
  endif

! -----------------------------------------------
! suggest always choose gpu_ibegin = 1
! so that transfer to GPU is in contiguous block
! -----------------------------------------------
  gpu_ibegin = 1
  gpu_iend = int( sml_gpu_ratio*sp%num )
! gpu_iend = min( gpu_iend, ubound(ptl_ph_gpu,1))
  cpu_ibegin = gpu_iend + 1
  cpu_iend = sp%num

! print *, 'gpu_ibegin', cpu_ibegin,cpu_iend
  dt=sml_dt*0.5_work_p/real(sml_ncycle_half)

if (gpu_ibegin <= gpu_iend .or. cpu_ibegin <= cpu_iend) then
 if(gpu_ibegin <= gpu_iend) then

  call t_startf("push_upd_dev_gpu")
! call push_update_device_gpu()
!if(sml_mype==0) print *, 'before update_device_ptl'
  call update_device_ptl( )
!if(sml_mype==0) print *, 'before update_device_diag'
  call update_device_diag()
!if(sml_mype==0) print *, 'after update_device_diag'
  call t_stopf("push_upd_dev_gpu")

  call t_startf("upd_dev_psn")
!print *, 'before update_device_psn_type'
  call update_device_psn_type(psn )
!print *, 'after update_device_psn_type'
  call t_stopf("upd_dev_psn")

  call t_startf("upd_dev_species_type")
! print *, 'before update_device_species_type'
  call update_device_species_type(sp, gpu_ibegin,gpu_iend)
  call update_device_sml()
  call update_device_neu()
  call t_stopf("upd_dev_species_type")
! print *, 'gpu_ibegin', cpu_ibegin,cpu_iend
  xmin = grid_guess_min(1)
  ymin = grid_guess_min(2)
  
  inv_dx = grid_inv_guess_d(1)
  inv_dy = grid_inv_guess_d(2)

  nx = size(grid_guess_table,1)
  ny = size(grid_guess_table,2)

  xydim = size(ptl_ph_gpu,1)
  mm = gpu_iend - gpu_ibegin + 1

  ! -------------------------
  ! launch kernel computation
  ! -------------------------
  tblock%x = block_dim
  tblock%y = 1
  tblock%z = 1
  tgrid%x = grid_dim
  tgrid%y = 1
  tgrid%z = 1

  if (diag_on) then
    call setval_gpu( size(diag_1d_f_pv1),diag_1d_f_pv1,dzero)
    if(diag_eflux_on) call setval_gpu( size(diag_1d_eflux_pv),diag_1d_eflux_pv,dzero)

    if (sml_deltaf) then
      call setval_gpu( size(diag_1d_df_pv1),diag_1d_df_pv1,dzero)
      if(diag_eflux_on) call setval_gpu( size(diag_2d_dflux_pv),diag_2d_dflux_pv,dzero)
    endif
  endif

  if(diag_heat_on) call setval_gpu(size(diag_heat_pv),diag_heat_pv,dzero)
  if(diag_heat_on) call setval_gpu(size(diag_heat_pv_psi),diag_heat_pv_psi,dzero)

  if (idebug >= 1) then
     if (sml_mype.eq.0) write(*,*) 'size diag_1d_f_pv1', diag_1d_npv1, nx, ny, inv_dx, inv_dy
  endif
 endif

  if (sml_perm_gpu_freq < 1) then
     perm_gpu = .false.
     perm_gpu_freq = ncycle
  else
     perm_gpu = .true.
     perm_gpu_freq = sml_perm_gpu_freq
  endif

 do icycle=1, ncycle
!
     if(gpu_ibegin <= gpu_iend) then
      if ((perm_gpu) .and. (mod(icycle,perm_gpu_freq).eq.0)) then

        call t_startf("gen_perm_gpu")
        allocate(iperm_gpu(gpu_ibegin:gpu_iend),xstart_gpu(nx*ny+1), stat=ierr)
        call assert(ierr.eq.0,'alloc(iperm_gpu) ',ierr)

        call gen_perm_gpu( nx, ny, xmin, ymin, &
                      inv_dx, inv_dy, &
                      gpu_iend-gpu_ibegin+1, ptl_gid_gpu, ptl_ph_gpu(:,1:2), xydim, iperm_gpu, xstart_gpu )
        call t_stopf("gen_perm_gpu")

        call t_startf("reordering_gpu")

        nn = size( ptl_ph_gpu,2)
        lld = size( ptl_ph_gpu,1)

        if (idebug >= 1) then
         if (sml_mype.eq.0) write(*,*) 'size ptl_ph_gpu',size( ptl_ph_gpu,2), size( ptl_ph_gpu,1)
        endif

        call reorder2d_gpu( mm,nn, ptl_ph_gpu, lld, iperm_gpu )

        call reorder1d_gpu( mm, ptl_gid_gpu, iperm_gpu )

        nn = size( ptl_ct_gpu,2)
        lld = size( ptl_ct_gpu,1)

        if (idebug >= 1) then
         if (sml_mype.eq.0) write(*,*) 'size ptl_ct_gpu',size( ptl_ct_gpu,2), size( ptl_ct_gpu,1)
        endif

        call reorder2d_gpu( mm,nn, ptl_ct_gpu, lld, iperm_gpu )

        nn = size( phase0_gpu,2)
        lld = size( phase0_gpu,1)
        call reorder2d_gpu( mm,nn, phase0_gpu, lld, iperm_gpu )

        if (idebug >= 1) then
         if (sml_mype.eq.0) write(*,*) 'size phase0_gpu',size( phase0_gpu,2), size( phase0_gpu,1)
        endif

        deallocate(iperm_gpu,xstart_gpu)

        call t_stopf("reordering_gpu")
      endif
     endif
     ! 
     do epc=1,sml_nrk
        call t_startf("electron_loop")

        sml_epc=epc
    
        select case(epc)
        case(1)
           dt_now=0.5_work_p*dt
        case(2)
           dt_now=dt
        end select

        if(ihybrid>1 .or. icycle >1 .or. epc > 1) diag_on=.false.

! if(gpu_ibegin <= gpu_iend) then
        call t_startf("gpu_processing")

        call t_startf("cuda_thrd_sync1")
        istatus = cudaThreadSynchronize()
        call t_stopf("cuda_thrd_sync1")

        if (idebug >= 1) then
! if (sml_mype.eq.0) write(*,*) 'before kernel', sp%ptl(5)%ph, sp%ptl(5)%ct
        endif

        if (idebug >= 1) then
! if (sml_mype == 0) write(*,*) 'before launch kernel',icycle, sizeof(phase0_gpu),istatus
        endif

        streamid = get_gpu_streamid()

        call t_startf("pushe_kernel_gpu")
        call pushe_kernel_gpu<<<tgrid,tblock,0,streamid>>>(istep, &
              epc,phase0_gpu,diag_on,dt_now, &
              gpu_ibegin,gpu_iend)
        call t_stopf("pushe_kernel_gpu")

        if (idebug >= 1) then
! if (sml_mype == 0) write(*,*) 'after launch kernel',icycle, epc, istep
        endif

        call t_stopf("gpu_processing")

        call t_startf("cpu_processing")

! endif
! --------------------------------------
! concurrently push particles on the cpu
! --------------------------------------
        if (cpu_ibegin <= cpu_iend) then
! if (idebug >= 1) then
! write(*,*) 'cpu_ibegin',cpu_ibegin, cpu_iend, sml_mype
! endif

           call split_indices(cpu_iend-cpu_ibegin+1, sml_nthreads, i_beg, i_end)
           i_beg = i_beg + (cpu_ibegin-1)
           i_end = i_end + (cpu_ibegin-1)

           if (idebug >= 1) then
              if (sml_mype.eq.0) write(*,*) 'before kernel',sp%ptl(5)%ph, sp%ptl(5)%ct
           endif

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, X, PHI, PHI_MID, &
!$OMP XFF, ITR, P, &
!$OMP NEW_PHASE, OLD_PHASE, RTN )
           do ith=1,sml_nthreads

              call t_startf("PUSHE_CPU_LOOP")
              do i=i_beg(ith),i_end(ith)
                 if(sp%ptl(i)%gid>0) then

                    ! get proper toroidal angle index and weight
                    x=sp%ptl(i)%ph(1:2)
                    phi=sp%ptl(i)%ph(3)
                    phi_mid=(floor(phi/grid%delta_phi) + 0.5_work_p) * grid%delta_phi

                    ! get field following posision at 1/2 angle
                    call field_following_pos2(x,phi,phi_mid,xff)

                    call search_tr2(grid,xff,itr,p)

                    !remove particle or sheath calculation
                    if(itr<0) then
                       if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                          call remove_particle(sp,i,-1,ith)
                       else

                         call sheath_calculation(grid,psn,sp,i,sp%type,itr,p,ith)
                       endif
                    endif

                    sp%tr_save(i)=itr
                    sp%p_save(:,i)=p

                    select case(epc)
                    case(1)
                       sp%phase0(:,i)=sp%ptl(i)%ph
                    end select

                    call pushe_single(grid,psn,sp,i,sp%phase0(:,i),new_phase,dt_now,ith,diag_on)

                    ! check r-z boundary validity and update psi variables
                    if((new_phase(1)<eq_min_r) .or. &
                       (new_phase(1)>eq_max_r) .or. &
                       (new_phase(2)<eq_min_z) .or. &
                       (new_phase(2)>eq_max_z)) then
                       call remove_particle(sp,i,-1,ith)
! if (idebug >= 1) then
! print *, 'particle eliminated due to rz_outside :', &
! i, sml_mype, sp%type, sp%ptl(i)%gid, &
! new_phase(1),new_phase(2)
! endif
                    else
                       ! bounce 
                       if(epc==sml_nrk .and. sml_bounce/=0) then
                          old_phase(:)=sp%phase0(:,i)
                          call bounce(new_phase,old_phase,rtn)
                          if(rtn<0) then
                             call remove_particle(sp,i,-2,ith)
                          endif
                       endif

                       !******************************************************
                       ! time advance one step
                       !******************************************************
                       sp%ptl(i)%ph= new_phase(:)

                    endif

                    if(sp%ptl(i)%ph(3)>= sml_2pi .or. sp%ptl(i)%ph(3)< 0D0 ) then
                       sp%ptl(i)%ph(3)=modulo(sp%ptl(i)%ph(3),sml_2pi)
                    endif

                 endif
              enddo !i - particle
              call t_stopf("PUSHE_CPU_LOOP")

           enddo ! ith

        endif
        call t_stopf("cpu_processing")

        call t_startf("cuda_thrd_sync2")
        istatus = cudaStreamSynchronize(streamid)
        call t_stopf("cuda_thrd_sync2")

        call t_stopf("electron_loop")
     enddo
  enddo

  ! ------------------------------------------
  ! copy data from GPU device back to CPU host
  ! ------------------------------------------
  if (idebug >= 1) then
! if (sml_mype == 0) print*,icycle, phase0_gpu(1,1),phase0_gpu(2,1)
  endif


 if(gpu_ibegin <= gpu_iend) then
  call t_startf("push_upd_host_gpu")
  call push_update_host_gpu( sp, &
                             psn, diag_on_input, &
                             gpu_ibegin, gpu_iend)
  call t_stopf("push_upd_host_gpu")
  diag_on_input = diag_on
 endif
endif

call t_startf("PUSHE_SEARCH_INDEX")
call chargee_search_index(grid,psn,sp)
call t_stopf("PUSHE_SEARCH_INDEX")

  if (idebug >= 1) then
    if (sml_mype == 0) then
    write(*,9090) sml_gpu_ratio, sum(neu_weight_sum_lost)
 9090 format(' sml_gpu_ratio, sum(neu_weight_sum_lost) ',2(1x,1pe14.5))
    endif
  endif

  sml_epc=1
end subroutine pushe_gpu


attributes(global) &
subroutine pushe_kernel_gpu(istep,ipc,phase0, &
                   diag_on,dt_now,gpu_ibegin,gpu_iend)
use cudafor
use dimensions_mod_gpu !, only : nthreads_dim
use sml_module_gpu, only : sml_deltaf, sml_nrk, sml_bounce, sml_2pi
use ptl_module_gpu, only : ptl_nphase, ptl_gid_gpu, ptl_ph_gpu
use grid_class_gpu, only : grid_delta_phi
use eq_module_gpu, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z
use precision_mod_gpu

implicit none

real(kind=work_p) :: phase0(gpu_ibegin:gpu_iend,ptl_nphase)
real(kind=work_p) :: phase0_i(ptl_nphase)
integer, value :: istep, ipc !! RK4 index
logical, value :: diag_on
real(kind=work_p), value :: dt_now
integer, value :: gpu_ibegin,gpu_iend
real(kind=work_p), dimension(ptl_nphase) :: old_phase, new_phase
integer :: i, k, iv, ith, rtn 
integer, parameter :: idebug = 0 

real (kind=work_p) :: phi_tmp
real (kind=work_p) :: phi_mid, x(2), phi, mu, rho, xff(2)
real (kind=work_p) :: p(3)
integer :: itr, ip 

  ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) + &
      ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
      (blockDim%x * blockDim%y)


! if (ith.eq.1) print*, 'threadIdx%x=', ptl_ph_gpu(10,2)


  if ((ith < 1).or.(ith > nthreads_dim)) then
      return
  endif


   do i=gpu_ibegin+(ith-1),gpu_iend,nthreads_dim
        ! for alive particles only
        if(ptl_gid_gpu(i)>0) then

           ! get proper toroidal angle index and weight
           x=ptl_ph_gpu(i,1:2)
           phi=ptl_ph_gpu(i,3)
           phi_mid=(floor(phi/grid_delta_phi) + 0.5_work_p) * grid_delta_phi
! mu=sp%ptl(i)%ct(pim)
! rho=gyro_radius(x,mu) !gyro radius
! sp%rhoi(i)=rho

             ! get field following posision at 1/2 angle
 
! if(idebug.eq.1)print *,'before field_following_pos2_gpu'

           call field_following_pos2_gpu(x,phi,phi_mid,xff)

           call search_tr2_gpu(xff,itr,p)

           !remove particle or sheath calculation
           if(itr<0) then
              if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                 call remove_particle_gpu(i,-1)
              else
                 call sheath_calculation_gpu(i,ipc,0,itr,p,ith)
              endif
           endif

! tr_save_gpu(i)=itr
! p_save_gpu(:,i)=p
 
          
           !******************************************************
           ! actual particle push 
           !******************************************************

       select case(ipc)
       case(1)
             phase0(i,:)=ptl_ph_gpu(i,:)
       end select


           phase0_i = phase0(i,:)
           call push_single_gpu(i,phase0_i, &
                                new_phase,dt_now,ith,diag_on,itr,p)
! phase0(:,i) = phase0_i 
           
           ! check r-z boundary validity and update psi variables
           if(new_phase(1)<eq_min_r .or. new_phase(1)>eq_max_r .or. &
              new_phase(2)<eq_min_z .or. new_phase(2)>eq_max_z)then

              call remove_particle_gpu(i,-1)
                if (idebug >= 1) then
! print *, 'particle eliminated due to rz_outside :' , &
! i, sml_mype, sp%type, sp%ptl(i)%gid, &
! new_phase(1),new_phase(2)
                endif
           else 
              ! bounce 
              if(ipc==sml_nrk .and. sml_bounce/=0) then
                 old_phase(:)=phase0(i,:)
                 call bounce_gpu(new_phase,old_phase,rtn)
                 if(rtn<0) then
                     call remove_particle_gpu( i,-2)
                 endif
              endif
              
              !******************************************************
              ! time advance one step
              !******************************************************
              ptl_ph_gpu(i,:)= new_phase(:)
              
           endif

           if(ptl_ph_gpu(i,3)>= sml_2pi .or. ptl_ph_gpu(i,3)< 0D0 ) then
              ptl_ph_gpu(i,3)=modulo(ptl_ph_gpu(i,3),sml_2pi)
           endif

        endif
    enddo
  
end subroutine pushe_kernel_gpu
       subroutine gen_perm( ilo,ihi,jlo,jhi,guess_min,inv_guess_d, &
     & n, gid, x,y, iperm )

       use sml_module, only : sml_nthreads
       use omp_module, only : split_indices

       implicit none
       integer, intent(in) :: ilo,ihi, jlo,jhi
       real(8), intent(in) :: guess_min(2)
       real(8), intent(in) :: inv_guess_d(2)

       integer, intent(in) :: n
       integer*8, intent(in) :: gid(n)
       real(8), intent(in) :: x(n), y(n)
       integer, intent(inout) :: iperm(n)

       integer :: ith,isize, i,j,k, ifree,ipos
       real(8) :: guess_min_1, guess_min_2
       real(8) :: inv_guess_d_1, inv_guess_d_2

       integer :: i_beg(sml_nthreads),i_end(sml_nthreads)

       integer, allocatable, dimension(:,:,:) :: icount

       integer :: nerr
       integer, parameter :: idebug = 0
       integer :: icheck(n)


       guess_min_1 = guess_min(1)
       guess_min_2 = guess_min(2)
       inv_guess_d_1 = inv_guess_d(1)
       inv_guess_d_2 = inv_guess_d(2)

       call split_indices( n, sml_nthreads, i_beg, i_end)

       allocate( icount(ilo:ihi,jlo:jhi,sml_nthreads ) )

!$omp parallel do private(ith)
       do ith=1,sml_nthreads
          icount(:,:,ith) = 0
       enddo


! ----------------------------------------
! first pass to count number of particles
! in each cell
! ----------------------------------------

!$omp parallel do private(ith,k,i,j)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
         if (gid(k) > 0) then
           i = max(ilo,min(ihi, &
     & int((x(k) - guess_min_1) * inv_guess_d_1 + 1)))

           j = max(jlo,min(jhi, &
     & int((y(k) - guess_min_2) * inv_guess_d_2 + 1)))
         else
           i = ihi
           j = jhi
         endif

         icount(i,j,ith) = icount(i,j,ith) + 1
        enddo
       enddo

! ---------------------
! setup pointer 
! similar to prefix sum
! ---------------------

       ifree = 1
       do j=jlo,jhi
       do i=ilo,ihi
         do ith=1,sml_nthreads
           isize = icount(i,j,ith)
           icount(i,j,ith) = ifree
           ifree = ifree + isize
         enddo
       enddo
       enddo

! ----------------------
! second pass to setup perm
! ----------------------

!$omp parallel do private(ith,k,i,j,ipos)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
         if (gid(k) > 0) then
           i = max(ilo,min(ihi, &
     & int((x(k) - guess_min_1) * inv_guess_d_1 + 1 )))
           j = max(jlo,min(jhi, &
     & int((y(k) - guess_min_2) * inv_guess_d_2 + 1 )))
         else
! ----------------------------------------
! try to map invalid particles at the end
! ----------------------------------------
           i = ihi
           j = jhi
         endif

         ipos = icount(i,j,ith)
         iperm(ipos) = k
         icount(i,j,ith) = ipos + 1
        enddo
       enddo

       deallocate( icount )

! ------------
! double check
! ------------

       if (idebug >= 1) then

       icheck(1:n) = 0
       nerr = 0
       do i=1,n
         if ((1 <= iperm(i)).and.(iperm(i) <= n)) then
             icheck(iperm(i)) = icheck(iperm(i)) + 1
         else
             write(*,*) 'iperm out of range '
             write(*,*) 'i, iperm(i) ', i, iperm(i)
             nerr = nerr + 1
         endif
       enddo

       do i=1,n
         if (icheck(i) .ne. 1) then
           write(*,*) 'invalid iperm value '
           write(*,*) 'i, icheck(i) ', i,icheck(i)
           nerr = nerr + 1
         endif
       enddo

       if (nerr.ne.0) then
         stop '** error ** '
       endif

       endif


       return
       end subroutine gen_perm
attributes(device) &
subroutine sheath_calculation_gpu(iptl,epc,type_gpu,itrout,pout,ith)


  use grid_class_gpu
  use psn_class_gpu, only : &
     psn_nwall, &
     psn_wall_nodes, &
     psn_sheath_lost, &
     psn_pot0

  use sml_module_gpu
  use ptl_module_gpu
! use diag_module_gpu
  use precision_mod_gpu
  use neu_module_gpu, only : neu_weight_sum_lost_gpu

  implicit none

  integer, intent(in) :: type_gpu, iptl, epc, ith
  integer :: itrout
  integer :: i,l, node,itr
  real (kind=work_p) :: pout(3)
  real (kind=work_p) :: p(3),psave(3)
  integer, parameter :: rgn_wall=100
  real (kind=work_p), parameter :: minus_val=-1D50
  real (kind=work_p) :: rho,b,en_para, x(2),phi, phi_mid, xff(2)
  real (kind=work_p) :: time_now, xn(2), dist_sqr, dist_min
  integer :: node_min
  real (kind=work_p) :: b_interpol_gpu

  real (kind=work_p) :: new_phase(ptl_nphase), dummy
  integer :: widx, ip
  real (kind=work_p) :: w1_change, en_perp
  real (kind=work_p) :: tmp_ph(ptl_nphase), tmp_ct(ptl_nconst)
  ! find nearest wall point

  new_phase = phase0_gpu(iptl,:)
  x = new_phase(1:2)
  phi=new_phase(3)
  phi_mid=(floor(phi/grid_delta_phi) + 0.5_work_p) * grid_delta_phi

  ! get field following posision at 1/2 angle
  call field_following_pos2_gpu(x,phi,phi_mid,xff)

  ! find position of previous time step

     call search_tr2_gpu(xff,itr,p)
 
 
  psave=p

  ! if old position is also outside of grid --> remove particle and return
  if(itr<0) then
     call remove_particle_gpu(iptl,-1)
     itrout=itr
     pout=p
     return
  endif

  ! search three nodes of the triangle and check if it is wall nodes
  do i=1, 3
! l=maxloc(p,1)

     l = 3
     if ((p(2) >= p(1)).and.(p(2) >= p(3))) l = 2
     if ((p(1) >= p(2)).and.(p(1) >= p(3))) l = 1

     node = grid_nd(l,itr)

     if(grid_rgn(node)==rgn_wall) then
        exit
     else
        p(l)=minus_val
     endif
  enddo

  !if the triangle does not have a wall node
  ! find nearest one using primitive scan
  if(grid_rgn(node)/=rgn_wall) then
     dist_min = 1D99
     do i=1, psn_nwall ! for all wall nodes
        ! check distance
        xn=grid_x(:,psn_wall_nodes(i))
        dist_sqr=(xn(1)-xff(1))**2 + (xn(2)-xff(2))**2
        ! check minimum
        if(dist_min > dist_sqr) then
           dist_min=dist_sqr
           node_min = i
        endif
     enddo
     node=psn_wall_nodes(node_min)
  endif


  !
  !check potential and energy

  rho=new_phase(4)
  b = b_interpol_gpu(x(1),x(2),0.0_work_p)

  en_para=ptl_c2_2m(type_gpu) * (rho*b)**2

  ! reflection
  new_phase(pirho)=-new_phase(pirho)

  widx=psn_node_to_wall(node)
  if(en_para < psn_sheath_pot(widx) *sml_ev2j .and. type_gpu==0) then
     ! do nothing -- just relection due to sheath potential 
  else
     w1_change = 1.0_work_p - new_phase(piw2) + new_phase(piw1)

     ! sum-up to wall node

! psn_sheath_lost(widx)=psn_sheath_lost(widx) + (1D0 + new_phase(piw2))*sp%charge*sp%ptl(iptl)%ct(piw0)
     if(sml_ipc==2 .and. epc==1) dummy = atomicAdd(psn_sheath_lost(widx), (w1_change)*ptl_charge(type_gpu)*ptl_ct_gpu(iptl,piw0))
     ! --- reflect with weight change
     new_phase(piw1)=-1.0_work_p+new_phase(piw2)
     ! w2 does not change

     ! for neutral
     if(sml_neutral) then
        if(sml_ipc==2) then
          if(epc==1) then
            ip = mod( iptl, size(neu_weight_sum_lost_gpu)) + 1 
            dummy = atomicAdd( neu_weight_sum_lost_gpu(ip), w1_change*ptl_ct_gpu(iptl,piw0))
! neu_weight_sum_lost(ith) = neu_weight_sum_lost(ith) + (w1_change)*sp%ptl(iptl)%ct(piw0)
          endif
        endif
     endif

     ! for heat load diagnosie
     if(sml_ipc==2 .and. diag_heat_on .and. epc==1) then
! if(ipc_gpu==2)then 
        en_perp = ptl_ct_gpu(iptl,pim)*b
        ! arguments:
        ! 1. weight_change, 2. potential, 3. en_para, 4. en_perp, 5. ct, 6. phase(old), 7. phase(new), 8. stype
        ! new_phase is old position
        tmp_ph = ptl_ph_gpu(iptl,:)
        tmp_ct = ptl_ct_gpu(iptl,:)
! print *, "call diag_heat_port_gpu"
        call diag_heat_port_gpu(w1_change,psn_sheath_pot(widx),en_para, en_perp, tmp_ct, new_phase, tmp_ph,grid_delta_phi, type_gpu,ith)
     endif

  endif

  ptl_ph_gpu(iptl,:)=new_phase

  if(ptl_gid_gpu(iptl)>0) then
     itrout=itr
     pout=psave
  endif

end subroutine sheath_calculation_gpu

  ! function evaluation
  attributes(device) &
    subroutine ignore_factor_gpu(x,y,factor)
    use precision_mod_gpu
    use sml_module_gpu, only: sml_ignore_drift_r0, sml_ignore_drift_z0, &
        sml_ignore_drift_slope1, sml_ignore_drift_slope2

    implicit none
    real(kind=work_p), intent(in) :: x, y
    real(kind=work_p), intent(out) :: factor 
    real(kind=work_p) :: y2
    real(kind=work_p) :: dx


    dx=x-sml_ignore_drift_r0

    if(dx<0_work_p) then
       y2= sml_ignore_drift_z0 + sml_ignore_drift_slope1*dx
    else
       y2= sml_ignore_drift_z0 + sml_ignore_drift_slope2*dx
    endif

    if(y>y2) then
       factor=1_work_p
    else
       factor=0_work_p
    endif

    end subroutine ignore_factor_gpu



end module push_mod_gpu



