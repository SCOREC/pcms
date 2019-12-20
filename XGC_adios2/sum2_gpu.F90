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

        ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) +  &
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


!        print*,'m = ', m, 'tgrid%x ',tgrid%x, 'tgrid%y ',tgrid%y


        streamid = 0
        if (present(streamid_in))  streamid = streamid_in
        call sum2_kernel_gpu<<<tgrid,tblock,0,streamid>>>(m,n,           &
     &                                   A_gpu,Asum_gpu)
        Asum = Asum_gpu
        return
        end subroutine sum2_gpu
