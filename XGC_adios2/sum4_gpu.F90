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
        i2 = (blockIdx%x  -1) + l2
        i3 = (blockIdx%y  -1) + l3

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
        subroutine sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4,A_gpu,              &
     &              Asum, streamid_in )
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

!        print*,'before sum4_kernel_gpu'
        streamid = 0
        if (present(streamid_in)) streamid = streamid_in
        call sum4_kernel_gpu<<<tgrid,tblock,0,streamid>>>(l1,u1,l2,u2,l3,u3,l4,u4, &
                                       A_gpu,Asum_gpu)
!        print*,'after sum4_kernel_gpu '
        Asum = Asum_gpu

!        print*,'after Asum  = Asum_gpu '

        deallocate( Asum_gpu )
!        print*,'after deallocate(Asum_gpu) '
        return
        end subroutine sum4_gpu
