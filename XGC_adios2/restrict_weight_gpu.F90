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

