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
!    if(i(1)<=0 .or. i(1)>grid%guess_n(1)) then
!       print *, 'Invaild number for guess table- R',i(1),x(1),grid%guess_min(1),grid%guess_max(1)
!    endif
!    if(i(2)<=0 .or. i(2)>grid%guess_n(2)) then
!       print *, 'Invaild number for guess table- Z',i(2),x(2),grid%guess_min(2),grid%guess_max(2)
!    endif

    i=min(max(i,1),grid_guess_n)
    init=grid_guess_table(i(1),i(2))

  end subroutine guess_gpu

