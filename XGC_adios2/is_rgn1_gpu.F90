  ! function evaluation
  attributes(device) &
  logical  function is_rgn1_gpu(r,z,psi)
    use precision_mod_gpu
    use eq_module_gpu, only : eq_x_psi, eq_x_r, eq_x_z, eq_x_slope, &
                              eq_x2_r, eq_x2_z, eq_x2_slope
    implicit none
    real (kind=work_p), intent(in) :: r,z,psi
    real (kind=work_p) , parameter :: epsil_psi =  1D-5

    if(psi <= eq_x_psi-epsil_psi .and. &
         -(r-eq_x_r)*eq_x_slope + (z-eq_x_z) > 0.0_work_p .and. -(r-eq_x2_r)*eq_x2_slope &
         + (z-eq_x2_z) < 0.0_work_p   ) then
       is_rgn1_gpu=.true.
    else
       is_rgn1_gpu=.false.
    endif
  end function is_rgn1_gpu

