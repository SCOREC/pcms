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
