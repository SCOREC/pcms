attributes(device)  &
real (kind=8)  function I_interpol_gpu(in_psi, ideriv,region)

    use one_d_cub_mod_gpu
    use itp_module_gpu, only : itp_min_psi,itp_max_psi
    use eq_module_gpu, only : eq_x_psi
    use sml_module_gpu, only : sml_bt_sign
    use precision_mod_gpu
    implicit none
    integer , intent(in) :: ideriv,region
    real (kind=work_p) , intent(in) :: in_psi
    real (kind=work_p)  :: psi
    integer :: ier
    real (kind=work_p) :: r8value


!    sign=-1D0 !-1D0 : cocurrent, 1D0 :counter current 
!    sml_bt_sign can be changed in setup.f90 2002/02/08
    
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
