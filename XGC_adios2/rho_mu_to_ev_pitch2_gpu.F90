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
!    if (idebug >= 1) print *, 'minus mu found :',rho,mu,b
    mu=0.0_work_p
  endif
  
  enj=(mu*b+ptl_c2_2m(sp_type)*(rho*b)**2)
  ev=enj*sml_j2ev
  v_pal=ptl_c_m(sp_type)*rho*b
  v_pep=SQRT(2.0_work_p*mu*b/ptl_mass(sp_type))
  pitch=v_pal/SQRT(v_pal**2+v_pep**2)
  
end subroutine rho_mu_to_ev_pitch2_gpu
