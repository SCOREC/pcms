attributes(device) &
subroutine derivs_single_with_e_gpu(ptli,dy,i,time,E_mag,tb)
  use sml_module_gpu, only : sml_mype, sml_n_vf_diag, sml_2pi_wedge_n
  use fld_module
  use ptl_module_gpu, only : type_gpu
  use precision_mod_gpu
  use grid_class_gpu, only : search_tr2_gpu, search_tr_check_guess_gpu
#ifdef RK4_ACCURATE_EFIELD
  use psn_class_gpu, only : efield_gk_elec_gpu
#endif
  implicit none

  integer,intent(in) :: ptli
  real (kind=work_p), intent(out) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (kind=work_p), intent(in) :: time 
  real (kind=work_p), intent(in) :: E_mag(3)
  type(tbuf), intent(inout) :: tb
!  integer, intent(in) :: ith
  !
  type(fld_type) :: fld
  logical :: rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)
  logical, parameter :: diag_on=.false.
  !
  real (kind=work_p) :: dpsi(2), E(3), bp, dtheta_norm(2), B
  real (kind=work_p) :: p(3), x(2), phi, phi_mid, xff(2)
  integer :: itr

  ! Save space information
  fld%r=tb%ph(1)
  fld%z=tb%ph(2)
  fld%phi=tb%ph(3)

  ! obtain B-field information 
  call field_gpu(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field from previous E-field

#ifdef RK4_ACCURATE_EFIELD
     x(1)=fld%r
     x(2)=fld%z
     if(fld%phi>= sml_2pi_wedge_n .or. fld%phi< 0D0 ) then
        fld%phi=modulo(fld%phi,sml_2pi_wedge_n)
     endif
     phi =fld%phi
     phi_mid=(floor(phi/grid_delta_phi) + 0.5_work_p) * grid_delta_phi

     call field_following_pos2_gpu(x,phi,phi_mid,xff)

#ifdef USE_TR_CHECK
     call search_tr_check_guess_gpu(xff,tb%tr,itr,p)
     if (itr < 0) then
       call search_tr2_gpu(xff,itr,p)
       tb%tr = itr
     endif
#else
     call search_tr2_gpu(xff,itr,p)
#endif

     if (itr .gt. 0) then
       call efield_gk_elec_gpu(i,fld,itr,p)
     else
#endif
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
#ifdef RK4_ACCURATE_EFIELD
     endif
#endif
!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp_elec_gpu(fld,time,ptli,dy,diag_on,vf_diag,tb)
     
  else
     call remove_particle_gpu(i,-1,tb)
!     call remove_particle_gpu( sp%ptl(i)%gid, &
!                           sp%ptl(i)%ph, &
!                           sp%ptl(i)%ct, &
!                           sp%type, &
!                           i, -1 )

!     print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  end if
end subroutine derivs_single_with_e_gpu
