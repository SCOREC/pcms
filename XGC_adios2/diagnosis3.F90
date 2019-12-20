subroutine diag_favg_port(ptli,derivs,sp_type,vd,ith)
  use ptl_module
  use diag_module
  use sml_module, only : sml_n_vf_diag, sml_deltaf, sml_ev2j, sml_mype, sml_inpsi
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z, eq_axis_r, eq_axis_z,eq_ftn, eq_tempi
  implicit none
  type(ptl_type), intent(in) :: ptli
  real (8), intent(in) :: derivs(ptl_nphase)
  integer, intent(in) :: sp_type
  real (8) :: vd(sml_n_vf_diag) ! 1: R_major, 2: B_toroidal, 3: B_total  4: radial ExB 5: v_para
  integer, intent(in) :: ith
  !
  real (8) :: psi,z, pn, wp, b, r, rho, mu, w, dw, vp, en, den, pe, we
  real (8) ::  diag_1d_de
  integer  :: ip, j, ie
  real (8) :: v(diag_favg_npv)

  if(ptli%gid > 0) then   ! only when gid is positive
     psi = vd(1)
     r=ptli%ph(1)
     z=ptli%ph(2)
     if(.not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) return ! EXIT for private region of lower divertor


     call search_psi_diag_favg(psi, ip, wp)

     ! local variables for readability
     b=vd(2)

     rho=ptli%ph(4)
     mu=ptli%ct(1)
     w=ptli%ct(2)
     dw=ptli%ph(5)*w
     vp=ptl_c_m(sp_type)*rho*B

     ! obtain variables
     v(1) = 1D0                     ! full weight
     v(2) = derivs(3)*r             ! g.c. toroidal velocity --  dzeta/dt * R_major
     v(3) = vd(7)                   ! g.c. poloidal velocity -- v dot Bp/Bp
     v(4) = vp                      ! v_|| parallel velocity
     v(5) = vp*r*vd(3)/vd(2)        ! v_zeta_|| / R --toroidal angular momentum of v_||   : v_|| * B_zeta /(R*B)
     v(6) = vd(5)                   ! radial drift  - psi_dot
     v(7) = ptl_c2_2m(sp_type)*(rho*B)**2 ! parallel mean energy
     v(8) = mu*B                    ! perp temperature
     v(9) = (v(7)+v(8))*v(6)        ! Energy radial flux
     v(10)= v(4)*vd(4)*vd(3)/vd(2)  ! V_exb * V_|| * B_phi/B
     v(11)= v(4)*vd(5)*vd(3)/vd(2)  ! V_r   * V_|| * B_phi/B
     v(12)= vd(4)                   ! radial drift by exb  - V_exb dot grad_psi
     v(13)= (v(7)+v(8))*vd(4)       ! heat flux by radial exb
     v(14)= vd(8)                   ! poloidal comp. of V_ExB
     v(15)= vd(6)                   ! grad_psi ^2

     diag_1d_f_pv1(:,ip  ,sp_type,ith)=v(:)*w*wp       +diag_1d_f_pv1(:,ip  ,sp_type,ith)
     diag_1d_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1D0-wp) +diag_1d_f_pv1(:,ip+1,sp_type,ith)


     if(ptl_deltaf_sp(sp_type)) then
        diag_1d_df_pv1(:,ip  ,sp_type,ith)=v(:)*dw*wp       +diag_1d_df_pv1(:,ip  ,sp_type,ith)
        diag_1d_df_pv1(:,ip+1,sp_type,ith)=v(:)*dw*(1D0-wp) +diag_1d_df_pv1(:,ip+1,sp_type,ith)
     endif

     


end subroutine diag_favg_port


subroutine diag_favg_output





end subroutine diag_favg_output
