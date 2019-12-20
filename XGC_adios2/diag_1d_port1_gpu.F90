! port1 called in push   
attributes(device) &
subroutine diag_1d_port1_gpu(ptli,derivs,sp_type,vd,ith,tb)

#ifdef USE_GPU_EMU
#define atomicAdd atomicAdd_d
#endif


  use ptl_module_gpu, only : ptl_gid_gpu
  use diag_module, only : diag_1d_npv1
  use sml_module_gpu, only : sml_n_vf_diag, sml_deltaf, sml_ev2j, sml_mype, sml_inpsi
  use eq_module_gpu, only :  eq_x_z, eq_x2_z, eq_axis_r, eq_axis_z, eq_tempi
  use diag_module_gpu !, only : diag_1d_npv1, diag_1d_f_pv1, diag_1d_df_pv1, diag_1d_pin, &
!                              diag_1d_dp_inv, diag_1d_npsi
  use precision_mod_gpu
  implicit none
  integer, intent(in) :: ptli
  real (kind=work_p), intent(in) :: derivs(ptl_nphase)
  integer, intent(in) :: sp_type
  real (kind=work_p) :: vd(sml_n_vf_diag) ! 1: R_major, 2: B_toroidal, 3: B_total  4: radial ExB 5: v_para
  type(tbuf), intent(in) :: tb

  integer, intent(in) :: ith
  !
  real (kind=work_p) :: psi,z, pn, wp, b, r, rho, mu, w, dw, vp, en, den, pe, we
  real (kind=work_p) ::  diag_1d_de
  integer  :: ip, j, ie
  real (kind=work_p) :: v(diag_1d_npv1)

  real (kind=work_p) :: dval0, dummy0,dval1,dummy1
  real (kind=work_p) :: dval2,dummy2,dval3,dummy3,dval4,dummy4,dval5,dummy5,dval6,dummy6,dval7,dummy7
  integer :: lb,ub,i,i0,ith0


  if(ptl_gid_gpu(ptli) > 0) then
     psi = vd(1)
     r=tb%ph(1)
     z=tb%ph(2)
     if(.not. is_rgn12_gpu(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) return ! EXIT for private region of lower divertor
     pn=(psi-diag_1d_pin)*diag_1d_dp_inv
     ip=floor(pn)+1
     if(ip <1 .or. diag_1d_npsi <= ip) return     ! EXIT for out of diag_1d_pin/pout range
     wp=1.0_work_p - pn + real(ip-1,kind=work_p) 
     
     ! local variables for readability 
     b=vd(2)
!     r=ptl_ph_gpu(ptli,1)
     rho=tb%ph(4)
     mu=tb%ct(1)
     w=tb%ct(2)
     dw=tb%ph(5)*w
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


#ifdef ORIGINAL
     diag_1d_f_pv1(:,ip  ,sp_type,ith)=v(:)*w*wp              +  diag_1d_f_pv1(:,ip  ,sp_type,ith)
     diag_1d_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1.0_work_p-wp) +  diag_1d_f_pv1(:,ip+1,sp_type,ith)
         
#else
     ith0 = 1 + mod(ith, size(diag_1d_f_pv1,4))
     lb = lbound(diag_1d_f_pv1,1)
     ub = ubound(diag_1d_f_pv1,1)

     do i=lb,ub
        i0 = (i-lb) + lbound(v,1)
        dval0 =  (w * wp)*v(i0)
        dummy0 = atomicAdd( diag_1d_f_pv1(i,ip,sp_type,ith0), dval0 )

        dval1 = (w * (1.0_work_p-wp)) * v(i0)
        dummy1 = atomicAdd( diag_1d_f_pv1(i,ip+1,sp_type,ith0), dval1)
      enddo

#endif

     if(ptl_deltaf_sp(sp_type)) then
#ifdef ORIGINAL

        diag_1d_df_pv1(:,ip  ,sp_type,ith)=v(:)*dw*wp         +  diag_1d_df_pv1(:,ip  ,sp_type,ith)
        diag_1d_df_pv1(:,ip+1,sp_type,ith)=v(:)*dw*(1.0_work_p-wp) +diag_1d_df_pv1(:,ip+1,sp_type,ith)
#else

         ith0 = 1 + mod(ith,size(diag_1d_df_pv1,4))
         lb = lbound(diag_1d_df_pv1,1)
         ub = ubound(diag_1d_df_pv1,1)

         do i=lb,ub
           i0 = (i-lb) + lbound(v,1)
           dval0 = (dw*wp)*v(i0)
           dummy0 = atomicAdd( diag_1d_df_pv1(i,ip  ,sp_type,ith0), dval0)

           dval1 = dw*(1.0_work_p-wp) * v(i0)
           dummy1 = atomicAdd(diag_1d_df_pv1(i,ip+1,sp_type,ith0), dval1)
         enddo


#endif

     endif

!     if(diag_tavg_on .and. diag_omid_on) then
!        if(  r-eq_axis_r > abs(z-eq_axis_z) ) then
!             diag_1d_omid_f_pv1(:,ip  ,sp_type,ith)=v(:)*w*wp       +diag_1d_omid_f_pv1(:,ip  ,sp_type,ith)
!             diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1.0_work_p-wp) +diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)
!        endif
!     endif
     if(diag_eflux_on)then
!     if(.true.)then
        diag_1d_emin = 0.0_work_p
        diag_1d_emax = 3.0_work_p*eq_ftn_gpu2(sml_inpsi,eq_axis_r,eq_axis_z,eq_tempi)*sml_ev2j
        diag_1d_de=(diag_1d_emax-diag_1d_emin)/diag_1d_ne

        en=v(7)+v(8)
        pe=(en-diag_1d_emin)/diag_1d_de
        ie=floor(pe)+1
        if(ie >=1 .and. diag_1d_ne > ie) then     ! EXIT for out of diag_1d_emin/emax range
          we=1.0_work_p - pe + real(ie-1,kind=work_p)
#ifdef ORIGINAL
          diag_1d_eflux_pv(1,ip  ,ie  ,sp_type,ith)=w*wp*we                          +diag_1d_eflux_pv(1,ip  ,ie  ,sp_type,ith)
          diag_1d_eflux_pv(1,ip+1,ie  ,sp_type,ith)=w*(1.0_work_p-wp)*we             +diag_1d_eflux_pv(1,ip+1,ie  ,sp_type,ith)
          diag_1d_eflux_pv(1,ip  ,ie+1,sp_type,ith)=w*wp*(1.0_work_p-we)             +diag_1d_eflux_pv(1,ip  ,ie+1,sp_type,ith)
          diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith)=w*(1.0_work_p-wp)*(1.0_work_p-we)+diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith)

          diag_1d_eflux_pv(2,ip  ,ie  ,sp_type,ith)=v(9)*w*wp*we                          +diag_1d_eflux_pv(2,ip  ,ie  ,sp_type,ith)
          diag_1d_eflux_pv(2,ip+1,ie  ,sp_type,ith)=v(9)*w*(1.0_work_p-wp)*we             +diag_1d_eflux_pv(2,ip+1,ie  ,sp_type,ith)
          diag_1d_eflux_pv(2,ip  ,ie+1,sp_type,ith)=v(9)*w*wp*(1.0_work_p-we)             +diag_1d_eflux_pv(2,ip  ,ie+1,sp_type,ith)
          diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith)=v(9)*w*(1.0_work_p-wp)*(1.0_work_p-we)+diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith)
#else
          ith0 = 1 + mod(ith, size(diag_1d_eflux_pv,5))

          dval0 =  w * wp * we
          dummy0 = atomicAdd( diag_1d_eflux_pv(1,ip,ie,sp_type,ith0), dval0)

          dval1 = w * (1.0_work_p-wp) * we
          dummy1 = atomicAdd( diag_1d_eflux_pv(1,ip+1,ie,sp_type,ith0), dval1)

          dval2 =  w * wp * (1.0_work_p-we)
          dummy2 = atomicAdd( diag_1d_eflux_pv(1,ip,ie+1,sp_type,ith0), dval2 )

          dval3 = w * (1.0_work_p-wp) * (1.0_work_p-we)
          dummy3 = atomicAdd( diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith0), dval3)

          dval4 =  w * wp * we * v(9)
          dummy4 = atomicAdd( diag_1d_eflux_pv(2,ip,ie,sp_type,ith0), dval4 )

          dval5 = w * (1.0_work_p-wp) * we * v(9)
          dummy5 = atomicAdd( diag_1d_eflux_pv(2,ip+1,ie,sp_type,ith0), dval5)

          dval6 =  w * wp * (1.0_work_p-we) * v(9)
          dummy6 = atomicAdd( diag_1d_eflux_pv(2,ip,ie+1,sp_type,ith0), dval6 )

          dval7 = w * (1.0_work_p-wp) * (1.0_work_p-we) * v(9)
          dummy7 = atomicAdd( diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith0), dval7)

#endif
          if(ptl_deltaf_sp(sp_type)) then
#ifdef ORIGINAL
            diag_2d_dflux_pv(1,ip  ,ie  ,sp_type,ith)=dw*wp*we                          +diag_2d_dflux_pv(1,ip  ,ie  ,sp_type,ith)
            diag_2d_dflux_pv(1,ip+1,ie  ,sp_type,ith)=dw*(1.0_work_p-wp)*we             +diag_2d_dflux_pv(1,ip+1,ie  ,sp_type,ith)
            diag_2d_dflux_pv(1,ip  ,ie+1,sp_type,ith)=dw*wp*(1.0_work_p-we)             +diag_2d_dflux_pv(1,ip  ,ie+1,sp_type,ith)
            diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith)=dw*(1.0_work_p-wp)*(1.0_work_p-we)+diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith)

            diag_2d_dflux_pv(2,ip  ,ie  ,sp_type,ith)=v(9)*dw*wp*we                          +diag_2d_dflux_pv(2,ip  ,ie  ,sp_type,ith)
            diag_2d_dflux_pv(2,ip+1,ie  ,sp_type,ith)=v(9)*dw*(1.0_work_p-wp)*we             +diag_2d_dflux_pv(2,ip+1,ie  ,sp_type,ith)
            diag_2d_dflux_pv(2,ip  ,ie+1,sp_type,ith)=v(9)*dw*wp*(1.0_work_p-we)             +diag_2d_dflux_pv(2,ip  ,ie+1,sp_type,ith)
            diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith)=v(9)*dw*(1.0_work_p-wp)*(1.0_work_p-we)+diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith)
#else
            ith0 = 1 + mod(ith, size(diag_2d_dflux_pv,5))
         
            dval0 =  dw * wp * we
            dummy0 = atomicAdd( diag_2d_dflux_pv(1,ip,ie,sp_type,ith0), dval0 )

            dval1 = dw * (1.0_work_p-wp) * we
            dummy1 = atomicAdd( diag_2d_dflux_pv(1,ip+1,ie,sp_type,ith0), dval1)

            dval2 =  dw * wp * (1.0_work_p-we)
            dummy2 = atomicAdd( diag_2d_dflux_pv(1,ip,ie+1,sp_type,ith0), dval2 )

            dval3 = dw * (1.0_work_p-wp) * (1.0_work_p-we) 
            dummy3 = atomicAdd( diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith0), dval3)

            dval4 =  dw * wp * we * v(9)
            dummy4 = atomicAdd( diag_2d_dflux_pv(2,ip,ie,sp_type,ith0), dval4 )

            dval5 = dw * (1.0_work_p-wp) * we * v(9)
            dummy5 = atomicAdd( diag_2d_dflux_pv(2,ip+1,ie,sp_type,ith0), dval5)

            dval6 =  dw * wp * (1.0_work_p-we) * v(9)
            dummy6 = atomicAdd( diag_2d_dflux_pv(2,ip,ie+1,sp_type,ith0), dval6 )

            dval7 = dw * (1.0_work_p-wp) * (1.0_work_p-we) * v(9)
            dummy7 = atomicAdd( diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith0), dval7)

#endif

          endif  ! ptl_deltaf_sp
        endif  ! ie
     endif   ! diag_eflux_on
  endif
  
end subroutine diag_1d_port1_gpu
