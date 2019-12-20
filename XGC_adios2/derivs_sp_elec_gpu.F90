attributes(device) &
subroutine derivs_sp_elec_gpu(fld,t, ptli,  &
               yprime, diag_on, vf_diag, tb)
  use fld_module, only : fld_type
  use sml_module_gpu, only : &
    sml_n_vf_diag, sml_inpsi, sml_outpsi, sml_deltaf, &
    sml_e_charge, sml_deltaf_f0_mode, sml_f0_1_lt_e, &
    sml_f0_1_lt, sml_f0_1_ln,  sml_dwdt_fix_bg,sml_dwdt_exb_only, sml_ignore_drift_near_wall
  use ptl_module_gpu, only : pim, pirho, piw1,piw2,&
      ptl_nphase, ptl_charge, ptl_mass, &
      ptl_nconst, ptl_deltaf_sp
  use eq_module_gpu, only : eq_den, eq_tempe, eq_tempi
  use precision_mod_gpu
  implicit none
  type(fld_type), intent(in) :: fld   ! field variable
  real (kind=work_p),       intent(in) :: t     ! time

  integer, intent(in) :: ptli  ! particle info
!  integer,  intent(in)  :: sp_type      ! particle species type (ion/elec)
  real (kind=work_p), intent(out) :: yprime(ptl_nphase)  !
  logical,  intent(in)  :: diag_on
  real (kind=work_p), intent(inout) :: vf_diag(sml_n_vf_diag)   ! variables for diagnosis 
  type(tbuf), intent(in) :: tb
  !
  real (kind=work_p) :: mass, charge,c_m  ! charge and mass 
  real (kind=work_p) :: r, z, phi, rho, mu,inv_r
  real (kind=work_p) :: B, B2, over_B, over_B2
  real (kind=work_p) :: D, nb_curl_nb
  real (kind=work_p) :: dbdr, dbdz, dbdphi
  real (kind=work_p) :: cmrho2, cmrho, murho2b, murho2b_c, vp
  real (kind=work_p) :: fr, fp, fz
  real (kind=work_p) :: fr_exb, fp_exb, fz_exb, yp_exb(3)

  ! for weight calculation
  real (kind=work_p) :: energy, vmag, pitch, dvpdt, denergy_dt, dpitch_dt
  real (kind=work_p) :: psi, den, dden, temp, dtemp, dfdp, tmp, envelop, df0(5)
  real (kind=work_p) :: one_m_w
  real (kind=work_p) :: total_ddpotdt    ! for electron -- from adiabatic reponse
  !
  real (kind=work_p) :: bp, f0_1_lt
  real (kind=work_p) :: i_factor
  ! some global parameter to be
  real (kind=work_p) :: sml_wdot_energy_max, sml_f0_psi_c, sml_f0_1_psi_w
  
  ! these are constant for deltaf method -- maybe should be moved to sml_module

  if(sml_deltaf_f0_mode==-1) then
     ! these are constant for deltaf method -- maybe should be moved to sml_module
     sml_wdot_energy_max=10D0
     sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
     sml_f0_1_psi_w=1D0/( 0.4*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )
  endif

  ! prepare constant
  mass=ptl_mass(0) !
  charge=ptl_charge(0) !-1.D0/ptl_charge
  c_m=charge/mass
    
  r= tb%ph(1)
  z= tb%ph(2)
  phi= tb%ph(3)
  rho= tb%ph(4)
  mu= tb%ct(pim)
  inv_r=1D0/r

  B = sqrt( fld%br**2 + fld%bz**2 + fld%bphi **2 )
  b2=b*b
  over_B=1/B
  over_B2=over_B*over_B
  
  ! normalized b dot curl of normalized b
  nb_curl_nb= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/r  + fld%dbpdr)*fld%bz )
  D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
      
  dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  dbdphi=0D0  ! no B perturbation


  ! bug fix -- rho**2 B grad B term added 2002/1/15
  !! optimization variables
  
  cmrho2=c_m*rho**2
  cmrho =c_m*rho
  murho2b=(mu+charge*cmrho2 *B)
  murho2b_c=murho2b/charge
  vp=cmrho*B

  ! F is grad H/q
  fr = fld%Er - (murho2B_c) * dbdr ! fr is force over q . -gradient of  Hamiltonian
  fp = fld%Ephi - (murho2B_c) * dbdphi/r  ! modified by shlee 5/30/2001 -- /r added
  fz = fld%Ez - (murho2B_c) * dbdz

  ! ignore all drifts when it close to divertor
  ! ion & electron --> pushe for electron
  ! now drift and parallel E-field are ignored.
  if(sml_ignore_drift_near_wall ) then
     call ignore_factor_gpu(r,z,i_factor)
     !ep_r = (fld%Er  *fld%br  )*fld%br  /b2
     !ep_p = (fld%Ephi*fld%bphi)*fld%bphi/b2 
     !ep_z = (fld%Ez  *fld%bz  )*fld%bz  /b2 

     fr = fr*i_factor !+ ep_r*(1D0-i_factor)
     fp = fp*i_factor !+ ep_p*(1D0-i_factor)
     fz = fz*i_factor !+ ep_z*(1D0-i_factor)     
  endif
  
  yprime(1)= D*( (fld%bz*Fp - fld%Bphi * Fz) * over_B2         &
       +  cmrho * fld%br                       &
       +  cmrho2 * (fld%dbzdp*inv_r - fld%dbpdz ) )
  yprime(2)= D*( (fld%bphi * fr - fld%br * fp ) * over_B2      &
       +  cmrho * fld%bz                       &
       +  cmrho2 * (fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r) )
  yprime(3)= D*( (fld%br * fz - fld%bz * fr) * over_B2         &
       +  cmrho * fld%bphi                     &
       +  cmrho2 * ( fld%dbrdz - fld%dbzdr) ) * inv_r


  fr_exb = fld%Er  ! fr is force over q . -gradient of  Hamiltonian
  fp_exb = fld%Ephi! modified by shlee 5/30/2001 -- /r added
  fz_exb = fld%Ez 
  
  yp_exb(1)= D*(fld%bz   *fp_exb - fld%Bphi * fz_exb) * over_B2  
  yp_exb(2)= D*(fld%bphi *fr_exb - fld%br   * fp_exb) * over_B2      
  yp_exb(3)= D*(fld%br   *fz_exb - fld%bz   * fr_exb) * over_B2 *inv_r

  yprime(pirho)=D*over_B2 *( &
       fld%br*fr + fld%bz*fz + fld%bphi*fp &
       + rho*( fr*(fld%dbzdp*inv_r-fld%dbpdz) + fz*(fld%bphi*inv_r+fld%dbpdr-fld%dbrdp*inv_r) + fp*(fld%dbrdz-fld%dbzdr)) &
       )
  
  if( ptl_deltaf_sp(0) .and. sml_extra_dwdt ) then
     energy=(0.5D0*mass*vp*vp + mu*b)
     vmag=sqrt(2D0*energy/mass)
     pitch = vp/vmag
     dvpdt=c_m*yprime(pirho)*B + (dbdr*yprime(1)+dbdz*yprime(2))*vp/B  !
     denergy_dt=mu*(dbdr*yprime(1)+dbdz*yprime(2)+dbdphi*yprime(3))+(mass*vp*dvpdt)
     dpitch_dt=dvpdt/vmag - 0.5D0*pitch*denergy_dt/energy


          ! dlog(f0) ------ need to be separate routine later     
     ! temp routine -- finish it later
     ! no space gradient
     psi=fld%psi

     !den=eq_ftn(psi,z,eq_den)
     den = eq_ftn_gpu2(psi,r,z,eq_den)


     !dden=eq_dftn(psi,z,eq_den)
     dden = eq_dftn_gpu2(psi,r,z,eq_den)

!     if(sp_type==0) then !electron
        !temp=eq_ftn(psi,z,eq_tempe)*sml_e_charge
        temp=eq_ftn_gpu2(psi,r,z,eq_tempe)*sml_e_charge
        ! dtemp=eq_dftn(psi,z,eq_tempe)*sml_e_charge
        dtemp=eq_dftn_gpu2(psi,r,z,eq_tempe)*sml_e_charge
!     else
        ! temp=eq_ftn(psi,z, eq_tempi)*sml_e_charge
!        temp=eq_ftn_gpu2(psi,r,z, eq_tempi)*sml_e_charge
        ! dtemp=eq_dftn(psi,z,eq_tempi)*sml_e_charge
!        dtemp=eq_dftn_gpu2(psi,r,z,eq_tempi)*sml_e_charge
!     endif
     
     if(sml_deltaf_f0_mode==-2) then  !consistent grad f with real profile
        dfdp=dden/den + dtemp/temp*(energy/temp - 1.5D0)
     elseif(sml_deltaf_f0_mode==-1) then
        tmp=1D0/sqrt(fld%dpsidr**2+fld%dpsidz**2)  ! drdpsi 
        envelop= exp( - ((sqrt(psi) - sml_f0_psi_c)*sml_f0_1_psi_w )**8 )
        !envelop=1D0
!        if(sp_type==0) then
           f0_1_Lt=sml_f0_1_Lt_e
!        else
!           f0_1_Lt=sml_f0_1_Lt
!        endif
        dfdp=-envelop*(sml_f0_1_Ln + f0_1_Lt*(energy/temp - 1.5D0))*tmp
     endif
     df0(1)=dfdp *fld%dpsidr !/(2D0*q+1D-99)
     df0(2)=dfdp *fld%dpsidz !/(2D0*q+1D-99)
     df0(3)=0D0  * r     ! r : yprime(3) is phi dot --> dimension matching
#ifdef DELTAF_MODE2
     df0(4)=-1D0/temp ! energy deriv
#else
     df0(4)=0D0 ! zero for extra dwdt
#endif
     df0(5)=0D0              ! pitch deriv

     if(sml_dwdt_fix_bg)then 
        one_m_w=1D0
     else
        one_m_w=1D0 - tb%ph(piw2)
     end if

     if( .not. sml_dwdt_exb_only ) then
        yprime(piw1)= -one_m_w* (&
             df0(1)*yprime(1) + df0(2)*yprime(2) + df0(3)*yprime(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )
#ifndef DELTAF_MODE2
     else if(sml_deltaf_f0_mode == -1) then
#else
     else if(sml_deltaf_f0_mode == -1 .or. sml_deltaf_f0_mode == -2) then  ! else-if below will be ignored.
#endif
        yprime(piw1)= -one_m_w* (&
             df0(1)*yp_exb(1) + df0(2)*yp_exb(2) + df0(3)*yp_exb(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )
     else if(sml_deltaf_f0_mode == -2) then  ! remove weight evolution from v_grad B  - total will be updated later
        yprime(piw1)= -one_m_w* (&
             df0(1)*(yp_exb(1)-yprime(1)) + df0(2)*(yp_exb(2)-yprime(2)) )

     endif
     ! electron -- minus adibatic response
!     if(sp_type==0) then
        total_ddpotdt=fld%ddpotdt - yprime(1)*(fld%Er-fld%Er00) - yprime(2)*(fld%Ez-fld%Ez00)  -r*yprime(3)*fld%Ephi

        yprime(piw1) = yprime(piw1) - one_m_w*total_ddpotdt/temp*sml_e_charge
!     endif
     yprime(piw2)=yprime(piw1)
  else
     yprime(piw1:piw2)=0D0
  endif
  if(diag_on) then
     vf_diag(1)=fld%psi
     vf_diag(2)=B
     vf_diag(3)=fld%Bphi
     vf_diag(4)=yp_exb(1)*fld%dpsidr + yp_exb(2)*fld%dpsidz   ! V_ExB dot grad psi
     vf_diag(5)=yprime(1)*fld%dpsidr + yprime(2)*fld%dpsidz   ! V_d dot grad psi
     vf_diag(6)=fld%dpsidr*fld%dpsidr + fld%dpsidz*fld%dpsidz  ! |grad psi|^2
     bp=sqrt( fld%br*fld%br + fld%bz*fld%bz ) 
     vf_diag(7)=(yprime(1)*fld%br + yprime(2)*fld%bz)/bp
     vf_diag(8)=(yp_exb(1)*fld%br + yp_exb(2)*fld%bz)/bp
  endif

end subroutine derivs_sp_elec_gpu
