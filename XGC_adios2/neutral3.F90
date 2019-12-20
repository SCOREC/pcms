subroutine neutral3(istep,grid,psn,sp)
  use grid_class
  use psn_class
  use ptl_module
  use neu_module
  use sml_module
  use eq_module
  use col_module, only: col_mode,col_varying_bg,col_vb_period
  use perf_monitor
  use f0_module
  use random_xgc
  implicit none
  type(species_type) :: sp
  type(grid_type) :: grid
  type(psn_type) :: psn  
  integer ,intent(in) :: istep
  
  integer :: inode, imu, ivp, ith, isp, idebug=1, inode2, j, k, itr, nsub=1, nd(3), i, itmp, &
             imup, ivpp, imu_h, ivp_h
  real (8) :: r, z, psi, Tn, n, b, mu_n, vp_n, P_ion, P_cx, en_neut, ene_neut, Te, enev 
  real (8) :: rate_neu_col, rate_neu_ion, sum_neui_ionized, sum_neue_ionized, factor, &
              sum_neu_weight_lost, neui_ionized, neue_ionized, smu, E_ion, new_n0, &
              r2, z2, psi2, n_ave, area_sum, n_tot, Te2, c1, c2, c3,  area, mu_np, vp_np, &
              Tc, Th, smup, smuh, en_hot, mu_nh, vp_nh, sum_neue_lost, sum_neue_h, sum_neue_c, &
              neue_cold, neue_hot, neue_loss, w1, w2, eperp, fsum, fesum, f, b_b0, mu_vol, &
              sum_ion_cx, sum_loss_cx

  real (8),external :: b_interpol,neutral_den,psi_interpol,te_neutral_col_3,neutral_temp_3, &
                       neu_psi_interpol
  real (8), allocatable :: maxwell_neut(:,:,:,:)
  real (8), allocatable :: delta_fi_neut_ion(:,:,:)
!  real (8), allocatable :: delta_fe_neut_ion(:,:,:)
  real (8), allocatable :: delta_fi_neut_cx(:,:,:)
!  real (8), allocatable :: delta_fe_neut_ion_c1(:,:,:)
  real (8), allocatable :: delta_fe_neut_ion_c2(:,:,:)
  real (8), allocatable :: delta_fe_neut_ion_h(:,:,:)
  real (8), allocatable :: delta_fe_neut_ion_l(:,:,:)
  real (8), allocatable :: delta_fi_neut_cx_M(:,:,:)
  real (8), allocatable :: sum_delta_fi(:)
  real (8), allocatable :: sum_delta_fi_M(:)
  real (8), allocatable :: n_p(:), n_neut(:), n_p2(:), fion(:), eperp1(:),epara(:), Te_ev(:)
  real (kind=8) :: dum1(3), dum2(3), xc(2), dx1(2), dx2(2)
  real (8) :: tmp, f0_f_tmp
  logical :: rgn_check
  logical :: first=.true.
  save first
  logical, external :: is_nan 
!  if(istep < neu_start_time) return ! no neutral effect before a couple profile evaluations
  E_ion=30D0 ! ev electron ionization energy (This is an assumption based on Daren's DEGAS2)

  if(sml_electron_on) then
    isp=0
  else
    isp=1
  endif


  if((neu_col_mode==2) .and. first) then
     ! neutral diffusion calculation

     call check_point("before neutral_step")
     call t_startf("NEUTRAL_STEP")
     call neutral_step(ptl_mass(1))
     call t_stopf("NEUTRAL_STEP")
     call check_point("after neutral_step")

     first=.false.  ! since plasma profile is not adjusted, yet
  endif

  if(neu_enforce_no_neutral) return ! This is for radiation cooling which requires neutral density

  !! Set up the Maxwellian distribution function for neutrals
  allocate(maxwell_neut(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2, &
         ptl_isp:ptl_nsp), n_p(f0_inode1:f0_inode2), &
         n_neut(f0_inode1:f0_inode2), &
         Te_ev(f0_inode1:f0_inode2), &
         n_p2(f0_inode1:f0_inode2), fion(f0_inode1:f0_inode2), &
!         delta_fe_neut_ion_c1(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2), &
         delta_fe_neut_ion_c2(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2), &
         delta_fe_neut_ion_h(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2), &
         delta_fe_neut_ion_l(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2))

! Get the plasma (electron) density 
  call t_startf("NEUTRAL_ELEC_DENSITY")
  do inode = f0_inode1, f0_inode2
     r = grid%x(1,inode) ! get major radius
     z = grid%x(2,inode)
     psi = grid%psi(inode)
     n_p(inode)=eq_ftn(psi,r,z,eq_den)
     Te_ev(inode)=eq_ftn(psi,r,z,eq_tempi)
     if(sml_electron_on)then
!       n_p(inode)=0D0
!       do imu = f0_imu1, f0_imu2
!          do ivp = -f0_nvp, f0_nvp
!             n_p(inode)=n_p(inode)+f0_f(ivp,inode,imu,0)*f0_grid_vol(inode,0)/grid%node_vol_nearest(inode)
!          enddo
!       enddo
        Te_ev(inode)=eq_ftn(psi,r,z,eq_tempe)
     endif
  enddo

! estimate the varied background electron density and temperature
! prepare eperp and epara
  allocate(eperp1(f0_imu1:f0_imu2), epara(-f0_nvp:f0_nvp))

  do imu=f0_imu1,f0_imu2
     smu=imu*f0_dsmu
     mu_n=smu*smu !!rh mu_n is either normalized mag. moment or (v_perp/v_th)^2
     eperp1(imu)=0.5D0*mu_n
  enddo

  do ivp=-f0_nvp, f0_nvp
     epara(ivp) = 0.5D0*(ivp*f0_dvp)**2
  enddo
  call t_stopf("NEUTRAL_ELEC_DENSITY")

  !get density and temperature of plasma for each grid point.
  call t_startf("NEUTRAL_COL_SNAPSHOT")
  do inode=f0_inode1, f0_inode2

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle

     ! simulation boundary check
     psi = grid%psi(inode)
     if(psi < sml_inpsi .or. sml_outpsi < psi) cycle

     tmp=n_p(inode)  ! original density

     ! obtain n & T
     !vsum=0D0
     fsum=0D0
     fesum=0D0
     b_b0=f0_b_b0(inode)
     do imu=f0_imu1, f0_imu2

        ! normalized volume element
        if(imu==0 .or. imu==f0_nmu) then
           mu_vol=0.5D0
        else
           mu_vol=1D0
        endif

        eperp=eperp1(imu)
        do ivp=-f0_nvp,f0_nvp
           f=f0_f(ivp,inode,imu,isp)
           f=max(f,0D0)
           !vsum  = vsum  + mu_vol
           fsum  = fsum  + mu_vol*f
           fesum = fesum + mu_vol*f*(eperp+epara(ivp))
        enddo  ! ivp
     enddo ! imu
     if (neu_col_mode /= 0) then
        if(col_varying_bg==1 ) then
          n_p(inode)=fsum*f0_grid_vol_vonly(inode,isp)
          ! set minimum n_p to be 1% of original density
          if(n_p(inode) < 1D-2*tmp) then
             print *, 'too low plasma density. using 1% of equilbrium density', tmp, n_p(inode)
             n_p(inode) = 1D-2*tmp
          endif

          Te_ev(inode)=2D0/3D0*fesum/fsum*f0_T_ev(inode,isp)
          ! set minimum Te_ev to be 1% of original density
          if(Te_ev(inode) < 1D-2*f0_T_ev(inode,isp)) then
             print *, 'too low plasma temperature. using 1% of equilbrium density', Te_ev(inode),f0_T_ev(inode,isp)
             Te_ev(inode) = 1D-2*f0_T_ev(inode,isp)
          endif
     !debug
!          if(sml_mype==0) print *, 'n & T', n_p(inode), Te_ev(inode)
        endif
     endif
  enddo !inode
  call t_stopf("NEUTRAL_COL_SNAPSHOT")

  !call check_point( "after maxwell_neut allocation")

  call t_startf("NEUTRAL_MAXWELLIAN")
  en_neut=0D0
  maxwell_neut=0D0
  do inode = f0_inode1, f0_inode2

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle

     r = grid%x(1,inode) ! get major radius
     z = grid%x(2,inode)
     psi = grid%psi(inode)
! Get neutral temperature
     Tn=neutral_temp_3(r,z,psi) ! in J
     Tn=Tn*sml_j2ev  ! convert to ev
!     n_neut(inode)=neutral_den(r,z,psi) ! normalized unit
     B=B_interpol(r,z,0D0) ! should be B* to be exact
!if(idebug==1) print *, 'maxwell_neut', n,  r,z, psi, Tn, sml_mype
!    return 
     do imu = f0_imu1, f0_imu2
        do ivp = -f0_nvp, f0_nvp
!           do isp = ptl_isp, ptl_nsp
! normalized
              smu=(imu*f0_dsmu) ! sqrt mu (normalized) or v_perp/vth
              mu_n = smu*smu  !!rh mu_n is either normalized mag. moment or (v_perp/v_th)^2
              vp_n = ivp*f0_dvp
              en_neut = 0.5D0*(mu_n + vp_n**2)*f0_T_ev(inode,1)/(Tn+1D-10) !normalized en_neut
              maxwell_neut(ivp,inode,imu,:) = 1D0/Tn*exp(-en_neut)*sqrt(mu_n*f0_T_ev(inode,1)/Tn)
!           enddo ! isp
        enddo ! ivp  
     enddo  ! imu
  enddo  ! inode
  call t_stopf("NEUTRAL_MAXWELLIAN")

!if(idebug==1) print *, 'maxwell_neut', n, mu_n,vp_n,en_neut, r,z, psi, Tn, Ti_ev, Te_ev

  !! Ionization ----------------------------------------------------------------
!  if(mod(istep,sml_f_source_period)==0) then
  ! Set up ionization probabilities for each grid in v-space
  allocate(delta_fi_neut_ion(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2))
!  call check_point("after delta_fi_neut_ion allocation")
!  if(sml_electron_on)then
!     allocate(delta_fe_neut_ion(-f0_nvp:f0_nvp, &
!         f0_inode1:f0_inode2,&
!         f0_imu1:f0_imu2))
!  endif

  sum_neui_ionized = 0D0
  sum_neue_ionized = 0D0
  neui_ionized = 0D0
  neue_ionized = 0D0
  sum_neu_weight_lost=0D0
  delta_fi_neut_ion=0D0
  fion=0D0
  n_neut=0D0

  call t_startf("NEUTRAL_IONIZE_PROB")
  do inode = f0_inode1, f0_inode2

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle

     r = grid%x(1,inode) ! get major radius
     z = grid%x(2,inode)
     psi = grid%psi(inode)
!     Te=Te_neutral_col_3(psi,r,z) ! eV
     Te = Te_ev(inode)
     n=neutral_den(r,z,psi) ! normalized unit
     area_sum=0D0
     n_ave=n
     n_tot = 1D0
!print *, 'before avergage loop'
     do i = 1, grid%num_t_node(inode)
!print *, 'before itr calculation'
        itr = grid%tr_node(i,inode)
!        area_sum=area_sum+grid%tr_area(itr)
!print *, 'before k loop'
        nd=grid%nd(:,itr)
        dx1= ( grid%x(:,nd(1)) - grid%x(:,nd(3)) )/real(nsub,8)
        dx2= ( grid%x(:,nd(2)) - grid%x(:,nd(3)) )/real(nsub,8)
        area= 0.5D0*abs( dx1(1)*dx2(2) - dx2(1)*dx1(2) )
        do j=1, nsub
           do k=1, 2*j-1 ! for all subtriangle
              itmp=(k-1)/2
              c1=(j-1)- itmp + real(mod(k,2)*2-1,8)/3D0
              itmp= k/2
              c2=itmp + real(mod(k,2)*2-1,8)/3D0
              c3=real(nsub,8) - c1 - c2

              xc=grid%x(:,nd(3)) + c1*dx1 + c2*dx2 ! center of subtriangle
              psi2=psi_interpol(xc(1),xc(2),0,0)
              if(psi2>eq_x_psi*0.7.and.psi2<=sml_outpsi)then
              n_ave=n_ave+neutral_den(xc(1),xc(2),psi2)
              n_tot=n_tot+1D0
              endif
           enddo
        enddo
     enddo
     n_neut(inode)=n_ave/n_tot

!if(idebug==1) print *, 'n_tot=', n, n_neut(inode), n_tot, psi, psi2
   
     P_ion = 0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS 
!     fion(inode)=n*0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS 
!     fion(inode)=n*degas2_atomic_data('D','e','hionize',      &
!                      'reaction_rate',n_p(inode),sml_e_charge,Te)
     rate_neu_ion = P_ion*sml_dt*sml_f_source_period
     fion(inode)=rate_neu_ion*n_neut(inode)
if(is_nan(fion(inode)))print *,'fion=', P_ion, n_neut(inode), Te, Te_ev(inode), f0_T_ev(inode,isp)
!     if(neu_ionize_mode/=2 .OR. psi < neu_ionize2_psi) then 
!     if(psi>eq_x_psi)then
!        neui_ionized=neui_ionized+fion(inode)*n_p(inode)*sml_dt*sml_f_source_period*grid%node_vol_nearest(inode)
!     endif
     if(psi>eq_x_psi*0.7.and.psi<=sml_outpsi)then
     do imu=f0_imu1, f0_imu2
        do ivp= -f0_nvp, f0_nvp
!           delta_fi_neut_ion(ivp,inode,imu) = rate_neu_ion*n_p(inode)*n_neut(inode)*maxwell_neut(ivp,inode,imu,1)
           delta_fi_neut_ion(ivp,inode,imu) = n_p(inode)*fion(inode)*maxwell_neut(ivp,inode,imu,1)
!if(is_nan(delta_fi_neut_ion(ivp,inode,imu)))print *, maxwell_neut(ivp,inode,imu,1),n_p(inode),fion(inode),n_neut(inode)
           neui_ionized=neui_ionized+delta_fi_neut_ion(ivp,inode,imu)*f0_grid_vol(inode,1)
        enddo !ivp
     enddo ! imu
     endif ! psi boundary setup

  enddo ! inode
  call t_stopf("NEUTRAL_IONIZE_PROB")

  call t_startf("NEUTRAL_IONIZE_PROB_RED")
!if(idebug==1) print *, 'delta_fi_neut_ion=', maxval(delta_fi_neut_ion), maxval(delta_fe_neut_ion) 
  do ith=2,sml_nthreads
     neu_weight_sum_lost(1)=neu_weight_sum_lost(1)+neu_weight_sum_lost(ith)
  enddo
  dum1(1) = neui_ionized
  dum1(2) = neue_ionized
  dum1(3) = neu_weight_sum_lost(1)

  call my_mpi_allreduce(dum1,dum2,3)
  sum_neui_ionized = dum2(1)
  sum_neue_ionized = dum2(2)
  sum_neu_weight_lost = dum2(3)
  neu_weight_sum_lost=0D0 !initilize for next use

  ! set maximum of lost particle numbers
  sum_neu_weight_lost = min(sum_neu_weight_lost, neu_weight_lost_max)
  call t_stopf("NEUTRAL_IONIZE_PROB_RED")

!  neu_actual_accumi_ionized = neu_actual_accumi_ionized + dum2(1)
!  neu_actual_accume_ionized = neu_actual_accume_ionized + dum2(2)
!  neu_weight_accum_ionized = neu_weight_accum_ionized + dum2(3)

  
  call t_startf("NEUTRAL_OTHER")
!  if(sml_mype==0)print *, 'neu_weight_lost,sum_neu_ionized=',sum_neui_ionized, sum_neue_ionized, sum_neu_weight_lost
     ! set new ionization normalization according to given recycle rate
     if(neu_adjust_n0==1) then
        ! Finding neutral density for given recycle rate
        new_n0=neu_base_den*sum_neu_weight_lost/(sum_neui_ionized+1D-10)*neu_recycle_rate
if(sml_mype==0.and.idebug==1)print *,'new_n0=', new_n0, sum_neui_ionized, sum_neu_weight_lost, neu_base_den

        !if(istep>=neu_start_time .and. istep < neu_start_time +neu_ion_period) then
        !   new_n0=min(1.d15, new_n0)
        !endif

        !if(sml_mype==0) then
        !   write(62,*) sml_dt*istep, neu_n0,new_n0&
        !        ,weight_sum_out/weight_sum_ionize*neu_recycle_rate,neu_mfp
        !endif
        n_neut=n_neut*new_n0/neu_base_den
        fion=fion*new_n0/neu_base_den
        neu_base_den= new_n0 + 1D-10 ! 1D-10 is for avoiding 'devide by zero'
     endif

  delta_fi_neut_ion=0D0
  neui_ionized=0D0
  dum1=0D0
  sum_neui_ionized=0D0
  do inode = f0_inode1, f0_inode2

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle

     r = grid%x(1,inode) ! get major radius
     z = grid%x(2,inode)
     psi = grid%psi(inode)
!     Te=Te_neutral_col_3(psi,r,z) ! eV
!     n=neutral_den(r,z,psi) ! normalized unit
!     Te = Te_ev(inode) ! ev
!     P_ion = 0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS 
!     rate_neu_ion = P_ion*sml_dt*sml_f_source_period
     if(psi>eq_x_psi*0.7.and.psi<=sml_outpsi)then
       do imu=f0_imu1, f0_imu2
          do ivp= -f0_nvp, f0_nvp
             delta_fi_neut_ion(ivp,inode,imu) = n_p(inode)*fion(inode)*maxwell_neut(ivp,inode,imu,1)
             neui_ionized=neui_ionized+delta_fi_neut_ion(ivp,inode,imu)*f0_grid_vol(inode,1)
          enddo
       enddo
     endif ! psi boundary
  enddo
  call t_stopf("NEUTRAL_OTHER")

  call t_startf("NEUTRAL_OTHER_RED")
  dum1(1) = neui_ionized
  call my_mpi_allreduce(dum1,dum2,3)
  sum_neui_ionized = dum2(1)
  if(sml_mype==0.and.idebug==1)print *,'ion ionization=', sum_neui_ionized
  call t_stopf("NEUTRAL_OTHER_RED")

!! electron from ionization
  if(sml_electron_on) then
    call t_startf("NEUTRAL_ELEC1")
    neue_loss=0D0
    neue_hot=0D0
    neue_cold=0D0
    sum_neue_lost=0D0
    sum_neue_h=0D0
    sum_neue_c=0D0
    delta_fe_neut_ion_h=0D0
    delta_fe_neut_ion_l=0D0
    delta_fe_neut_ion_c2=0D0
    en_neut=0D0
    ene_neut=0D0

!    norm_c=0D0
!$OMP PARALLEL DO &
!$OMP PRIVATE ( inode, rgn_check, r, z, psi, imu, ivp, smu, mu_n, vp_n, ene_neut, imup, ivpp, smup, mu_np, vp_np, &
!$OMP           en_neut, Tc, en_hot, mu_nh, smuh, imu_h, vp_nh, ivp_h, w1, w2,f0_f_tmp)
    do inode = f0_inode1, f0_inode2

       !private region check
       rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
       if(.not. rgn_check) cycle

       !sml_inpsi, sml_outpsi, wall check
       if(sml_inpsi > grid%psi(inode) .or. sml_outpsi < grid%psi(inode) &
          .or. grid%rgn(inode)==grid_rgn_wall ) cycle


       r = grid%x(1,inode) ! get major radius
       z = grid%x(2,inode)
       psi = grid%psi(inode)
!       Te=Te_neutral_col_3(psi,r,z) ! eV
!       P_ion = 0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS

       if(psi>eq_x_psi*0.7.and.psi<=sml_outpsi)then
       do imu=f0_imu1, f0_imu2
          do ivp= -f0_nvp, f0_nvp
              smu=(imu*f0_dsmu) ! mu or sqrt mu (normalized)
              mu_n = smu*smu
              vp_n = ivp*f0_dvp
              ene_neut = 0.5D0*(mu_n+vp_n**2)*f0_T_ev(inode,0)

!! cold electron product
              
              do imup=f0_imu1, f0_imu2
                 do ivpp= -f0_nvp, f0_nvp
                    smup=(imup*f0_dsmu) ! mu or sqrt mu (normalized)
                    mu_np = smup*smup
                    vp_np = ivpp*f0_dvp
                    en_neut = 0.5D0*(mu_np + vp_np**2)*f0_T_ev(inode,0)
                    if(en_neut>=1.5*E_ion)then
                      Tc=E_ion/6D0
                    else
                      Tc=0.5D0*en_neut-E_ion/3D0
                    endif
                    Tc=max(Tc,3D0) ! set a minmum temperature >=3.0 ev
                    if(ene_neut<6D0*Tc)then
                      f0_f_tmp=max(f0_f(ivpp,inode,imup,0),0D0)

                      delta_fe_neut_ion_c2(ivp,inode,imu)=delta_fe_neut_ion_c2(ivp,inode,imu)+fion(inode)* &
                        f0_f_tmp*exp(-ene_neut/Tc)*sqrt(mu_n*f0_T_ev(inode,0)/Tc)/Tc* &
                        f0_grid_vol(inode,0)/grid%node_vol_nearest(inode)
!                       norm_c(ivpp,inode,imup)=norm_c(ivpp,inode,imup)+exp(-ene_neut/Tc)*2D0*sqrt(mu_n*f0_T_ev(inode,0)/Tc)/Tc*f0_grid_vol(inode,0)/grid%node_vol_nearest(inode)
                    endif
                 enddo ! ivpp
              enddo ! imup

!! warm electron product
              en_hot = ene_neut-E_ion  ! need to use degas2 data
              en_hot=max(en_hot, 1D0)
              en_hot=min(en_hot, ene_neut)
              vp_nh = sqrt(en_hot*2D0/f0_T_ev(inode,0))*(2D0*ranx()-1D0)
              mu_nh = (en_hot/f0_T_ev(inode,0)-0.5D0*vp_nh**2)*2D0 !!rh this is now (v_perp/v_th)^2
              smuh = sqrt(mu_nh)
              imu_h = floor(smuh/f0_dsmu)
              ivp_h = floor(vp_nh/f0_dvp)

              if(imu_h<f0_imu1.or.ivp_h<-f0_nvp.or.imu_h+1>f0_imu2.or.ivp_h+1>f0_nvp) then
                 smuh=sqrt(en_hot/(ene_neut+1D-10)*mu_n) 
                 imu_h=floor(smuh/f0_dsmu)
                 vp_nh=sqrt(en_hot/(ene_neut+1D-10))*vp_n 
                 ivp_h=floor(vp_nh/f0_dvp)
              endif

              w1=(smuh-imu_h*f0_dsmu)/f0_dsmu
              w2=(vp_nh-ivp_h*f0_dvp)/f0_dvp
              f0_f_tmp=max(f0_f(ivp,inode,imu,0),0D0)
              delta_fe_neut_ion_h(ivp_h,inode,imu_h)=delta_fe_neut_ion_h(ivp_h,inode,imu_h)+fion(inode) &
                 *f0_f_tmp*(1D0-w1)*(1D0-w2)
              delta_fe_neut_ion_h(ivp_h+1,inode,imu_h)=delta_fe_neut_ion_h(ivp_h+1,inode,imu_h)+fion(inode) &
                 *f0_f_tmp*(1D0-w1)*w2
              delta_fe_neut_ion_h(ivp_h,inode,imu_h+1)=delta_fe_neut_ion_h(ivp_h,inode,imu_h+1)+fion(inode) &
                 *f0_f_tmp*w1*(1D0-w2)
              delta_fe_neut_ion_h(ivp_h+1,inode,imu_h+1)=delta_fe_neut_ion_h(ivp_h+1,inode,imu_h+1)+fion(inode) &
                 *f0_f_tmp*w1*w2

!! Lost electrons
             delta_fe_neut_ion_l(ivp,inode,imu)=fion(inode)*f0_f_tmp

          enddo ! imu
       enddo ! ivp

       endif ! psi
    enddo ! inode

! For diagnosis: to sum up the three parts of electrons individually for checkup
    do inode = f0_inode1, f0_inode2

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle
     !sml_inpsi, sml_outpsi, wall check
     if(sml_inpsi > grid%psi(inode) .or. sml_outpsi < grid%psi(inode) &
         .or. grid%rgn(inode)==grid_rgn_wall ) cycle

       do imu=f0_imu1, f0_imu2
          do ivp= -f0_nvp, f0_nvp
             neue_cold=neue_cold+delta_fe_neut_ion_c2(ivp,inode,imu)*f0_grid_vol(inode,0)
             neue_hot=neue_hot+delta_fe_neut_ion_h(ivp,inode,imu)*f0_grid_vol(inode,0)
             neue_loss=neue_loss+delta_fe_neut_ion_l(ivp,inode,imu)*f0_grid_vol(inode,0)
          enddo ! imu
       enddo ! ivp
    enddo ! inode
    call t_stopf("NEUTRAL_ELEC1")

    call t_startf("NEUTRAL_ELEC1_RED")
    dum1(1) = neue_cold
    dum1(2) = neue_hot
    dum1(3) = neue_loss
    call my_mpi_allreduce(dum1,dum2,3)
    sum_neue_c = dum2(1)
    sum_neue_h = dum2(2)
    sum_neue_lost = dum2(3)
    call t_stopf("NEUTRAL_ELEC1_RED")

if(sml_mype==0.and.idebug==1)print *,'neue_ionization', sum_neue_c, sum_neue_h, sum_neue_lost

! Scale the ionization distribution to be the same for ions and electrons (three parts)
   delta_fi_neut_ion=delta_fi_neut_ion*sum_neue_lost/sum_neui_ionized
   delta_fe_neut_ion_c2=delta_fe_neut_ion_c2*sum_neue_lost/sum_neue_c
         
 endif ! sml_electron

  call t_startf("NEUTRAL_CX1")
  !! Charge exchange----------------------------------------------------------------
!  if(mod(istep,sml_f_source_period)==0) then
  allocate(delta_fi_neut_cx(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2), sum_delta_fi(f0_inode1:f0_inode2), &
         delta_fi_neut_cx_M(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2), sum_delta_fi_M(f0_inode1:f0_inode2))

  delta_fi_neut_cx=0D0
  delta_fi_neut_cx_M=0D0

  sum_ion_cx=0D0
  sum_loss_cx=0D0

  sum_delta_fi=0D0
  sum_delta_fi_M=0D0

  do inode = f0_inode1, f0_inode2  

     !private region check
     rgn_check=(grid%rgn(inode)==1 .or. grid%rgn(inode)==2)
     if(.not. rgn_check) cycle

     r = grid%x(1,inode) ! get major radius
     z = grid%x(2,inode)
     psi = grid%psi(inode)
!     Te=Te_neutral_col_3(psi,r,z) ! eV
     n=neutral_den(r,z,psi) ! normalized unit
     B=B_interpol(r,z,0D0) 
     ! Get ion/electron temperature
     if(sml_electron_on)then
     enev = 1.5D0*Te_ev(inode) !psn%tempe_ev(inode)
     else
     enev = 1.5D0*f0_T_ev(inode,1) !psn%tempi_ev(inode)
     endif
     
     ! probability of charge exchange per sec
     P_cx= n*1.1D-8* (enev)**0.3 /sqrt(ptl_mass_au) *1D-6 
     rate_neu_col=P_cx*sml_dt*sml_f_source_period
!     sum_delta_fi(inode)=0D0
!     sum_delta_fi_M(inode)=0D0
     ! lost ions by charge exchange
     do imu=f0_imu1, f0_imu2
        do ivp= -f0_nvp, f0_nvp
           f0_f_tmp=max(f0_f(ivp,inode,imu,1),0D0)
           delta_fi_neut_cx(ivp,inode,imu) = -rate_neu_col*f0_f_tmp
           sum_delta_fi(inode)=sum_delta_fi(inode)+delta_fi_neut_cx(ivp,inode,imu)*f0_grid_vol(inode,1)
!           if(z>eq_x_z) then
             sum_ion_cx = sum_ion_cx+delta_fi_neut_cx(ivp,inode,imu)*f0_grid_vol(inode,1)
!           endif
        enddo
     enddo
     ! ion birth by charge exchange -- it has neutral energy distribution
     do imu=f0_imu1, f0_imu2
        do ivp= -f0_nvp, f0_nvp
           delta_fi_neut_cx_M(ivp,inode,imu) = rate_neu_col*n_p(inode)*maxwell_neut(ivp,inode,imu,1)
           sum_delta_fi_M(inode)=sum_delta_fi_M(inode)+delta_fi_neut_cx_M(ivp,inode,imu)*f0_grid_vol(inode,1)
!           if(z>eq_x_z) then
             sum_loss_cx = sum_loss_cx+delta_fi_neut_cx_M(ivp,inode,imu)*f0_grid_vol(inode,1)
!           endif
        enddo
     enddo
    ! enforce ion birth by CX to be the same as ion lost by CX for each space cell
     delta_fi_neut_cx_M(:,inode,:) = -delta_fi_neut_cx_M(:,inode,:)*sum_delta_fi(inode)/(sum_delta_fi_M(inode)+1D-45)
  enddo 
  delta_fi_neut_cx = delta_fi_neut_cx + delta_fi_neut_cx_M
  call t_stopf("NEUTRAL_CX1")

  call t_startf("NEUTRAL_CX1_RED")
  dum1(1) = sum_ion_cx
  dum1(2) = sum_loss_cx
  call my_mpi_allreduce(dum1,dum2,3)
  sum_ion_cx = dum2(1)
  sum_loss_cx = dum2(2)
  if(sml_mype==0)print *,'charge exchange', sum_ion_cx, sum_loss_cx
  call t_stopf("NEUTRAL_CX1_RED")
  
  ! update df0g after neutral collisions
  f0_df0g(:,:,:,1) = f0_df0g(:,:,:,1) +delta_fi_neut_ion(:,:,f0_imu1:f0_imu2)+ delta_fi_neut_cx(:,:,f0_imu1:f0_imu2)
  if(sml_electron_on) then
     f0_df0g(:,:,:,0) = f0_df0g(:,:,:,0) -delta_fe_neut_ion_l(:,:,:)+delta_fe_neut_ion_h(:,:,:)+ &
                       delta_fe_neut_ion_c2(:,:,:)
!     deallocate(delta_fe_neut_ion)
  endif

  deallocate(delta_fi_neut_ion,delta_fi_neut_cx,sum_delta_fi, maxwell_neut,  &
             delta_fi_neut_cx_M,sum_delta_fi_M, n_p, n_p2, delta_fe_neut_ion_l, &
             delta_fe_neut_ion_h, delta_fe_neut_ion_c2, fion, eperp1, epara, Te_ev )


end subroutine neutral3



real (kind=8) function te_neutral_col_3(psi,r,z)
! electron temperature in eV
  use eq_module
  implicit none
  real (kind=8), external :: tempe_ev
  real (kind=8) :: psi,r,z,xd

     te_neutral_col_3=eq_ftn(psi,r,z,eq_tempe) !tempe_ev(psi,z,0)

end function te_neutral_col_3


real (kind=8) function neutral_den1_3(r,z) !return nuetral density in normalized unit
  use neu_module
  use eq_module, only : eq_x_psi, eq_x_z, eq_axis_r,eq_axis_z
  use sml_module, only : sml_2pi
  implicit none
  real (kind=8) , intent(in) :: r,z
  real (kind=8) , external :: psi_interpol
  real (kind=8) :: theta, n_theta, psi,n_r, alpha, r_s, z_s, d,del
  integer :: i

  ! n_theta
  ! finding ther value of (r,z) position
  theta=acos((r-eq_axis_r)/dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2))
  if(z<eq_axis_z) then
     theta=sml_2pi-theta
  endif
  ! Del value -- angle distance from x-point
  del=min(abs(theta-neu_peak_theta),min( abs(theta-neu_peak_theta+sml_2pi), abs(theta-neu_peak_theta-sml_2pi)))
  ! theta dependance of density
  n_theta=neu_base_den*(    1D0 + (neu_delta_n-1D0)*exp( -(del/neu_delta_theta)**2 )    )

  ! n_r
  ! radial dependance of density - distance from separatrix
  psi=psi_interpol(r,z,0,0)
  if(psi>eq_x_psi) then
!     n_r=0.9D0*neu_nr !n_r=1
     n_r=1
  else  ! distance from separatrix, ignoring non-orthogonal effect
     if(z>eq_x_z) then
        i=int(theta/sml_2pi*real(neu_sep_mtheta))+1
        i=min(neu_sep_mtheta-1,max(1,i))
        alpha=theta/sml_2pi*real(neu_sep_mtheta) + 1 - i
        r_s= (1D0-alpha)*neu_sep_r(i) + alpha*neu_sep_r(i+1)
        z_s= (1D0-alpha)*neu_sep_z(i) + alpha*neu_sep_z(i+1)
        d=dsqrt( (r-r_s)**2 + (z-z_s)**2 )
        n_r=dexp( - d/neu_mfp  )
     else
!        n_r=0.9D0*neu_nr !n_r=1
        n_r=neu_nr
     endif
  endif

  neutral_den1_3=n_theta*n_r
end function neutral_den1_3




real (kind=8) function neutral_den_3(r,z,psi_in)
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8) ,external :: neutral_den1,neu_psi_interpol
  real (kind=8), intent(in) :: r , z ,psi_in
  integer :: ipsi,itheta,itheta_p1
  real (kind=8) :: a1,a2,theta,psi

  if(neu_grid_mode>0) then
     psi=neu_psi_interpol(r,z,psi_in)  ! evaluate effective psi for neutral simulation
  else
     psi=psi_in
  endif

  if(neu_col_mode==1) then
     neutral_den_3=neutral_den1(r,z)
  else
     if((neu_grid_min_psi<psi).and.(psi<neu_grid_max_psi).and.((z>eq_x_z).or. &
        (neu_grid_mode>0))) then
        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z<eq_axis_z) then
           theta= sml_2pi - theta
        endif
        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        a1=max(0D0,min(1D0,a1))
        itheta= int(theta/neu_grid_dtheta) + 1
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1

        neutral_den_3=neu_grid_den(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_den(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_den(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_den(ipsi+1,itheta_p1)*a1*a2
        neutral_den_3=neutral_den_3*neu_base_den  ! propotional to neu_n0

     elseif((psi<=neu_grid_min_psi).and.(z>eq_x_z)) then
        neutral_den_3=1D-10  ! very small number instead of zero
     else
        neutral_den_3=neutral_den1(r,z)
     endif
     neutral_den_3=max(neutral_den_3,1D-10)  ! very small lower bound for safety
  endif
end function neutral_den_3



real (kind=8) function neutral_temp_3(r,z,psi_in)
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8),external :: init_tempi_ev,tempi_ev,neu_psi_interpol
  real (kind=8),intent(in) :: r,z,psi_in
  integer :: ipsi,itheta,itheta_p1
  real (kind=8) :: a1,a2,theta,psi

  if(neu_grid_mode>0) then
     psi=neu_psi_interpol(r,z,psi_in)  ! evaluate effective psi for neutral simulation
  else
     psi=psi_in
  endif

  if(neu_col_mode==1) then
     neutral_temp_3=eq_ftn(psi,r,z,eq_tempi)*neu_temp_factor*sml_ev2j
  else ! mode 2
     if((psi>neu_grid_min_psi).and.(psi<neu_grid_max_psi).and.((z>eq_x_z).or. &
        (neu_grid_mode>0))) then
        theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
        if(z < eq_axis_z) then
           theta= sml_2pi - theta
        endif
        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        a1=max(0D0,min(1D0,a1))

        itheta= int(theta/neu_grid_dtheta) + 1
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1

        neutral_temp_3=neu_grid_temp(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_temp(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_temp(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_temp(ipsi+1,itheta_p1)*a1*a2
     else
!        neutral_temp=neu_temp0*sml_ev2norme
        neutral_temp_3=eq_ftn(psi,r,z,eq_tempi)*sml_ev2j
     endif
  endif

end function neutral_temp_3






