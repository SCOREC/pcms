!*****************************************************************************
! neutral collision (charge exchange, ionization)
! first created 2002/09/26
! last modified 2011/11/23
! This routine assumes initial maximum number and electron number is the same.
! shodul be generalized later.
!*****************************************************************************

subroutine neutral_col(istep,grid,psn,sp,sp2)
  use grid_class
  use psn_class
  use ptl_module
  use neu_module
  use sml_module
  use eq_module, only : eq_x_psi,eq_x_z
  use col_module, only: col_mode,col_varying_bg,col_vb_period
  use perf_monitor
  use random_xgc
  implicit none
  include 'mpif.h'
  type(species_type):: sp, sp2
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type) ::  ptl(sp%maxnum) 
  integer ,intent(in) :: istep
  real (kind=8),external :: b_interpol,neutral_den,psi_interpol,te_neutral_col,neutral_temp,neutral_flow
!  real (kind=8),external :: outptl
  real (kind=8) :: r,z,phi,rho,mu,b,n,f_cx,enev,pitch,Te,psi,t0,rate_neu_col
  real (kind=8) ,allocatable::  f_ion(:)
!  integer , allocatable ::pindex(:)
  integer (kind=8) :: dumi, num_ionized, newgid_start, dum4(neu_grid_mpsi),i,k,l,ith
  integer ::  old_ptl_num,old_ptl_ion_num,  ierror, ipsi
  integer (kind=8), allocatable :: global_ionized(:)
  real (kind=8) :: zmax, zdum, v, dum , dum1(3), dum2(3), dum3(neu_grid_mpsi), &
       weight_sum_out,weight_sum_ionize,new_n0,weight_sum_inside,weight_actual_ionize
  real (kind=8) :: coll,agg,ekin
  real (kind=8) :: t00,t11
  real (kind=8) :: En_ev,alpfac,sqrt_alpfac,en_j
  target ptl

  ! evaluate self-consistent neutral profiles by tracking neutral particles in time under evolved plasma profiles

  if (neu_col_mode /= 0) then
     if(col_varying_bg==1 .and. (istep==1 .or.  mod(istep,neu_mode2_period)==0)) then
        call t_startf("NEUTRAL_COL_SNAPSHOT")
        call col_snapshot(sp)
        call t_stopf("NEUTRAL_COL_SNAPSHOT")
     endif
  endif

  if((neu_col_mode==2).AND.(mod(istep,neu_mode2_period)==0 .OR. istep==1)) then
     ! neutral diffusion calculation

     call t_startf("NEUTRAL_STEP")
     call neutral_step(ptl_mass(1))
     call t_stopf("NEUTRAL_STEP")
  endif


  if(istep < neu_start_time) return ! no neutral effect before a couple profile evaluations

  ! neu_ion_period and neu_col_period and neu_cx_period should be the same #### make change later


  ! algorithm 2
  !! Ionization ----------------------------------------------------------------
  if(mod(istep,neu_ion_period)==0) then
     ! print *,'ionize',istep,neu_ion_period
     ! ionization due to electron - neutral collision
     ! recycling gid==-gid ptl
     ! choose a ptl randomly
!     call check_point('inside the ionization 0')
!     print *, 'memory_cleaning', sml_mype
!     call memory_cleaning  ! eliminate all gid(i)<0 particles and find the sum of weights
!     print *, 'memory_cleaning end', sml_mype
!call check_point('before ionization')
     if(neu_varying_mfp==1) then  ! Change mean free path
        call change_neu_mfp(sp) ! set new mean free path
     endif

     ! initialize
!     outptl_num=0
!     empty_num=0
     allocate(f_ion(sp2%num))
     weight_sum_ionize=0D0 ! weight sum of ionize ptl (expectation value)
     weight_actual_ionize=0D0 ! weight sum of ionize ptl (actual value)
     f_ion=0D0     ! ionization probability
     weight_sum_inside=0D0
     neu_actual_ionized=0D0
     ! calculate ionization probability
     do i=1, sp2%num
        if(sp2%ptl(i)%gid>0) then
           ! weight summation of ionize ptl (expectation value)
           r=sp2%ptl(i)%ph(1)
           z=sp2%ptl(i)%ph(2)
           psi=psi_interpol(r,z,0,0)
           if(psi<neu_ionize2_psi) weight_sum_inside=weight_sum_inside+sp2%ptl(i)%ct(piw0)
           Te=Te_neutral_col(psi,r,z) ! eV
           n=neutral_den(r,z,psi) ! normalized unit
           f_ion(i)= n * 0.8D-8*sqrt(Te)*(exp(-13.56/Te)/(1D0+0.01D0*Te))*1D-6  ! MKS !bugfix 1D-6
           !          print *, i,n,sqrt(Te), exp(-13.56/Te)
           f_ion(i)=f_ion(i) ! s^-1 --> normalized unit  , unit conversion
           if(neu_ionize_mode/=2 .OR. psi < neu_ionize2_psi) then ! exclude ionize_mode==2 .and. psi > neu_ionize2_psi
              weight_sum_ionize=weight_sum_ionize+(1D0-dexp(-f_ion(i)*sml_dt*neu_ion_period))*sp2%ptl(i)%ct(piw0)
           endif

        endif
     enddo
     ! MPI summation
     call t_startf("NEUTRAL_COL_SUM")
     !OpenMP summation
     do ith=2,sml_nthreads
        neu_weight_sum_lost(1)=neu_weight_sum_lost(1)+neu_weight_sum_lost(ith)
     enddo
     if(neu_ionize_mode==1) then        
        dum1(1)=neu_weight_sum_lost(1)
        dum1(3)=dum1(1) ! for diagnosis only
        neu_weight_sum_lost=0D0 !initilize for next use
     else  ! ionize mode 2
        dum1(1)=max(neu_old_inside_weight-weight_sum_inside,1D-50)
        neu_old_inside_weight=weight_sum_inside
        dum1(3)=neu_weight_sum_lost(1) !for diagnosis ---
        neu_weight_sum_lost=0D0
     endif
     dum1(2)=weight_sum_ionize
     call my_mpi_allreduce(dum1,dum2,3)
     weight_sum_out=dum2(1)  !if neu_ionize_mode /=1 weight_sum_out is 'difference of weight sum inside sep'
     neu_weight_accum_lost=neu_weight_accum_lost+weight_sum_out
     weight_sum_ionize=dum2(2)
     neu_weight_accum_ionized=neu_weight_accum_ionized+weight_sum_ionize
     call t_stopf("NEUTRAL_COL_SUM")

     ! set new ionization probability f_ion according to given recycle rate
     if(neu_adjust_n0==1) then
        ! Finding neutral density for given recycle rate
        new_n0=neu_base_den*weight_sum_out/(weight_sum_ionize+1D-10)*neu_recycle_rate
        if(istep>=neu_start_time .and. istep < neu_start_time +neu_ion_period) then
           new_n0=min(1.d15, new_n0)
        endif
        !if(sml_mype==0) then
        !   write(62,*) sml_dt*istep, neu_n0,new_n0&
        !        ,weight_sum_out/weight_sum_ionize*neu_recycle_rate,neu_mfp
        !endif

        ! set new ionization probability
        f_ion=f_ion*new_n0/neu_base_den
        neu_base_den=new_n0 + 1D-10 ! 1D-10 is for avoiding 'devide by zero'
     endif

     ! create electron-ion pair according to the ionization probability calculated above
     i=1 ! ptl index
     old_ptl_num=sp2%num
     if(sml_electron_on) old_ptl_ion_num=sp%num
     do while(i<=old_ptl_num ) ! gid=-gid ptl -> f_ion(i)=0
        rate_neu_col=1D0-dexp(-f_ion(i)*sml_dt*neu_ion_period)
        if(rate_neu_col>ranx()) then
           if((f_ion(i)*sml_dt*neu_ion_period)>1.5D0 .AND. sml_mype==0 .AND. sml_mstep/=1) then
              print *, sml_mype, f_ion(i)*sml_dt*neu_ion_period , 'Too large dt in ionization'
           endif

           l=sp2%num+1
           sp2%num=sp2%num+1
           if (sp2%num > sp2%maxnum) then
             print *, sml_mype , 'not enough memory space: increase ptl_maxnum'
             print *, real(i)/real(sp2%num),neu_base_den
             stop
           endif

           r=sp2%ptl(i)%ph(1)
           z=sp2%ptl(i)%ph(2)
           psi=psi_interpol(r,z,0,0)
           phi=sp2%ptl(i)%ph(3)
           rho=sp2%ptl(i)%ph(4)
           mu=sp2%ptl(i)%ct(pim)
           B=b_interpol(r,z,0D0) !### use phi instad of 0D0? 


           ! kinetic energy cooling for ionizing electrons
           if(sml_electron_on) then
              En_j=mu*B + ptl_c2_2m(sp2%type)* (rho*B)**2
              En_ev=En_j*sml_j2ev
              alpfac=1D0-neu_ionize_elec_cooling_ev/En_ev
              alpfac=max(alpfac,1D-2)
              sqrt_alpfac=sqrt(alpfac)
              sp2%ptl(i)%ph(4)=sp2%ptl(i)%ph(4)*sqrt_alpfac
              sp2%ptl(i)%ct(pim)=sp2%ptl(i)%ct(pim)*alpfac
           endif

           ! set newly born ptl's position & weight
           sp2%ptl(l)%ph(1)=r
           sp2%ptl(l)%ph(2)=z
           sp2%ptl(l)%ph(3)=phi
           sp2%ptl(l)%ct(piw0)=sp2%ptl(i)%ct(piw0)
           if(neu_ionize_mode/=2 .or. psi < neu_ionize2_psi) then
              neu_actual_ionized=neu_actual_ionized+sp2%ptl(l)%ct(piw0)
           endif
           sp2%ptl(l)%ph(piw1:piw2)=0D0  ! delta-f scheme is not implimented for neutral collision
  
           ! ionize -- 2002/10/09 same energy ionization inside separatrix
           b=b_interpol(r,z,phi)
           if((psi<eq_x_psi).and.(neu_col_mode==1)) then
              call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,sp%type)
              pitch=2D0*ranx() - 1D0
              if(sml_mype==0.and.mu.lt.0)print *, 'charge ionization -mu found'
              call ev_pitch_to_rho_mu2(enev,pitch,b,rho,mu,sp%type)
           else  ! Maxwellian - outside of separatrix or mode2
              if(sml_electron_on) then
                 t0=neutral_temp(r,z,psi)
                 ! perp velocity
                 zmax=1.d0 - dexp(-7.d0)
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0/ptl_mass(sp2%type)*dlog(1.d0-zdum)*t0)
                 mu=0.5D0*ptl_mass(sp2%type)*(v**2)/B
                 ! parallel velocity
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0/ptl_mass(sp2%type)*dlog(1.d0-zdum)*t0)
                 v= v*dcos(sml_pi*ranx())
                 ! newly generated ions are given neutral parallel flow
                 if(neu_flow_on==1) v=v+neutral_flow(r,z,psi)
                 rho=v/b*ptl_mass(sp2%type)/ptl_charge(sp2%type)
              else
                 ! set new electron energy corresponding to cooling energy of the ionizing electrons
                 t0=neu_ionize_elec_cooling_ev-13.56D0 ! ionization energy is subtracted
                 pitch=2D0*ranx() - 1D0
                 call ev_pitch_to_rho_mu2(t0,pitch,b,rho,mu,sp2%type)
              endif
           endif
	   
           sp2%ptl(l)%ph(4)=rho
           sp2%ptl(l)%ct(pim)=mu
           
           !neutral diagnosis - ionization
           call diag_neu_ionization(sp2%ptl(l),psi,b,sp2%type)

!           if(sml_electron_on) sp2%reflec_on(i)=0
!           if(sml_electron_on) sp2%reflec_on(l)=0

           ! create new ion with same energy and temp as ionized neutral.
           if(sml_electron_on) then
              ! perp velocity
              t0=neutral_temp(r,z,psi)
              zmax=1.d0 - dexp(-7.d0)
              zdum=zmax*ranx()
              v= dsqrt(-2.d0/ptl_mass(1)*dlog(1.d0-zdum)*t0)
              mu=0.5D0*ptl_mass(1)*v**2/B
              ! parallel velocity
              zdum=zmax*ranx()
              v= dsqrt(-2.d0/ptl_mass(sp2%type)*dlog(1.d0-zdum)*t0)
              v=v*dcos(sml_pi*ranx())
              if(neu_flow_on==1) v=v+neutral_flow(r,z,psi) ! newly generated ions are given neutral parallel flow
              rho=v/b*ptl_mass(sp2%type)/ptl_charge(sp2%type)

              sp%num=sp%num+1
              sp%ptl(sp%num)%ph(:)=sp%ptl(l)%ph(:)  !ptl%ion%phase(:,ptl%ion%num)=ptl%elec%phase(:,l)
              sp%ptl(sp%num)%ph(4)=rho
              sp%ptl(sp%num)%ct(pim)=mu
              sp%ptl(l)%ct(piw0)=sp%ptl(i)%ct(piw0)
!              ptl%ion%phase(4,ptl%elec%num)=ptl%elec%phase(4,l)*sml_charge_rel(2)/sqrt(sml_mass_rel(2))
           endif

!           neu_actual_accum_ionized=neu_actual_accum_ionized+sp2%ptl(l)%ct(piw0)
           weight_actual_ionize=weight_actual_ionize+sp2%ptl(l)%ct(piw0)
           ! increase the total number of ionization events by one
           if(z>eq_x_z) then
              ipsi= int((psi-sml_inpsi)/neu_grid_dpsi2) + 1
              ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
              neu_ion_events(ipsi)=neu_ion_events(ipsi)+1
              neu_ion_events_weights(ipsi)=neu_ion_events_weights(ipsi)+sp2%ptl(i)%ct(piw0)
           endif
        endif
        i=i+1
     enddo

     if(neu_ionize_mode==2) then
        neu_old_inside_weight=neu_old_inside_weight+neu_actual_ionized
     endif

     deallocate(f_ion)
     call check_point('during ionization')
     ! print out n0
     if(sml_mype==0 .and.sml_istep >neu_ion_period  ) then
        write(62,*) sml_time/sml_tran, neu_base_den,dum2(3)*sml_e_charge/(sml_dt*neu_ion_period)
     endif

     call t_startf("NEUTRAL_RED1")
     call MPI_ALLREDUCE(neu_ion_events,dum4,neu_grid_mpsi,MPI_INTEGER8,MPI_SUM,sml_comm,ierror)
     call my_mpi_allreduce(neu_ion_events_weights,dum3,neu_grid_mpsi)
!     call my_mpi_allreduce(neu_actual_accum_ionized,dum,1)
     call my_mpi_allreduce(weight_actual_ionize,dum,1)
     call t_stopf("NEUTRAL_RED1")
     neu_actual_accum_ionized=neu_actual_accum_ionized+dum
     ! print out accumulated total number and weights of ionized particles
     if((sml_mype==0).and.(mod(istep,2*neu_ion_period)==0)) then
        write(63,*) '#',sml_time,sml_tran, sml_time/sml_tran ! insert empty line for gnuplot  - time separate
        write(63,*)
        do i=1, neu_grid_mpsi-1
           write(63,*) (sml_inpsi+(real(i-1)+0.5D0)*neu_grid_dpsi2)/eq_x_psi, dum4(i), dum3(i)
        enddo
     endif

     ! print out total weights of ionized and lost particles
     if((sml_mype==0).and.(mod(istep,2*neu_ion_period)==0)) then
        write(64,1000) sml_time/sml_tran, neu_actual_accum_ionized, neu_weight_accum_ionized, neu_weight_accum_lost, &
              dum/real(neu_ion_period)/(sml_dt), weight_sum_out/real(neu_ion_period)/(sml_dt)
     endif

     ! set global IDs for new particles and increment max global ID
     call t_startf("NEUTRAL_GATSCAT")
     num_ionized = sp2%num - old_ptl_num
     allocate(global_ionized(sml_totalpe))
     ! gather number of new particles on each processor
     call MPI_GATHER(num_ionized,1,MPI_INTEGER8,global_ionized,1,MPI_INTEGER8,0,sml_comm,ierror)
     if (sml_mype==0) then
        dumi = 0
        do i=1,sml_totalpe
           dumi = dumi + global_ionized(i)
        enddo
        sp2%maxgid = sp2%maxgid + dumi
     ! refill global array with newgid_start values and scatter
        global_ionized(sml_totalpe) = sp2%maxgid - global_ionized(sml_totalpe)
        do i=sml_totalpe-1,1,-1
           global_ionized(i) = global_ionized(i+1) - global_ionized(i)
        enddo
     endif
     call MPI_SCATTER(global_ionized,1,MPI_INTEGER8,newgid_start,1,MPI_INTEGER8,0,sml_comm,ierror)
     deallocate(global_ionized)
     ! assign global IDs to new particles
     do i=1,num_ionized
        sp2%ptl(old_ptl_num+i)%gid = newgid_start + i
        if(sml_electron_on) then
           sp%ptl(old_ptl_ion_num+i)%gid = newgid_start + i
        endif
     enddo
     call t_stopf("NEUTRAL_GATSCAT")
     if(neu_update_elec) then
        call t_startf("NEUTRALE")
        call neutrale(grid,psn,sp,old_ptl_num)
        call t_stopf("NEUTRALE")
     endif
  endif ! end of ionization

  !! Elastic Collisions -------------------
  if((mod(istep,neu_col_period)==0).and.(neu_elastic_col_on==1)) then
     ! Pitch angle scattering by neutrals -> change mu and rho
     do i=1, sp%num
        if(sp%ptl(i)%gid > 0) then
           r=sp%ptl(i)%ph(1)
           z=sp%ptl(i)%ph(2)
           psi=psi_interpol(r,z,0,0)
           rho=sp%ptl(i)%ph(4)
           mu=sp%ptl(i)%ct(pim)
           b=b_interpol(r,z,0D0)
           n=neutral_den(r,z,psi)
           if(sml_mype==0.and.mu.lt.0)print *, 'elastic collision -mu found'
           call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,sp%type)
           ekin=enev*sml_ev2j
           !cccccccccccccccccccccccccccccccccccccc
           !ccc   pitch angle scattering for small value of ion-neutral
           !collision frequency, coll
           !ccccc Boozer and Kuo-Petravic Phys. Fluids 24, 851 (1981)
           !C============
           coll=1.55D-13*n*dsqrt(enev/1D3/ptl_mass_au)
           !coll = col_col*(sml_en_order/(ekin))**1.5
           agg = ranx() - .5D0
           dum = 1.D0 - pitch**2
           dum = max(0.D0,dum)
           agg = dsign(1.D0,agg)*dsqrt(dum*coll*sml_dt*neu_col_period*0.5D0)
           pitch = pitch*(1.D0 - coll*sml_dt*neu_col_period*.5D0) + agg  !2002/09
           pitch = min( 0.9999999D0,pitch)
           pitch = max(-0.9999999D0,pitch)

           !cc   change rho ,rmu,  due to scattering

           rho = pitch*dsqrt(2.D0*ekin*ptl_mass(sp%type))/b/ptl_charge(sp%type)
           mu = ekin/b - ptl_c2_2m(sp%type)*rho*rho*b

           !neutral diagnosis - elastic
           call diag_neu_elastic(sp%ptl(i),rho,mu,psi,b,sp%type)


           sp%ptl(i)%ph(4) = rho
           sp%ptl(i)%ct(pim) =mu

        endif
     enddo
  endif

  !! Charge exchange----------------------------------------------------------------
  if(mod(istep,neu_cx_period)==0) then
     !  print *, 'cx',istep,neu_cx_period
     ! charge exchange. -> change ion mu, rho

     do i=1, sp%num
        if(sp%ptl(i)%gid > 0) then
           r=sp%ptl(i)%ph(1)
           z=sp%ptl(i)%ph(2)
           psi=psi_interpol(r,z,0,0)
           rho=sp%ptl(i)%ph(4)
           mu=sp%ptl(i)%ct(pim)
           b=b_interpol(r,z,0D0)
           n=neutral_den(r,z,psi)
           call rho_mu_to_ev_pitch2(rho,mu,b,enev,pitch,sp%type)
           ! probability of charge exchange per sec
           f_cx=( n * 1.1D-8* (enev)**0.3 /sqrt(ptl_mass_au) *1D-6 )
           ! 1.1 K_i(eV)^0.3 * 10^-8 / Sqrt(Mi(AU))  cm s^-1   -> 1D-6 m s^-1
           f_cx = f_cx  ! s^-1 --> normalized unit.
           rate_neu_col=1D0-dexp(-f_cx*sml_dt*neu_cx_period)
!           if(sml_dt*neu_cx_period*f_cx > ranx() ) then
           if(rate_neu_col>ranx()) then
              !Charge exchange -- give random change of pitch, Ti=Tn 2002/10/09 inside separatrix  -> same energy
              ! Ti=/= Tn exactly, because of temp. dependance of f_cx
              if(psi<eq_x_psi .and. neu_col_mode==1) then
                 pitch=2D0*ranx() - 1D0
                 call ev_pitch_to_rho_mu2(enev,pitch,b,rho,mu,sp%type)
              else  ! outside of seapratrix , Maxwellian dist
                 t0=neutral_temp(r,z,psi)
                 !perp velocity
                 zmax=1.d0 - dexp(-7.d0)
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0/ptl_mass(sp%type))
                 mu=0.5D0*v**2*ptl_mass(sp%type)/b
                 ! parallel velocity
                 zdum=zmax*ranx()
                 v= dsqrt(-2.d0*dlog(1.d0-zdum)*t0/ptl_mass(1))
                 v= v*dcos(sml_pi*ranx())
                 if(neu_flow_on==1) v=v+neutral_flow(r,z,psi) ! newly generated ions are given neutral parallel flow
                 rho=ptl_mass(sp%type)/ptl_charge(sp%type)*v/b
              endif


              !neutral diagnosis - charge exchange
              call diag_neu_cx(sp%ptl(i),rho,mu,psi,b,sp%type)

              sp%ptl(i)%ph(4)=rho
              sp%ptl(i)%ct(pim)=mu

              ! increase the total number of charge exchange events by one
              if(z>eq_x_z) then
                 ipsi= int((psi-sml_inpsi)/neu_grid_dpsi2) + 1
                 ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
                 neu_cx_events(ipsi)=neu_cx_events(ipsi)+1
                 neu_cx_events_weights(ipsi)=neu_cx_events_weights(ipsi)+sp2%ptl(i)%ct(piw0)
              endif
           endif

        endif
     enddo

     call t_startf("NEUTRAL_RED2")
     call MPI_ALLREDUCE(neu_cx_events,dum4,neu_grid_mpsi,MPI_INTEGER8,MPI_SUM,sml_comm,ierror)
     call my_mpi_allreduce(neu_cx_events_weights,dum3,neu_grid_mpsi)
     call t_stopf("NEUTRAL_RED2")
     ! print out accumulated total number and weights of charge exchange events
     if((sml_mype==0).and.(mod(istep,2*neu_cx_period)==0)) then
        write(65,*) '#', sml_time/sml_tran ! insert empty line for gnuplot  - time separate
        write(65,*)
        do i=1, neu_grid_mpsi-1
           write(65,*) (sml_inpsi+(real(i-1)+0.5D0)*neu_grid_dpsi2)/eq_x_psi, dum4(i), dum3(i)
        enddo
     endif
  endif

#if defined(ADIOS)
  ! write out neutral diag data
  if(mod(istep,neu_cx_period)==0) then
     
     call diag_neu_output(new_n0)
     
  endif
#endif

1000 format(6(e19.13,1x))
end subroutine neutral_col

real (kind=8) function neutral_den(r,z,psi_in)
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
     neutral_den=neutral_den1(r,z)
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

        neutral_den=neu_grid_den(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_den(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_den(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_den(ipsi+1,itheta_p1)*a1*a2
        neutral_den=neutral_den*neu_base_den  ! propotional to neu_n0

     elseif((psi<=neu_grid_min_psi).and.(z>eq_x_z)) then
        neutral_den=1D-10  ! very small number instead of zero
     else
        neutral_den=neutral_den1(r,z)
     endif
     neutral_den=max(neutral_den,1D-10)  ! very small lower bound for safety
  endif
end function neutral_den

real (kind=8) function neutral_temp(r,z,psi_in)
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
     neutral_temp=eq_ftn(psi,r,z,eq_tempi)*neu_temp_factor*sml_ev2j
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

        neutral_temp=neu_grid_temp(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_temp(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_temp(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_temp(ipsi+1,itheta_p1)*a1*a2
     else
!        neutral_temp=neu_temp0*sml_ev2norme
        neutral_temp=eq_ftn(psi,r,z,eq_tempi)*sml_ev2j
     endif
  endif

end function neutral_temp

real (kind=8) function neutral_flow(r,z,psi_in)
  use neu_module
  use sml_module
  use eq_module
  implicit none
  real (kind=8),external :: neu_psi_interpol
  real (kind=8),intent(in) :: r,z,psi_in
  integer :: ipsi,itheta,itheta_p1
  real (kind=8) :: a1,a2,theta,psi

  if(neu_grid_mode>0) then
     psi=neu_psi_interpol(r,z,psi_in)  ! evaluate effective psi for neutral simulation
  else
     psi=psi_in
  endif

  if(neu_col_mode==1) then
     neutral_flow=0D0
  else ! mode 2
     if((psi<neu_grid_max_psi).and.((z>eq_x_z).or.(neu_grid_mode>0))) then
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

        neutral_flow=neu_grid_flow(ipsi,itheta)*(1D0-a1)*(1D0-a2) + &
             neu_grid_flow(ipsi+1,itheta)*a1*(1D0-a2) + &
             neu_grid_flow(ipsi,itheta_p1)*(1D0-a1)*a2 + &
             neu_grid_flow(ipsi+1,itheta_p1)*a1*a2
     else
        neutral_flow=0D0
     endif
  endif

end function neutral_flow

real (kind=8) function neutral_den1(r,z) !return nuetral density in normalized unit
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

  neutral_den1=n_theta*n_r
end function neutral_den1

real (kind=8) function te_neutral_col(psi,r,z)
! electron temperature in eV
  use eq_module
  implicit none
  real (kind=8), external :: tempe_ev
  real (kind=8) :: psi,r,z,xd

     te_neutral_col=eq_ftn(psi,r,z,eq_tempe) !tempe_ev(psi,z,0)

end function te_neutral_col

subroutine change_neu_mfp(sp)  !set new mean free path
  use ptl_module
  use diag_module 
  use eq_module !, only : eq_x_z,eq_x_psi,eq_den_v1
  use neu_module
  use sml_module, only: sml_mype
  implicit none
  type(ptl_type) :: ptl
  type(species_type) :: sp
  integer :: in, out,i,j
  real (kind=8) :: den(diag_flow_npsi),r,z,psi,dum(diag_flow_npsi),maxden
  real (kind=8) :: psi_interpol

  ! finding density as a function of psi
  den=0D0  ! den is weight sum for a while
  do i=1, sp%num
     r=sp%ptl(i)%ph(1)
     z=sp%ptl(i)%ph(2)

     if(sp%ptl(i)%gid > 0) then
        psi=psi_interpol(r,z,0,0)
        if(diag_1d_pin<psi .AND. psi<diag_1d_pout .AND. z>eq_x_z ) then
           j=int((psi-diag_1d_pin)/diag_1d_dp)+1
           j=min(diag_flow_npsi,max(1,j))

           den(j)=den(j)+sp%ptl(i)%ct(piw0) 

        endif
     endif
  enddo

  call my_mpi_allreduce(den,dum,diag_flow_npsi)
  ! den is density.
  den=dum/(diag_1d_vol) *1D-6  ! cm-3

  in=int( (0.98*eq_x_psi-diag_1d_pin)/diag_1d_dp ) + 1
  in=min(diag_flow_npsi,max(1,in))
  out=int( (eq_x_psi-diag_1d_pin)/diag_1d_dp   ) + 1
  out=min(diag_flow_npsi,max(1,out))

  maxden=maxval( den(in:out) ) ! finding maxmum density from 0.98*Psi_x ~ Psi_x

!  if(sml_mype==1 ) then
!     print *, neu_mfp, maxden, col_denb_edge,col_denb_edge/maxden,in,out,flow_stat_n
!  endif

  neu_mfp=neu_mfp0*eq_den_v1/maxden


end subroutine change_neu_mfp


!!$subroutine neutral_inside_weight
!!$  use ptl_module
!!$  use neu_module
!!$  real (kind=8) :: psi, psi_interpol
!!$  integer :: i
!!$  neu_old_inside_weight=0D0
!!$
!!$  do i=1, ptl_num
!!$     psi=psi_interpol(ptl_phase(i,1), ptl_phase(i,2),0,0)
!!$     if(psi < eq_x_psi) neu_old_inside_weight=neu_old_inside_weight+ptl_weight(i)
!!$  enddo
!!$end subroutine neutral_inside_weight

subroutine neutrale(grid,psn,sp,old_ptl_num)
  use grid_class
  use psn_class
  use sml_module
  use eq_module, only : eq_x_z, eq_x_psi
  use ptl_module
  use smooth_module
  use random_xgc
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer :: old_ptl_num

  call t_startf("NEUTRALE_SEARCH_INDEX")
  call neutrale_search_index(grid,psn,sp,old_ptl_num)
  call t_stopf("NEUTRALE_SEARCH_INDEX")

  call t_startf("NEUTRALE_SCATTER")
  call neutrale_scatter         ! sub-subroutine
  call t_stopf("NEUTRALE_SCATTER")

  call t_startf("NEUTRALE_MPISUM")
  call neutrale_mpisum(grid,psn,sp)   ! corresponding chargei_gyro_average
  call t_stopf("NEUTRALE_MPISUM")
  psn%eden00_1d=psn%eden00_1d+psn%n2ne00_1d

contains

  subroutine neutrale_scatter
    implicit none
    logical, save :: first = .true.
    real (8),allocatable, save :: n2ne00(:,:)
    real (8):: inv_delta_phi ! for fast calculation -- const
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: nn_beg(sml_nthreads), nn_end(sml_nthreads)
    ! in loop
    integer ::  i,irho, j, node
    real (8) :: phi, wphi(0:1), particle_weight
    real (8) :: wp
    ! for 1D
    integer :: larmor, ip
    real (8) :: x(2), psi, pn
    !
    integer :: jth, iphi
    real (8), external :: psi_interpol

    if (first) then
        allocate(n2ne00(grid%npsi00,sml_nthreads))
       first = .false.
    endif

    inv_delta_phi=1d0/grid%delta_phi

    ! for reduction, divide particles among OpenMP threads
    ! and calculate contributions to nodes for each subset
    ! of particles independently. This will introduce a
    ! roundoff difference in the results for different numbers
    ! of threads.
    call split_indices((sp%num-old_ptl_num), sml_nthreads, i_beg, i_end)

!!$OMP PARALLEL DO &
!!$OMP PRIVATE( ITH, I, IRHO, J, NODE, &
!!$OMP          PHI, WPHI, PARTICLE_WEIGHT, WP, &
!!$OMP          IP, X, PSI, PN )
    do ith=1,sml_nthreads
       call t_startf("NEUTRALE_SCATTER_LOOP")

       n2ne00(:,ith) =0D0

       do i=old_ptl_num+i_beg(ith),old_ptl_num+i_end(ith)
          if(sp%ptl(i)%gid<=0) cycle

          if((sp%tr_save(i)) > 0) then

             x=sp%ptl(i)%ph(1:2)  ! for 1D potential

             ! phi weight
             phi=sp%ptl(i)%ph(3)
             wphi(1) = (phi*inv_delta_phi) - grid%iphi_offset ! larger index weight
             wphi(0) = 1D0 - wphi(1)  ! smaller index weight

             ! particle weight
             if(sml_deltaf) then
                particle_weight=sp%ptl(i)%ph(piw1)*sp%ptl(i)%ct(piw0)
             else
                particle_weight=sp%ptl(i)%ct(piw0) ! for full f simulation only
             endif

             do j=1, 3
                node=grid%nd(j,sp%tr_save(i))
                wp=sp%p_save(j,i)
             enddo

             !1D ---

             psi = psi_interpol(x(1),x(2),0,0)
             pn=(psi-grid%psi00min)/grid%dpsi00
             ip=floor(pn)+1
             if(0 < ip .and. ip < grid%npsi00 .and. (x(2) > eq_x_z .or. psi > eq_x_psi) .and. (.not. sml_00_xz_up .or. x(2)>eq_x_z)) then
                wp=1D0 - ( pn - real(ip-1,8) )  ! recylce of wp -- differenet meaning from above : w of p(1:3) vs w of psi
                n2ne00(ip  ,ith)=n2ne00(ip  ,ith) + particle_weight* wp
                n2ne00(ip+1,ith)=n2ne00(ip+1,ith) + particle_weight*(1D0-wp)
             endif

          else

             !eliminate particle
             call remove_particle(sp,i,-1,ith)

          endif

       enddo

       call t_stopf("NEUTRALE_SCATTER_LOOP")
    enddo  ! end of particle-thread loop

    ! open mp sum

!!$OMP PARALLEL DO &
!!$OMP PRIVATE( JTH, ITH, NODE, IRHO, IPHI )

    do ith=2, sml_nthreads
       n2ne00(:,1)=n2ne00(:,1)+n2ne00(:,ith)
    enddo
    psn%n2ne00_1d=n2ne00(:,1)

  end subroutine neutrale_scatter
end subroutine neutrale

subroutine neutrale_search_index(grid,psn,sp,old_ptl_num)
  use grid_class
  use psn_class
  use ptl_module
  use sml_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  !
  integer :: ith, i,itr, ip
  integer :: i_beg(sml_nthreads), i_end(sml_nthreads)
  real (8) :: phi_mid, x(2), phi, mu, rho, xff(2)
  real (8) :: p(3)
  integer ::   old_ptl_num
  real (kind=8), external :: gyro_radius
  logical, parameter :: USE_SEARCH_TR2 = .true.

  phi_mid=0.5D0*(grid%phimin+grid%phimax)

  call split_indices((sp%num-old_ptl_num), sml_nthreads, i_beg, i_end)

!!$OMP PARALLEL DO &
!!$OMP PRIVATE( ITH, I, PHI,  &
!!$OMP          X, PHI, MU, RHO, XFF,  &
!!$OMP          P, ITR, IP )
  do ith=1,sml_nthreads
     call t_startf("NEUTRALE_SRCHLOOP") 
     do i=old_ptl_num+i_beg(ith),old_ptl_num+i_end(ith)
        if(sp%ptl(i)%gid>0)  then

           ! get proper toroidal angle index and weight
           x=sp%ptl(i)%ph(1:2)
           phi=sp%ptl(i)%ph(3)
!           mu=sp%ptl(i)%ct(pim)
!           rho=gyro_radius(x,mu)  !gyro radius
!           sp%rhoi(i)=rho

             ! get field following posision at 1/2 angle
           call field_following_pos2(x,phi,phi_mid,xff)

           ! find triangle
           if (USE_SEARCH_TR2) then
              call search_tr2(grid,xff,itr,p)
           else
              call search_tr(grid,xff,itr,p)
           endif

           sp%tr_save(i)=itr
           sp%p_save(:,i)=p
           ! compare performance
           !do ip=1,3
           !   sp%p_save(ip,i)=p(ip)
           !enddo
        endif
     enddo
     call t_stopf("NEUTRALE_SRCHLOOP") 
  enddo
end subroutine neutrale_search_index

subroutine neutrale_mpisum(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  include 'mpif.h'
  !
!  integer :: nn, ierr
!  integer :: irho, iphi, ipe
!  integer :: idest, isendtag, isource, irecvtag
!  integer, dimension(MPI_STATUS_SIZE) :: istatus
!  real (8) :: tmp(grid%nnode,0:1)
  real (8) :: dum00(grid%npsi00)
!  real (8) :: inv_nphi_total
  !
!  character (len=30) :: filename
!  integer :: i


  ! 00 mode
  call my_mpi_allreduce(psn%n2ne00_1d,dum00,grid%npsi00)
  psn%n2ne00_1d=dum00/psn%vol00

end subroutine neutrale_mpisum
