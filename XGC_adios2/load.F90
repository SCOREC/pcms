! Initialize Phase variabls of ions.
!!>For single particle simulation, call load_single
!!<For V-space hole simulation, call load_special3
subroutine load(grid,psn,spall)
  use ptl_module
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  implicit none
  !use grid_class
  !use psn_class
  !use eq_module
  !use random_xgc
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i, isp
  integer (kind=8) :: ioff, idum1,idum2
  real (kind=8) :: r,z,phi,rho,mu,b,psi
  real (kind=8) :: psi_interpol, b_interpol,ni,ti_ev
  real (kind=8) :: ni_c, ti_ev_c, psi_c, en_ev,br,bz,bphi,ti_ev_vertual
#ifdef FORCE_ZEROCHARGE
  real (kind=8) :: totpert,corpert
  integer (kind=8) :: nptot,npcor,ierr
#endif

#ifdef FORCE_ZEROCHARGE_NEW
  real (kind=8), dimension(:), allocatable :: totpert_binpsi,corpert_binpsi
  integer (kind=8), dimension(:), allocatable :: totnp_binpsi, cornp_binpsi, binpsi 
  integer (kind=8) :: nbin_psi,ibinpsi
  real (kind=8) :: dpsi,w
#endif
  ! loading for single particle simulation
  if(sml_special==1) then
     call load_single(spall)
     return
  endif
  
  ! restart from previous runs -----------------------------
  if(sml_restart)then
#if defined(ADIOS) 
     call restart_read(grid,psn,spall)  ! sml_restart is set to zero if no restart file (timestep.dat)
    if(sml_restart) then
       return
    endif
#else
     stop '** error ** Compile with adios or binpack support for restart!!'
#endif
  endif
    
  if(.not. sml_restart) then
     sml_istep=0
     sml_gstep=0
     sml_time=0.D0
     sml_run_count=1
  endif
  !---------------------------------------------------------

  ! First, set particle position (r,z,phi) uniformly in simulation region
  call uniform_space_dist(grid,psn,spall)
!  if(.not. sml_deltaf .and. sml_electron_on .and. spall(1)%num==spall(0)%num) then
!     !debug
!     spall(0)%phase(1:3,1:spall(1)%num) = spall(1)%phase(1:3,1:spall(1)%num)
!     spall(1)%phase(6:9,1:spall(1)%num) = spall(1)%phase(6:9,1:spall(1)%num)
!  endif

#ifdef MAXWELL_V_DIST_ORG
  call maxwell_v_dist(grid,spall)
#else
  call maxwell_v_dist2(grid,spall)
#endif


  do isp=ptl_isp, ptl_nsp

!     ioff = sml_mype*sp%num
     idum1=sml_mype
     idum2=spall(isp)%num 
     ioff = idum1 * idum2
     do i=1,spall(isp)%num
        spall(isp)%ptl(i)%gid=ioff+i
     enddo

!     sp%maxgid=sml_totalpe*sp%num
     idum1=sml_totalpe
     idum2=spall(isp)%num
     spall(isp)%maxgid=idum1*idum2
  enddo

#ifdef FORCE_ZEROCHARGE
  !Julien force the total perturbation charge of each species to be zero
  !delta f only
  do isp=ptl_isp, ptl_nsp
    if(isp==1)then !Ion
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: Cancelling n=0 m=0 <>_r  mode'
    endif
    totpert=0D0
    corpert=0D0
    nptot=0
    npcor=0
    do i=1,spall(isp)%num
      totpert=totpert+spall(isp)%ptl(i)%ct(piw0)*spall(isp)%ptl(i)%ph(piw1)
      npcor=npcor+1
    enddo
    call MPI_ALLREDUCE(totpert,corpert,1,MPI_REAL8   ,MPI_SUM,sml_comm,ierr)
    call MPI_ALLREDUCE(npcor  ,nptot  ,1,MPI_INTEGER8,MPI_SUM,sml_comm,ierr)
    totpert=corpert
    corpert = corpert / REAL(nptot,KIND(corpert))
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: ',totpert,corpert,nptot
    endif
    totpert=0D0
    do i=1,spall(isp)%num
      spall(isp)%ptl(i)%ph(piw1)=spall(isp)%ptl(i)%ph(piw1)-corpert/spall(isp)%ptl(i)%ct(piw0)
      totpert=totpert+spall(isp)%ptl(i)%ct(piw0)*spall(isp)%ptl(i)%ph(piw1)
    enddo
    call MPI_ALLREDUCE(totpert,corpert,1,MPI_REAL8 ,MPI_SUM,sml_comm,ierr)
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: 0==',corpert
    endif
    endif
  enddo
#endif

#ifdef FORCE_ZEROCHARGE_NEW
  !Julien force the total perturbation charge of each species to be zero
  !delta f only
  nbin_psi=sml_00_npsi !300
  dpsi    =eq_x_psi/dble(nbin_psi)
  allocate(totnp_binpsi(1:nbin_psi))
  allocate(totpert_binpsi(1:nbin_psi))
  allocate(corpert_binpsi(1:nbin_psi))
  allocate(cornp_binpsi(1:nbin_psi))
  do isp=ptl_isp, ptl_nsp
    if(isp==1)then !Ion
      allocate(binpsi(1:spall(isp)%num))
      if(sml_mype==0)then
        print *,'FORCE_ZEROCHARGE00: Cancelling n=0 m=0 mode'
      endif
      totpert_binpsi=0D0
      corpert_binpsi=0D0
      totnp_binpsi=0
      cornp_binpsi=0
      do i=1,spall(isp)%num
        r=spall(isp)%ptl(i)%ph(pir)
        z=spall(isp)%ptl(i)%ph(piz)
        psi=psi_interpol(r,z,0,0)
        ibinpsi=max(1,min(floor(psi/dpsi),nbin_psi))
        binpsi(i)=ibinpsi
        w=spall(isp)%ptl(i)%ct(piw0)*spall(isp)%ptl(i)%ph(piw1)
        totpert_binpsi(ibinpsi)=totpert_binpsi(ibinpsi)+w
        totnp_binpsi(ibinpsi)=totnp_binpsi(ibinpsi)+1
      enddo
!=================      !if(sml_mype==0) print *,'ALL REDUCE',sml_mype   
      call MPI_ALLREDUCE(totpert_binpsi,corpert_binpsi,nbin_psi,MPI_REAL8   ,MPI_SUM,sml_comm,ierr)
      call MPI_ALLREDUCE(totnp_binpsi  ,cornp_binpsi  ,nbin_psi,MPI_INTEGER8,MPI_SUM,sml_comm,ierr)
      do i=1,nbin_psi
        if(cornp_binpsi(i)>0)then
          corpert_binpsi(i) = corpert / REAL(npcor,KIND=8)
        else
          corpert_binpsi(i) = 0D0
        endif
      enddo
!call mpi_barrier(sml_comm)
      if(sml_mype==0)then
        print *,'FORCE_ZEROCHARGE00: '!,totpert,corpert,nptot
      endif
      do i=1,spall(isp)%num
        if(cornp_binpsi(binpsi(i))>0)then
          w=corpert_binpsi(binpsi(i))/spall(isp)%ptl(i)%ct(piw0)
          spall(isp)%ptl(i)%ph(piw1)=spall(isp)%ptl(i)%ph(piw1)-w
        endif
      enddo
!call mpi_barrier(sml_comm)
      if(sml_mype==0)then
        print *,'FORCE_ZEROCHARGE00: OUT'
      endif
    endif!EOF isp==1
  enddo
!call mpi_barrier(sml_comm)
  if(allocated(binpsi)) deallocate(binpsi)
  if(allocated(totpert_binpsi)) deallocate(totpert_binpsi)
  if(allocated(totnp_binpsi)) deallocate(totnp_binpsi)
  if(allocated(corpert_binpsi)) deallocate(corpert_binpsi)
  if(allocated(cornp_binpsi)) deallocate(cornp_binpsi)
#endif

end subroutine load

#ifdef REFACTOR_ZER_CHARGE
subroutine force_zero_charge(spall,choice)
  use ptl_module
  use sml_module
!  use grid_class
!  use psn_class
  implicit none
  include 'mpif.h'
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i, isp
  real (kind=8) :: totpert,corpert
  integer (kind=8) :: nptot,npcor,ierr,choice

  if(choice==0)then
  !Julien force the total perturbation charge of each species to be zero
  !delta f only
  do isp=ptl_isp, ptl_nsp
    if(isp==1)then !Ion
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: Cancelling n=0 m=0 <>_r  mode'
    endif
    totpert=0D0
    corpert=0D0
    nptot=0
    npcor=0
    do i=1,spall(isp)%num
      totpert=totpert+spall(isp)%ptl(i)%ct(piw0)*spall(isp)%ptl(i)%ph(piw1)
      npcor=npcor+1
    enddo
    call MPI_ALLREDUCE(totpert,corpert,1,MPI_REAL8   ,MPI_SUM,sml_comm,ierr)
    call MPI_ALLREDUCE(npcor  ,nptot  ,1,MPI_INTEGER8,MPI_SUM,sml_comm,ierr)
    totpert=corpert
    corpert = corpert / REAL(nptot,KIND(corpert))
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: ',totpert,corpert,nptot
    endif
    totpert=0D0
    do i=1,spall(isp)%num
      spall(isp)%ptl(i)%ph(piw1)=spall(isp)%ptl(i)%ph(piw1)-corpert/spall(isp)%ptl(i)%ct(piw0)
      totpert=totpert+spall(isp)%ptl(i)%ct(piw0)*spall(isp)%ptl(i)%ph(piw1)
    enddo
    call MPI_ALLREDUCE(totpert,corpert,1,MPI_REAL8 ,MPI_SUM,sml_comm,ierr)
    if(sml_mype==0)then
      print *,'FORCE_ZEROCHARGE: 0==',corpert
    endif
    endif
  enddo
  elseif(choice==2)then

  elseif(choice==3)then

  endif
 
end subroutine force_zero_charge
#endif

! obtain maxwell distribution function with w0 weight
! actual distribution function is not maxwellian
subroutine maxwell_v_dist2(grid,spall)
  use grid_class
  use sml_module
  use ptl_module
  use eq_module
  use random_xgc
  use initial_perturbation
  implicit none
  type(grid_type), intent(in)::grid
  type(species_type) :: spall(0:ptl_nsp_max)
  real (kind=8) :: r,z,phi,b,psi,n,t_ev,t_ev_vertual,rho,mu,n_c,t_ev_c,bphi,br,bz,en_ev,psi_c,maxe
  real (kind=8) :: rho_shift, up
  integer :: i,isp
  real (kind=8),external :: b_interpol, psi_interpol!, init_ipara_flow, init_epara_flow
  real (kind=8) :: marker_den, w0_adjust, theta
#ifdef INIT_GENE_PERT
  !Julien, April 2017
  real (kind=8) :: dN, phig
#endif
  ! set initial density particle weight -- assuming maxwellian distribution function
!  do isp=ptl_isp, ptl_nsp
!     if(isp==1) then    ! main ion
!        marker_den=sml_marker_den
!     elseif(isp==0) then ! electron
!        marker_den=sml_marker_den*real(spall(0)%num)/real(spall(1)%num)
!     endif

!     do i=1, spall(isp)%num
        !setting particle weight          
!        spall(isp)%ptl(i)%ct(piw0)=n/marker_den !sml_marker_den=(ptl number per unit volume) -- background profile, 2002/05/30
!     enddo
!  enddo
#ifdef INIT_GENE_PERT
  call ip_initialize()
#endif

  ! obtain rho and mu
  do isp=ptl_isp, ptl_nsp
     if(isp==1) then    ! main ion
        marker_den=sml_marker_den
     elseif(isp==0) then ! electron
        marker_den=sml_marker_den*real(spall(0)%num)/real(spall(1)%num)
     endif

     do i=1, spall(isp)%num
        !retrive position variable from the uniform_space_dist
        r=spall(isp)%ptl(i)%ph(1)
        z=spall(isp)%ptl(i)%ph(2)
        phi=spall(isp)%ptl(i)%ph(3)
        b=b_interpol(r,z,phi)
        psi=psi_interpol(r,z,0,0)
     
        
        !get density and temp
        if(isp==1) then
           !n=init_den_wz(psi,z)  !get density
           ! t_ev=init_tempi_ev_wz(psi,z)
           n=eq_ftn(psi,r,z,eq_den)
           t_ev=eq_ftn(psi,r,z,eq_tempi)
        else
           !n=init_den_wz(psi,z)
           !t_ev=tempe_ev_wz(psi,z)
           n=eq_ftn(psi,r,z,eq_den)
           t_ev=eq_ftn(psi,r,z,eq_tempe)
        endif

        !setting particle weight
        spall(isp)%ptl(i)%ct(piw0)=n/marker_den !sml_marker_den=(ptl number per unit volume) -- background profile, 2002/05/30


        ! Parallel velocity offset - initial flow
        if(sml_initial_flow) then
           if(ptl_deltaf_sp(isp) .and. sml_mype==0) then
              print *, 'sml_initial_flow is not implimented with delta-f, yet'
           endif
           if(isp==1) then ! ion
              up=eq_ftn(psi,r,z,eq_flowi) !init_ipara_flow(psi,z)
           else
              up=eq_ftn(psi,r,z,eq_flowe) !init_epara_flow(psi,z)
              if(sml_mype==0) then
                 print *, 'sml_initial_flow is not implimented with electron, yet'
              endif
           endif
           up=eq_axis_r/r*up ! constant angular velocity assuming toroidal flow
           !rho=rho + up/(ptl_c_m(isp)*b)
        else
           up=0D0
        endif
        

        !**************************************************************
        !  actual routine that calculate v and w0 adjustments
        !
        if(.not. sml_flat_marker) then
           call load_v_single(t_ev, up, b, ptl_mass(isp), ptl_charge(isp), rho, mu, w0_adjust)
        else
           call load_flat_v_single(t_ev, up, b, ptl_mass(isp), ptl_charge(isp), rho, mu, w0_adjust)
        endif
        !
        !**************************************************************


        spall(isp)%ptl(i)%ph(pirho)=rho                          ! set rho and mu to phase variable
        spall(isp)%ptl(i)%ct(pim)=mu


        rho_shift=rho - up/(ptl_c_m(isp)*b)
        en_ev=(mu*B + 0.5D0*(ptl_charge(isp)*rho_shift*B)**2/ptl_mass(isp))*sml_j2ev
        spall(isp)%ptl(i)%ct(pif0)=n*sqrt(1D0/t_ev**3)*exp(-en_ev/t_ev)

        spall(isp)%ptl(i)%ct(piw0)=spall(isp)%ptl(i)%ct(piw0)*w0_adjust


        ! applying initial noise
        if(ptl_deltaf_sp(isp)) then
#ifndef INIT_GENE_PERT           
          theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
          if(z<eq_axis_z) then
           theta= - theta
          endif

          
          spall(isp)%ptl(i)%ph(piw1:piw2) = exp(-0.5* ((psi/0.2059-0.383)/0.06)**2  - 0.5* (theta*4/sml_pi)**2) &
            * sml_initial_deltaf_noise*2D0*(ranx()-0.5D0) 

           ! for more exact calculation, adjust ct(pif0) accordingly  ( when sml_initial_deltaf_noise > 1D-4 )
           ! ct(pif0) should be ct(pif0) / (1-w1) or ct(pif0) * (1-w1)
#else
!         phig = modulo(phi+sml_pi,2D0*sml_pi)-sml_pi
         phig = phi
         call ip_eval(dN,r,z,phig)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#warning 'code is totally broken in load.F90:384'
!         dN = 1d0
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if(sml_initial_deltaf_noise.gt.0D0)then
           spall(isp)%ptl(i)%ph(piw1:piw2) = dN * (pertw + sml_initial_deltaf_noise*2D0*(ranx()-0.5D0))
         else
            spall(isp)%ptl(i)%ph(piw1) = dN * pertw
            spall(isp)%ptl(i)%ph(piw2) = 0.0
         endif
!         if(sml_initial_deltaf_noise>0D0)then
!           spall(isp)%ptl(i)%ph(piw1:piw2) = spall(isp)%ptl(i)%ph(piw1:piw2) * (pertw + sml_initial_deltaf_noise*2D0*(ranx()-0.5D0))
!         endif
#endif

#ifdef GAM_TEST
           spall(isp)%ptl(i)%ph(piw1:piw2)=gam_weight(r,z,psi)
#endif
        else
           spall(isp)%ptl(i)%ph(piw1:piw2)=0D0
        endif

     enddo
  enddo

#ifdef INIT_GENE_PERT
  call ip_destroy()
#endif

end subroutine maxwell_v_dist2


!  flat marker distribution function
! The distribution funtion has a form of
! g = C    when v <= v_a
!   = C exp(-(v-v_a)/w)   when v > v_a
! The PDF (integrated) is
! G(v) = C v     when v <= v_a
!      = C[ v_a + w (1-exp(-(v-v_a)/w))]   when v > v_a
!      where C = 1/[ v_a + w ( 1- exp(-(v_c - v_a)/w))]
!
! The inverse of G is
! G_inv(A) = A/C     when A < C v_a
!          = v_a - w ln [1 + (v_a - A/C)/w]   when A > C v_a
!
subroutine load_flat_v_single(t_ev, up, b, mass, charge, rho, mu, w0_adjust)
  use sml_module
  use random_xgc
  implicit none
  real (8) , intent(in) :: t_ev, up, b, mass, charge
  real (8) , intent(out) :: rho, mu, w0_adjust
  !
  real (8) :: temp, vth, maxwell_norm,  vshift
  real (8) :: va, w, vc, v, g
  ! C :: normalized constant
  ! va : start of exp decay
  ! w : width of exp decay
  ! vc : cut-off v   - no particle after v_c

  temp=t_ev*sml_ev2j
  vth=sqrt(temp/mass)
  maxwell_norm=sqrt(mass/sml_2pi/temp)

  ! parallel velocity distribution
  va=sml_flat_marker_decay_start1 *vth
  vc=sml_flat_marker_cutoff1      *vth
  w =sml_flat_marker_width1       *vth

  call get_v_and_dist

  v=v*sign(1D0, ranx()-0.5D0)
  vshift= v + up ! shift of velocity
  rho=vshift/b*mass/charge !*mass/charge

  
  w0_adjust = maxwell_norm*dexp(-0.5D0*mass*v*v/temp) / (g*0.5D0)   ! actual g is half due to -v and v direction
  

  ! perp velocity distribution
  
  va=sml_flat_marker_decay_start2 *vth
  vc=sml_flat_marker_cutoff2      *vth
  w =sml_flat_marker_width2       *vth
  
  call get_v_and_dist
  
  mu=0.5D0*mass*v*v/b
  
  w0_adjust = w0_adjust * &
       mass/temp * dexp( -0.5D0*mass*v*v/temp) * v / g

contains
  subroutine get_v_and_dist  ! input:: v_a, v_c, w , output:: v, g
    implicit none
    real (8) :: c   ! c is constant for each direction - can be stored as global variable
    real (8) :: A

    c=1D0/(va + w*(1D0 - dexp(-(vc-va)/w)))
    
    A=ranx()
    if(A < c*va) then
       v=A/c
       g=c
    else
       v=va - w*dlog(1D0 + (va-A/c)/w)
       g=c*dexp(-(v-va)/w)
    endif

  end subroutine get_v_and_dist

end subroutine load_flat_v_single



subroutine load_v_single(t_ev, up, b, mass, charge, rho, mu, w0_adjust)
  use sml_module
  use random_xgc
  implicit none
  real (8) , intent(in) :: t_ev, up, b, mass, charge
  real (8) , intent(out) :: rho, mu, w0_adjust
  !
  real (8) :: t_ev_vertual, maxe, t, zmax, zdum, v, vshift, t1, t2, en_ev
 
  
  ! Parallel Velocity distribution ***************************************************************
  t_ev_vertual=t_ev*sml_marker_temp_factor  !set marker particle temperature for parallel
  
  !adjust maker particle temperature
  if(t_ev_vertual < sml_marker_min_temp) then
     t_ev_vertual = sml_marker_min_temp
  endif

  ! cut-off energy
  maxe= sml_load_maxe*t_ev/t_ev_vertual
  
  ! temperature in Joule 
  t= t_ev_vertual *sml_ev2j
  
  ! v_parap from random number
  zmax=1D0 - dexp(-maxe)
  zdum=zmax*ranx()
  v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*t)
  v= v*cos(sml_pi*ranx())
  vshift= v + up ! shift of velocity
  rho=vshift/b*mass/charge !*mass/charge
  
  en_ev = 0.5D0*mass*v*v*sml_j2ev

  ! w0_adjust
  w0_adjust =    exp( - en_ev*(1D0/t_ev - 1D0/t_ev_vertual) ) &
                * (t_ev_vertual/t_ev*(1D0-exp(-maxe))/(1D0-exp(-sml_load_maxe)))**0.5D0
  

  ! Perp velocity distribution *******************************************************************
  t_ev_vertual=t_ev*sml_marker_temp_factor2  ! for low mu particle population

  !adjust maker particle temperature
  if(t_ev_vertual < sml_marker_min_temp) then
     t_ev_vertual = sml_marker_min_temp
  endif

  ! cut-off energy
  maxe= sml_load_maxe*t_ev/t_ev_vertual
  
  ! temperature in Joule 
  t1=t_ev_vertual*sml_ev2j
  t2=t_ev*sml_marker_temp_factor3*sml_ev2j
  
  
  if(ranx()>sml_low_mu_fill_population) then  ! main distribution
     t=t1

     zmax=1D0 - dexp(-maxe)
     zdum=zmax*ranx()
     v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*t)
     mu=0.5D0*mass*v**2/b
     en_ev = mu*b*sml_j2ev
     
     ! w0_adjust
     w0_adjust = w0_adjust* exp( - en_ev*(1D0/t_ev - 1D0/t_ev_vertual) ) &
          * (t_ev_vertual/t_ev*(1D0-exp(-maxe))/(1D0-exp(-sml_load_maxe)))!**1D0
     
  else
     t=t2
     
     zmax=1D0 - dexp(-sml_load_maxe)
     zdum=zmax*ranx()
     v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*t)
     mu=0.5D0*mass*v**2/b
     en_ev = mu*b*sml_j2ev
     
     ! w0_adjust
     w0_adjust = w0_adjust* exp( - en_ev*(1D0/t_ev - 1D0/(t*sml_j2ev)) ) &
          * ((t*sml_j2ev)/t_ev*(1D0-exp(-sml_load_maxe))/(1D0-exp(-sml_load_maxe)))!**1D0
  endif
  w0_adjust=1.0
     
end subroutine load_v_single



subroutine maxwell_v_dist(grid,spall)
  use grid_class
  use sml_module
  use ptl_module
  use eq_module
  use random_xgc
  implicit none
  type(grid_type), intent(in)::grid
  type(species_type) :: spall(0:ptl_nsp_max)
  real (kind=8) :: r,z,phi,b,psi,n,t_ev,t_ev_vertual,rho,mu,n_c,t_ev_c,bphi,br,bz,en_ev,psi_c,maxe
  integer :: i,isp
  real (kind=8),external :: b_interpol, psi_interpol!, init_ipara_flow, init_epara_flow
  real (kind=8) :: marker_den
  real (8) :: up
!  integer, pointer :: pnum
!  real (kind=8), pointer :: phase(:,:)
#ifdef GAM_TEST
  real (kind=8), external :: gam_weight
#endif
  real (8) :: enp


  do isp=ptl_isp, ptl_nsp
     if(isp==1) then    ! main ion
        marker_den=sml_marker_den
     elseif(isp==0) then ! electron
        marker_den=sml_marker_den*real(spall(0)%num)/real(spall(1)%num)
     endif

     do i=1, spall(isp)%num
        !retrive position variable from the uniform_space_dist
        r=spall(isp)%ptl(i)%ph(1)
        z=spall(isp)%ptl(i)%ph(2)
        phi=spall(isp)%ptl(i)%ph(3)
        b=b_interpol(r,z,phi)
        psi=psi_interpol(r,z,0,0)
     
     
        !get density and temp
        if(isp==1) then
           !n=init_den_wz(psi,z)  !get density
           ! t_ev=init_tempi_ev_wz(psi,z)
           n=eq_ftn(psi,r,z,eq_den)
           t_ev=eq_ftn(psi,r,z,eq_tempi)
        else
           !n=init_den_wz(psi,z)
           !t_ev=tempe_ev_wz(psi,z)
           n=eq_ftn(psi,r,z,eq_den)
           t_ev=eq_ftn(psi,r,z,eq_tempe)
        endif

        !setting particle weight          
        spall(isp)%ptl(i)%ct(piw0)=n/marker_den !sml_marker_den=(ptl number per unit volume) -- background profile, 2002/05/30
        ! set very small or zero value for delta-f weight
        if(ptl_deltaf_sp(isp)) then
!           if(sml_deltaf_f0_mode==1) then
!              spall(isp)%ptl(i)%ph(piw1:piw2)= (n-f0_den2(grid,r,z,psi))/n  
!           else
              spall(isp)%ptl(i)%ph(piw1:piw2)=0D0
!           endif
           if(abs(spall(isp)%ptl(i)%ph(piw1)) <= sml_initial_deltaf_noise ) then
              spall(isp)%ptl(i)%ph(piw1:piw2) = sml_initial_deltaf_noise*2D0*(ranx()-0.5D0)*spall(isp)%ptl(i)%ph(1)/eq_axis_r ! ballooning mode shape
           endif
#ifdef GAM_TEST
           spall(isp)%ptl(i)%ph(piw1:piw2)=gam_weight(r,z,psi)
#endif
        else
           spall(isp)%ptl(i)%ph(piw1:piw2)=0D0
        endif
     
        ! Velocity distribution
        t_ev_vertual=t_ev*sml_marker_temp_factor  !set marker particle temperature
        !maxe=12D0/sml_marker_temp_factor

        !adjust maker particle temperature
        if(t_ev_vertual < sml_marker_min_temp) then
           t_ev_vertual = sml_marker_min_temp
        endif
        maxe= sml_load_maxe*t_ev/t_ev_vertual

        call maxwell_dist(t_ev_vertual,b,rho,mu,ptl_mass(isp),ptl_charge(isp),maxe)   !get rho and mu from Monte-Carlo method
        
        ! Parallel velocity offset - initial flow
        if(sml_initial_flow) then
           if(ptl_deltaf_sp(isp) .and. sml_mype==0) then
              print *, 'sml_initial_flow is not implimented with delta-f, yet'
           endif
           if(isp==1) then ! ion
              up=eq_ftn(psi,r,z,eq_flowi) !init_ipara_flow(psi,z)
           else
              up=eq_ftn(psi,r,z,eq_flowe) !init_epara_flow(psi,z)
              if(sml_mype==0) then
                 print *, 'sml_initial_flow is not implimented with electron, yet'
              endif
           endif
           up=eq_axis_r/r*up ! constant angular velocity assuming toroidal flow
           rho=rho + up/(ptl_c_m(isp)*b)
        endif
        spall(isp)%ptl(i)%ph(pirho)=rho                          ! set rho and mu to phase variable
        spall(isp)%ptl(i)%ct(pim)=mu 

        if(sml_marker_temp_factor/=1D0) then
           spall(isp)%ptl(i)%ph(pirho)=rho                          ! set rho and mu to phase variable
           spall(isp)%ptl(i)%ct(pim)=mu
           en_ev=(ptl_c2_2m(isp)*(rho*b)**2 + mu*B)*sml_j2ev
           spall(isp)%ptl(i)%ct(piw0)=spall(isp)%ptl(i)%ct(piw0)*exp( - en_ev*(1D0/t_ev - 1D0/t_ev_vertual) ) &
                * (t_ev_vertual/t_ev*(1D0-exp(-maxe))/(1D0-exp(-sml_load_maxe)))**1.5D0  ! modification due to energy cut-off are added 1/3/2013 -- (1 - exp(-E_max) factor
        endif
        
        ! canonical maxwellian distribution
        ! Adjust weight.
        ! Not good for even distribution - Algorithm need to be changed later

        !set initial f0  - ignore initial flow for quick coding - add later
!        en_ev=(ptl_c2_2m(isp)*(rho*b)**2 + mu*B)*sml_j2ev        
        en_ev=(mu*B + 0.5D0*(ptl_charge(isp)*rho*B)**2/ptl_mass(isp))*sml_j2ev
!rh #ifndef MU_LINEAR
!rh         enp=mu*eq_axis_b*sml_j2ev/t_ev
!rh         ! Do not use sqrt(enp), it can cause problems in weight update
!rh         !spall(isp)%ptl(i)%ct(pif0)=n/t_ev * exp(-en_ev/t_ev) * 2D0*sqrt(enp)
!rh         ! Use this instead -->
!rh         spall(isp)%ptl(i)%ct(pif0)=n*sqrt(1D0/t_ev**3)*exp(-en_ev/t_ev)
!#else
        spall(isp)%ptl(i)%ct(pif0)=n*sqrt(1D0/t_ev**3)*exp(-en_ev/t_ev)
!rh #endif
     enddo

  enddo
     
end subroutine maxwell_v_dist


!!  Distribute particle uniformly in space. (inpsi < psi < outpsi)
!!> The r,z,phi,psi of phase variable are set
!!< 
subroutine uniform_space_dist(grid,psn,spall)
  use grid_class
  use psn_class
  use ptl_module
  use sml_module
  use eq_module
  use random_xgc
  implicit none
  include 'mpif.h'
  type(grid_type), intent(in) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: isp
  integer :: valid, ierr  !! number of particle that generated inside the region inpsi<psi<outpsi
  integer (8) :: total  !! number of total particle generated 
  real (kind=8) :: rdim, zdim   !! R,Z dimension of whole simulation region
  real (kind=8) :: roffset, zoffset !! offset (start position) of simulation region
  real (kind=8) :: r,z,psi,phi
  real (kind=8) , external :: psi_interpol
  real (kind=8) :: x(2),p(3), valid2, dum, xff(2), phi_mid
  integer :: itr

! simulation boundary is imposed 2001/01/24
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z
  
  do isp=ptl_isp, ptl_nsp
     valid=0
     valid2=0
     total=0

     ! generate particle until # of valid particles become ptl_num 
     do while(valid<spall(isp)%num)
        !generate r,z in simulation region

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#warning 'The code is totally broken in uniform_space_dist' 
        r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) !2002/05/27
        z=zdim*ranx() + zoffset
!! test for distribution
!        r=2d0
!        z=0d0
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psi=psi_interpol(r,z,0,0)
        total=total+1

        phi=(ranx()+ grid%iphi_offset)*grid%delta_phi   ! set toroidal angle
        if(sml_2pi_wedge_n < phi) phi=sml_2pi_wedge_n   ! just for in case
        if(phi<0D0) phi=0D0             ! just for in case

        !check psi validity
        if(sml_inpsi < psi .AND. psi < sml_outpsi) then
           valid2=valid2+1D0


           ! check inside wall validity
           x(1)=r
           x(2)=z
           phi_mid=(floor(phi/grid%delta_phi) + 0.5D0) * grid%delta_phi
  
           ! get field following posision at 1/2 angle
           call field_following_pos2(x,phi,phi_mid,xff)

           ! find position of previous time step
           !if (USE_SEARCH_TR2) then
              call search_tr2(grid,xff,itr,p)
           !else
           !   call search_tr(grid,xff,itr,p)
           !endif
                      
           if(itr>0)then
              valid=valid+1
              !set phase variables to global storage
              spall(isp)%ptl(valid)%ph(1:3)=(/ r, z, phi /)
           end if
        endif
     enddo
  enddo
 
  !get sml_marker_den correctly
  call mpi_allreduce(valid2,dum,1, mpi_real8, mpi_sum, sml_comm,ierr)
  sml_marker_den=dum/real(sml_totalpe)/real(sml_monte_num)*sml_marker_den

  ! check validity of toroidal anagle range
  call shift_check(grid,spall(1))
  if(sml_electron_on) call shift_check(grid,spall(0))
  

end subroutine uniform_space_dist



subroutine maxwell_dist(ti_ev,b,rho,mu,mass,charge,maxe)
  use sml_module
  use random_xgc
  implicit none
  real (kind=8) , intent(in) :: b,ti_ev,mass,charge,maxe
  real (kind=8) , intent(out) :: rho,mu
  real (kind=8) :: theta,r,z,pitch,phi,zdum,zmax,v,ti

  ti=ti_ev*sml_ev2j
  
  ! perp velocity
  zmax=1.d0 - dexp(-maxe)
  zdum=zmax*ranx()
!  v= sqrt(-2.d0*dlog(1.d0-zdum)*ti/mass)
  v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*ti)
  mu=0.5D0*mass*v**2/b   ! mu is indep to mass and charge for given temp
  ! parallel velocity
  zdum=zmax*ranx()
!  v= sqrt(-2.d0*dlog(1.d0-zdum)*ti/mass)
  v= sqrt(-2.d0/mass*dlog(1.d0-zdum)*ti)

  v= v*cos(sml_pi*ranx())
  rho=v/b*mass/charge !*mass/charge
  
endsubroutine


!!************************************************************
!!>Calculating the voulme of shells
!!
!! first created 2002/05/27
!! adopted from get_marker_den -> renaming and modifying
!! r dependant loading
!!<***********************************************************
subroutine get_volume(grid,psn,spall)
  use grid_class
  use psn_class
  use diag_module
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  use random_xgc
  use perf_monitor
  implicit none
  type(grid_type) ::grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  integer, parameter :: npos_1d=2
  integer ::  valid    !! number of valid particle 
  integer (8) ::  total    !! number of total particle generated
  integer :: j,ierror, sp, i,l
  real (8) :: rdim, zdim  !! R,Z dimension of whole simulation region
  real (8) :: roffset, zoffset !! offset (start position) of simulation region
  real (8) :: mc_den
  real (8) :: tvolume, vvolume !! Total volume and valid volume
  real (8) :: r,z,psi,theta,phi,dum1,dum2
  !
  real (8) :: dum_1d(diag_1d_npsi)  ! for diag_1d
  real (8) :: dum_neu(diag_neu_npsi) 
  real (8) :: dum_00(grid%npsi00)   ! for 00 charge density
  real (8) :: dum_c(col_2_mtheta, 0:col_2_m)   ! for conserving collision
  real (8) :: dum_vb(col_2_mtheta,0:col_2_m)   ! for bg update
  !
  real (8) , external :: psi_interpol
  real (8) :: aa,bb   !! Linear interpolation weight
  ! node_vol_ff
  logical, parameter :: USE_SEARCH_TR2=.true.
  real (8), allocatable :: dum_ff(:,:)
  real (8) :: x(2), xff(2), phi_dest, p(3), wphi
  integer :: itr, node, ml(1)


! simulation boundary is imposed 2001/01/24
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  valid=0
  total=0
  
  ! initialize volume
  diag_1d_vol=1D-99
  diag_neu_vol=1D-99
  psn%vol00=1D-99
  grid%node_vol_ff=1D-99
  if(sml_f0_grid) grid%node_vol_nearest=1D-99
  if(col_mode/=0) then
      if(col_varying_bg==1) col_vb_vol=1D-99
      if(col_mode==2) col_2_vol=1D-99
  endif
  
  phi_dest = 0.5D0*grid%delta_phi
  ! mem allocation
  allocate(dum_ff(grid%nnode,2))

!  print *, 'get_vol'
  do while(valid<sml_monte_num)
     if (sml_cylindrical) then
       ! Cylindrical limit with periodic boundary conditions
       r=rdim*ranx() + roffset
     else
       ! Toroidal geometry
       r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) )
     endif
     z=zdim*ranx() + zoffset
     phi=ranx()*grid%delta_phi
     psi=psi_interpol(r,z,0,0)
     total=total+1
     
     ! marker density
     if(sml_inpsi < psi .AND. psi < sml_outpsi) then
        valid=valid+1
     endif
     
     if(grid%psi00min < psi .and. psi < grid%psi00max .and. is_rgn12(r,z,psi) ) then
        j=int((psi-grid%psi00min)/grid%dpsi00)+1
        j=min(grid%npsi00-1,max(1,j))
        bb=(psi-grid%psi00min)/grid%dpsi00 + 1 -j
        aa= 1D0-bb
        psn%vol00(j)  =psn%vol00(j  )+1D0*aa
        psn%vol00(j+1)=psn%vol00(j+1)+1D0*bb
     endif
     
     if(diag_1d_pin < psi .and. psi < diag_1d_pout .and. (z>eq_x_z) .and. z<eq_x2_z .and. is_rgn12(r,z,psi)) then ! diagnosis excludes region 2 z<eq_x_z
        j=int((psi-diag_1d_pin)/diag_1d_dp)+1
        j=min(diag_1d_npsi-1,max(1,j))
        bb=(psi-diag_1d_pin)/diag_1d_dp + 1 -j
        aa=1D0-bb
        diag_1d_vol(j)=diag_1d_vol(j)+1D0*aa
        diag_1d_vol(j+1)=diag_1d_vol(j+1) +1D0*bb 
     endif

     ! neutral diagnosis
     if(diag_neu_pin < psi .and. psi < diag_neu_pout .and. (z>eq_x_z) .and. z<eq_x2_z .and. is_rgn12(r,z,psi)) then ! diagnosis excludes region 2 z<eq_x_z
        j=int((psi-diag_neu_pin)*diag_neu_dp_inv)+1
        j=min(diag_neu_npsi-1,max(1,j))
        bb=(psi-diag_neu_pin)*diag_neu_dp_inv + 1 -j
        aa=1D0-bb
        diag_neu_vol(j)=diag_neu_vol(j)+1D0*aa
        diag_neu_vol(j+1)=diag_neu_vol(j+1) +1D0*bb 
     endif

     if(col_mode/=0) then
        !for conserving collision
        if(col_mode==2) then
           if(col_2_pin < psi .and. psi < col_2_pout .and. z>eq_x_z .and. z<eq_x2_z .and. is_rgn12(r,z,psi)) then ! collision excludes region 2 z<eq_x_z : think more
              theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
              if(z < 0D0) then
                 theta= sml_2pi - theta    
              endif

              j=nint((psi-col_2_pin)*col_2_inv_dp)
              l=mod(nint(theta*col_2_inv_dtheta),col_2_mtheta)+1
              col_2_vol(l,j)=col_2_vol(l,j)+1D0
           endif
        endif

        !for background density calculation
        if(col_varying_bg==1) then
           if(col_vb_pin < psi .and. psi < col_vb_pout .and. z>eq_x_z .and. z<eq_x2_z .and. is_rgn12(r,z,psi)) then  ! region 2 z<eq_x_z : think more
              theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
              if(z < 0D0) then
                 theta= sml_2pi - theta    
              endif

              j=int((psi-col_vb_pin)*col_vb_inv_dp)
              j=min(col_vb_m-1,max(0,j))
              l=mod(nint(theta*col_vb_inv_dtheta),col_vb_mtheta)+1
              bb=(psi-col_vb_pin)*col_vb_inv_dp - j
              aa=1D0-bb
              col_vb_vol(l,j)=col_vb_vol(l,j)+1D0*aa
              col_vb_vol(l,j+1)=col_vb_vol(l,j+1)+1D0*bb
           endif
        endif
     endif
  

     if(.true.) then
!     if(sml_inpsi < psi .AND. psi < sml_outpsi .and. z>eq_x_z) then ! for all particles
        x(1)=r
        x(2)=z
        ! coordinate transform
        call field_following_pos2(x,phi,phi_dest,xff)
        

        !triangle search
        if (USE_SEARCH_TR2) then
           call search_tr2(grid,xff,itr,p)
        else
           call search_tr(grid,xff,itr,p)
        endif
        

        if(itr>0) then
           wphi=1D0 - (phi/grid%delta_phi)

           do j=1,3
              node=grid%nd(j,itr)
              grid%node_vol_ff(node,0) =  &
                   grid%node_vol_ff(node,0) + &
                   p(j)*wphi
              grid%node_vol_ff(node,1) =  &
                   grid%node_vol_ff(node,1) + &
                   p(j)*(1D0-wphi)
           enddo
           
           if(sml_f0_grid) then         
              ml=maxloc(p)
              node= grid%nd(ml(1),itr)
              grid%node_vol_nearest(node)= &
                   grid%node_vol_nearest(node) + 1D0
           endif
        endif
     endif
     
     
  enddo


  ! total sum
  dum1=real(total)
  
  call t_startf("GET_VOLUME_RED")
  call my_mpi_allreduce(dum1,dum2,1)
  call my_mpi_allreduce(psn%vol00,dum_00,grid%npsi00)
  psn%vol00=dum_00

  call my_mpi_allreduce(diag_1d_vol,dum_1d,diag_1d_npsi)
  diag_1d_vol=dum_1d

  call my_mpi_allreduce(diag_neu_vol,dum_neu,diag_neu_npsi)
  diag_neu_vol=dum_neu

  call my_mpi_allreduce(grid%node_vol_ff,dum_ff,grid%nnode*2)
  grid%node_vol_ff = dum_ff

  if(sml_f0_grid) then
     call my_mpi_allreduce(grid%node_vol_nearest, dum_ff(:,1), grid%nnode)
     grid%node_vol_nearest=dum_ff(:,1)
  endif

  if(col_mode/=0) then
      if(col_varying_bg==1) then
          call my_mpi_allreduce(col_vb_vol,dum_vb,col_vb_mtheta*(col_vb_m+1))
          col_vb_vol=dum_vb
      endif
      if(col_mode==2) then
          call my_mpi_allreduce(col_2_vol,dum_c,col_2_mtheta*(col_2_m+1))
          col_2_vol=dum_c
      endif
  endif
  call t_stopf("GET_VOLUME_RED")

  ! marker_den = density of simulation ptl = ptl_num*total_pe/(real volume * monte_num *totla_pe/sum of total )
  ! = ptl_num*sum of total /(real_volume*monte_num)

  if (sml_cylindrical) then
    mc_den=dum2/ (  rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) )
  else
    mc_den=dum2/ (  rdim*zdim*sml_2pi_wedge_n*(roffset+0.5D0*rdim) )
  endif
  sml_marker_den= mc_den  ! this is monte-carlo density -- real marker density will be set at the uniform_dist
  psn%vol00=psn%vol00 / mc_den

  !diag_1d
  diag_1d_vol=diag_1d_vol / mc_den
  
  !diag_neu
  diag_neu_vol=diag_neu_vol / mc_den

  if(col_mode/=0) then
      if(col_varying_bg==1) then
          col_vb_vol=col_vb_vol / mc_den
      endif
      if(col_mode==2) then
          col_2_vol=col_2_vol / mc_den
      endif
  endif

  ! node_vol_ff
  grid%node_vol_ff=grid%node_vol_ff/mc_den/real(sml_nphi_total)
  if(sml_f0_grid) then
     grid%node_vol_nearest=grid%node_vol_nearest/mc_den/real(sml_nphi_total)
  endif
  
  ! w0 calculation for local maxwell-----------
  
  deallocate(dum_ff)
end subroutine get_volume

! for single particle simulation. Load a particle in a given position and velocity.
subroutine load_single(spall)
  use ptl_module
  use sml_module
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  real (kind=8) :: rho,mu,en,pitch,b
  real (kind=8), external :: b_interpol,psi_interpol


  b=b_interpol(ptl_special_r,ptl_special_z,0D0)
  
  call ev_pitch_to_rho_mu2(ptl_special_en_ev,ptl_special_pitch,b,rho,mu,1)
  call rho_mu_to_ev_pitch2(rho,mu,b,en,pitch,1)
  print *, en,ptl_special_en_ev
  print *,rho,mu
  spall(1)%ptl(1)%ph(1)=ptl_special_r
  spall(1)%ptl(1)%ph(2)=ptl_special_z
  spall(1)%ptl(1)%ph(3)=ptl_special_phi
  spall(1)%ptl(1)%ph(4)=rho
  spall(1)%ptl(1)%ct(pim)=mu
  spall(1)%ptl(1)%ph(piw1:piw2)=1D-10
  spall(1)%ptl(1)%ct(piw0)=1D0
  spall(1)%ptl(1)%gid=sml_mype+1
  spall(1)%maxgid=sml_totalpe
  spall(1)%num=1

  if(sml_electron_on) then
     spall(0)%ptl(1)%ph(:)=spall(1)%ptl(1)%ph(:)
     spall(0)%ptl(1)%ph(4)=spall(1)%ptl(1)%ph(4)*sqrt(ptl_mass(0))/ptl_charge(0)
     spall(0)%ptl(1)%gid=sml_mype+1
     spall(0)%num=1
     spall(0)%maxgid=sml_totalpe
  endif
end subroutine load_single

!! ******************************************************
!!>functions for energy/pitch <-> rho/mu conversion
!!<*****************************************************
! ^^^^GPU related routine^^^^
subroutine rho_mu_to_ev_pitch2(rho,mu,b,ev,pitch,sp_type)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8), intent(inout) :: rho,mu,b
  real (kind=8), intent(out) :: ev,pitch
  real (kind=8) :: enj,v_pal,v_pep
  integer :: sp_type

  if(mu<0) then
     print *, 'minus mu found :',rho,mu,b
     mu=0D0
  endif
  
  enj=(mu*b+ptl_c2_2m(sp_type)*(rho*b)**2)
  ev=enj*sml_j2ev
  v_pal=ptl_c_m(sp_type)*rho*b
  v_pep=SQRT(2.d0*mu*b/ptl_mass(sp_type))
  pitch=v_pal/SQRT(v_pal**2+v_pep**2)
  
end subroutine rho_mu_to_ev_pitch2

subroutine ev_pitch_to_rho_mu2(ev,pitch,b,rho,mu,sp_type)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8), intent(in) :: ev,pitch,b
  real (kind=8), intent(out) :: rho,mu
  real (kind=8) :: en
  integer :: sp_type
  
  en=ev*sml_ev2j
  !       rho=pitch*DSQRT(2d0*en/mass)/b*mass/charge
  rho=pitch*SQRT(en/ptl_c2_2m(sp_type)) /B
  mu=max(0D0,(1d0-pitch**2)*en/b)
       
end subroutine ev_pitch_to_rho_mu2

!real (kind=8) function init_den_wz(psi,z)
!  use eq_module
!  use sml_module
!  use random
!  implicit none
!  real (kind=8), intent(in) :: psi,z
!  real (kind=8), external :: init_den,ranx
!  real (kind=8) :: tmp
!  if(z>eq_x_z .or. psi > eq_x_psi) then
!     tmp=psi
!  else
!     tmp=eq_x_psi
!  endif
!  init_den_wz=init_den(tmp)
!end function init_den_wz

!real (kind=8) function init_tempi_ev_wz(psi,z)
!  use eq_module
!  implicit none
!  real (kind=8), intent(in) :: psi,z
!  real (kind=8), external :: init_tempi_ev
!  real (kind=8) :: tmp
!  if(z>eq_x_z .or. psi > eq_x_psi) then
!     tmp=psi
!  else
!     tmp=eq_x_psi
!  endif
!  init_tempi_ev_wz=init_tempi_ev(tmp)
!end function init_tempi_ev_wz

!!$real (kind=8) function tempe_ev_wz(psi,z)
!!$  use eq_module
!!$  implicit none
!!$  real (kind=8), intent(in) :: psi,z
!!$  real (kind=8), external :: tempe_ev
!!$  real (kind=8) :: tmp
!!$  if(z>eq_x_z .or. psi > eq_x_psi) then
!!$     tmp=psi
!!$  else
!!$     tmp=eq_x_psi
!!$  endif
!!$  tempe_ev_wz=tempe_ev(tmp)
!!$
!!$end function tempe_ev_wz
!!$
!!$real (8) function init_ipara_flow(psi,z)
!!$  use eq_module
!!$  implicit none
!!$  real (8), intent(in) :: psi, z
!!$  real (8) :: tpsi
!!$  real (8) :: a,b
!!$  
!!$  if(z>eq_x_z .or. psi > eq_x_psi) then
!!$     tpsi=psi
!!$  else
!!$     tpsi=eq_x_psi
!!$  endif
!!$
!!$  a=0.5D0 *(eq_ipara_flow_edge-eq_ipara_flow_out)
!!$  b=0.5D0 *(eq_ipara_flow_edge+eq_ipara_flow_out)
!!$  init_ipara_flow= a*dtanh(2D0*(eq_ipara_flow_ped_c-tpsi)/eq_ipara_flow_ped_width) +b
!!$
!!$end function init_ipara_flow
!!$
!!$real (8) function init_epara_flow(psi,z)
!!$  use eq_module
!!$  implicit none
!!$  real (8), intent(in) :: psi, z
!!$  real (8) :: tpsi
!!$  real (8) :: a,b
!!$  
!!$  if(z>eq_x_z .or. psi > eq_x_psi) then
!!$     tpsi=psi
!!$  else
!!$     tpsi=eq_x_psi
!!$  endif
!!$  
!!$  a=0.5D0 *(eq_epara_flow_edge-eq_epara_flow_out)
!!$  b=0.5D0 *(eq_epara_flow_edge+eq_epara_flow_out)
!!$  init_epara_flow= a*dtanh(2D0*(eq_epara_flow_ped_c-tpsi)/eq_epara_flow_ped_width) +b
!!$
!!$end function init_epara_flow

#ifdef GAM_TEST
#ifndef POL_MODE
real (kind=8) function gam_weight(r,z,psi)
  use sml_module
  use eq_module
  implicit none
  real (kind=8) :: r, z, psi
  real (kind=8) :: a,b,c,k, r_minor2

  a=300D0
  b=sml_inpsi
  k=0.137/3.6D-3
  r_minor2=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  gam_weight=sml_initial_deltaf_noise *sin( k * sqrt(r_minor2) )


!  gam_weight=sml_initial_deltaf_noise*sin(a*(psi-b))
!  gam_weight=sml_initial_deltaf_noise

end function gam_weight


#else
real (kind=8) function gam_weight(r,z,psi)
  use sml_module
  use eq_module
  implicit none
  real (kind=8) :: r, z, psi
  real (kind=8) :: kr, kz


  kr=1D2
  kz=2D2

  gam_weight=sml_initial_deltaf_noise *sin( kr*r ) * sin( kz*z )


!  gam_weight=sml_initial_deltaf_noise*sin(a*(psi-b))
!  gam_weight=sml_initial_deltaf_noise

end function gam_weight
#endif


#endif


