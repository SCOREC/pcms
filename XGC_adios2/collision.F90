
subroutine collision(istep,sp)
  use ptl_module
  use col_module
  use perf_monitor
  implicit none
  integer, intent(in) :: istep
  type(species_type) :: sp
  integer :: flag


  if (col_mode /= 0) then
     if(col_varying_bg==1 .and. (istep==1 .or. mod(istep,col_vb_period)==0)) then
        call t_startf("COL_SNAPSHOT") 
        call col_snapshot(sp)
        call t_stopf("COL_SNAPSHOT") 
     endif

     if(mod(istep,col_period)==0) then
        flag=2* 2**2 - 1 ! all collision
        if(col_mode==1) then
           call t_startf("COLLISION1") 
           call collision1(sp, flag)
           call t_stopf("COLLISION1") 
        elseif (col_mode==2) then
           call t_startf("CONSERVING_COLLISION") 
           call conserving_collision(sp,istep,1)
           call t_stopf("CONSERVING_COLLISION")
        endif
     endif
  endif

end subroutine collision

!ion collision with no conservation
subroutine collision1(sp, iflag)
  use sml_module, only : sml_dt, sml_ev2j, sml_j2ev, sml_2pi, sml_nthreads
  use eq_module, only : eq_axis_r, eq_axis_z
  use ptl_module
  use omp_module, only : split_indices
  use col_module
  IMPLICIT NONE

  type(species_type) :: sp
  integer, intent(in) :: iflag

  real (kind=8), external :: b_interpol, psi_interpol

  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  integer :: iflag2, isp, i
  real (kind=8) :: c_m_sp, s_mass, s_charge
  real (kind=8) :: r, z, phi, psi, b, dt
  real (kind=8) :: dnb, ti_ev, up, theta
  real (kind=8) :: rho, mu, rho_m, accel_factor, ekmin
  real (kind=8) :: ekin, pitch

  dt = sml_dt * col_period
  iflag2= 1 + col_en_col_on*2    ! +1, pitch angle, +2 energy col
    
  do isp=ptl_isp , ptl_nsp
      c_m_sp=ptl_c_m(isp)
      s_mass=sp%mass
      s_charge=sp%charge
      ekmin = 1.d-3 * sml_ev2j  ! 1 eV minimum temperature for collision
      call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PHI, B, PSI, RHO, MU, &
!$OMP&         THETA, DNB, TI_EV, UP, RHO_M, EKIN, PITCH, &
!$OMP&         ACCEL_FACTOR) 
      do ith=1,sml_nthreads
         do i=i_beg(ith),i_end(ith)
             r=sp%ptl(i)%ph(1)
             z=sp%ptl(i)%ph(2)
             phi=sp%ptl(i)%ph(3)
             psi=psi_interpol(r,z,0,0)
             b=b_interpol(r,z,phi)
             rho=sp%ptl(i)%ph(pirho)
             mu=sp%ptl(i)%ct(pim)

             theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
             if(z < 0D0) then
                theta= sml_2pi - theta
             endif

             call background_ion_profile(theta,r,z,psi,dnb,ti_ev,up)

             if (col_moving_frame) then
                rho_m=rho - up/(c_m_sp*b)
             else
                rho_m=rho
             endif

             !test ptl collision
             call rho_mu_to_ev_pitch2(rho_m,mu,b,ekin,pitch,sp%type) ! ekin is in eV unit
             ekin=ekin*sml_ev2j         ! conversion to SI unit
             ekin = max(ekmin, ekin)        ! 2002/09/17 added for preventing small energy
                
             if(col_accel) then
                call col_accel_factor(psi,accel_factor)
             else
                accel_factor=1D0
             endif
                
             call scatr_one(ekin, pitch, s_mass, s_charge, dnb, ti_ev, ptl_mass_au, ptl_charge_eu,&
                            accel_factor, dt, ekmin, iflag2)
                 
             call ev_pitch_to_rho_mu2(ekin*sml_j2ev, pitch, b, rho_m, mu,sp%type)

             if (col_moving_frame) then
                sp%ptl(i)%ph(pirho)=rho_m + up/(c_m_sp*b)
             else
                sp%ptl(i)%ph(pirho)=rho_m
             endif
             sp%ptl(i)%ct(pim)=mu
          enddo
      enddo
  enddo

end subroutine collision1

subroutine scatr_one(ekin, pitch, massa, chargea, denb, tempb_ev, massb_au, chargeb_eu,&
                     accel, dt, ekmin, iflag2)
  use sml_module, only : sml_ev2j
  use random_xgc
  IMPLICIT NONE

  real (kind=8) :: ekin, pitch, massa, chargea
  real (kind=8) :: denb, tempb_ev, massb_au, chargeb_eu
  real (kind=8) :: accel, dt, ekmin
  integer, intent(in) :: iflag2

  real (kind=8) :: colb, colbs, fac0b
  real (kind=8) :: dum,agg, del_pitch, esig

! iflag2: +3 pitch angle + energy

  call find_freq(ekin, massa, chargea, denb, tempb_ev, &
                massb_au, chargeb_eu, colb, colbs, fac0b,accel)

  if(mod(iflag2,2)==1) then  !pitch angle scattering
     !!print*, 'ion--pitch called'
     agg = ranx() - .5d0
     dum = 1.d0 - pitch**2
     dum = max(0.d0, dum)
     del_pitch = dsign(1.d0,agg)*sqrt(dum*colb*dt*0.5d0)
     pitch = pitch*(1.d0-colb*dt*0.5d0) + del_pitch
  endif
  if(mod(iflag2/2,2)==1) then !energy collision
     !!print*, 'ion--energy called'
     agg = ranx()- .5d0
     esig = dsign(1.d0,agg)*sqrt(2*ekin*tempb_ev*sml_ev2j*colbs*dt)
     ekin = ekin - colbs*dt*fac0b + esig
     ekin = max(ekin,ekmin)
     if(ekin > 1D0 .or. ekin < 2D0 )  then
           
     else
        print *,'colbs',colbs, fac0b, esig
        ekin=ekmin
        stop
     endif
  endif
   
  pitch = min(1.D0,pitch)
  pitch = max(-1.D0,pitch)
end subroutine scatr_one


! Warning : This routine uses cgs unit partially
subroutine find_freq(en_a, mass, charge, dn_b, en_b_ev, mass_b, charge_b, freq_scat, freq_slow, freq_fac0,accel)
  ! calculate collision frequencies and convert into MKS unit
  ! approx to psi function
  ! pitch angle scattering, small value of col!  (alpha)
  ! Boozer and Kuo-Petravic Phys. Fluids 24, 851 (1981)
  ! profiles-density (cm-3) dnb (background), dni (impurity), temp (kev)
  ! psi(x) = (2/sqrt(pi)Int[dt*sqrt(t)*exp(-t)]   
  ! psi(x) = Phi(v/vth) = Phi(sqrt(x)) in the paper
  ! psi(x) = 1 - 1/(1 + p), approximatly
  ! p = x**1.5[a0 + a1*x +a2*x**2 +a3*x**3 +a4*x**4]
  ! R.White,M.Redi Jan (2000) error in psi(x) < 1.D-3, asymptotically correct
  ! relative error dpsi/psi < 1.D-2
  ! en_a(J), mass(kg), charge(C), dn_b (m^-3), en_b_ev (eV), mass_b (Atomic Unit), charge_b (Electron Charge unit)
  ! 
  use sml_module !, only : sml_prot, sml_zprt, sml_en_order_kev, sml_en_order
  implicit none
  real(kind=8), intent(in) :: en_a, mass, charge, mass_b, charge_b, dn_b, en_b_ev,accel !, mass_b, charge_b
  real(kind=8), intent(out) :: freq_scat, freq_slow, freq_fac0
  real(kind=8) :: dn, dumb, dd,d3, dum, ap0, ap1, ap2, ap3, ap4, ee_ev, cnst, vprt, &
       & vt, ap_psi, f, g, gp, bmax, bmin1, massab, bmin2, bmin, clog, dnu_b
  data ap0 /.75225/,ap1 /-.238/,ap2 /.311/,ap3 /-.0956/,ap4 /.0156/
  real (8) :: mass_au, charge_eu

  vprt = sqrt(2D0*en_a/mass)  ! MKS - m/s
  mass_au = mass/sml_prot_mass
  charge_eu=charge/sml_e_charge
!!  vt = dsqrt(2d0*(en_b_ev/1.d3*sml_en_order/sml_en_order_kev)*mass/mass_b)
  vt = sqrt(2.d0*(en_b_ev*sml_ev2j)/(mass_b*sml_prot_mass))  ! MKS - m/s

!!  norm_r_cgs = nc_norm_r * 100.d0 !!! conversion MKS to cgs
!!  cnst = 2.41D11*(charge/mass)**2/(norm_r_cgs/nc_norm_t*vprt)**3  !100 is conversion MKS to CGS
!!  cnst = cnst*col_accel

  cnst = 2.41D11*(charge_eu/mass_au)**2/(100d0*vprt)**3  !100 is conversion MKS to CGS
  ! cnst = mu_b / x^3 /Clog /den / some order 1 const
  ! mu_b=Braginskii collision frequency
  ! x = v/vth, Clog = Coulumb Logarithm
  ! 2.4D11 = e(cgs)^4 / m(cgs)^2 * 4*pi
  cnst = cnst*accel
  dn = dn_b/1d6             ! dn - CGS cm^-3



  ee_ev = en_a*sml_j2ev   !!! test ptl's kinetic energy in [eV]

  ! calculate psi(x), f, g, gp
  dumb = vprt/vt          ! MKS/MKS  ,  v/v_th
  dd = dumb**2   ! dd = x, dumb = sqrt(x)
  d3 = dumb**3
  dum = d3*(ap0 +ap1*dd +ap2*dd**2 +ap3*dd**3 +ap4*dd**4)
  ap_psi = 1.D0 - 1.D0/(1.D0 + dum) !psi(x) = Phi(v/vth)

  f = (2.D0 - 1.D0/dd)*ap_psi + 2.257D0*dumb*dexp(-dd) 
  g = 2D0*mass_au*ap_psi/mass_b
  gp = 2.257D0*dexp(-dd)*(mass_au/mass_b)*dumb

  ! find Coulomb logarithm
!!  bmax  = 7.4d-3*dsqrt(en_b_ev/1000d0*1.d13/(charge**2*dn_b))
!!  bmin1 = 1.4d-10*charge*charge_b/(ee_ev + en_b_ev)*1000d0
!!  massab = mass_au*mass_b/(mass_au+mass_b)
!!  bmin2 = 6.2d-10/(dsqrt(massab*1836.d0*en_b_ev/1000d0))


  bmax  = 7.4d-3*sqrt(en_b_ev*1.d10/(charge_eu**2*dn))
  bmin1 = 1.4d-7*charge_eu*charge_b/(ee_ev + en_b_ev)
  massab = mass_au*mass_b/(mass_au+mass_b)
  bmin2 = 6.2d-10/(sqrt(massab*1836.d0*en_b_ev/1.d3))
  bmin = max(bmin1,bmin2)
  clog = dlog(bmax/bmin)

  ! collision frequencies - scattering, slowing down
  dnu_b = cnst*dn*charge_b**2
  freq_scat = dabs(clog*dnu_b*f)
  freq_slow = dabs(clog*dnu_b*g)
  freq_fac0 = en_a*(1.d0 - mass_b*gp/(g*mass_au)) 

  return

end subroutine find_freq

subroutine background_ion_profile(theta,r,z,psi,deni,tempi,up)
  use sml_module, only : sml_2pi
  use eq_module
  use col_module
  implicit none

  real (kind=8), intent(in) :: theta,r,z,psi
  real (kind=8), intent(inout) :: deni, tempi, up

  integer :: j, l, isp
  real (kind=8) :: aa, bb, tmp

  isp=1 !for ion species

  if (col_varying_bg==1 .and. col_vb_pin<psi .and. psi<col_vb_pout .and. z>eq_x_z .and. is_rgn12(r,z,psi)) then
      tmp=(psi-col_vb_pin)*col_vb_inv_dp
      j=int(tmp)
      j=min(col_vb_m-1,max(0,j))
      l=mod(nint(theta*col_vb_inv_dtheta),col_vb_mtheta)+1
      bb=tmp - j
      aa=1D0-bb

      deni=col_vb(1,l,j,isp)*aa+col_vb(1,l,j+1,isp)*bb
      tempi=col_vb(2,l,j,isp)*aa+col_vb(2,l,j+1,isp)*bb

      j=nint((psi-col_vb_pin)*col_vb_inv_dp)
      up=col_vb(3,l,j,isp)
  else 
      deni=eq_ftn(psi,r,z,eq_den)
      tempi=eq_ftn(psi,r,z,eq_tempi)
      up=eq_ftn(psi,r,z,eq_flowi)
  endif

end subroutine background_ion_profile

!get ion density and temperature for collision routine
subroutine col_snapshot(sp)
  use sml_module,only : sml_mype, sml_j2ev, sml_2pi, sml_deltaf, sml_nthreads, sml_istep
  use eq_module, only : eq_x_psi, eq_axis_r, eq_axis_z, is_rgn12
  use ptl_module
  use col_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  type(species_type) :: sp
  integer :: i,j,l, isp
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  real (kind=8) :: dum(5,col_vb_mtheta,0:col_vb_m,ptl_isp:ptl_nsp,sml_nthreads)
  real (kind=8) :: dum2(5,col_vb_mtheta,0:col_vb_m,ptl_isp:ptl_nsp)
  real (kind=8) :: B,c_m_ion, c2_2m_ion, v(3)
  real (kind=8) :: psi,r,z,weight,phi,aa,bb, theta, tmp
  real (kind=8),external :: psi_interpol, b_interpol

  isp =1 ! for ion
  col_vb(:,:,:,isp)=1D-99  ! to avoid divide by zero
  dum(:,:,:,isp,:)=1D-99  ! to avoid divide by zero
  c_m_ion=ptl_c_m(isp)
  c2_2m_ion=ptl_c2_2m(isp)

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PHI, &
!$OMP          PSI, WEIGHT, B, &
!$OMP          V, THETA, TMP, J, L, &
!$OMP          AA, BB)
  do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
          r=sp%ptl(i)%ph(1)
          z=sp%ptl(i)%ph(2)
          phi=sp%ptl(i)%ph(3)

          if(sp%ptl(i)%gid > 0) then
              psi=psi_interpol(r,z,0,0)
              if(col_vb_pin<psi .AND. psi<col_vb_pout .AND. is_rgn12(r,z,psi) ) then ! 2002/10/10 , region 2 condition modified
                  if(ptl_deltaf_sp(sp%type) .and. col_varying_bg==1 .and. .not. col_mode .eq. 3) then !col_mode==3 is under construction & not sure it needs col_snapshot 
                      print *,'deltaf collision is not implemented yet.'
                      stop
!                      weight=ptl%ion%phase(6,i)*ptl%ion%phase(8,i)
                  else
                      weight=sp%ptl(i)%ct(piw0)
                  endif

                  !energy calculation
                  B=b_interpol(r,z,phi)

                  v(1)=weight
                  v(3)=sp%ptl(i)%ph(pirho)*B
                  v(2)=weight*(c2_2m_ion*v(3)**2 + sp%ptl(i)%ct(pim)*B)
                  v(3)=weight*v(3)

                  theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
                  if(z < 0D0) then
                      theta= sml_2pi - theta
                  endif

                  tmp=(psi-col_vb_pin)*col_vb_inv_dp
                  j=int(tmp)
                  j=min(col_vb_m-1,max(0,j))
                  l=mod(nint(theta*col_vb_inv_dtheta),col_vb_mtheta)+1
                  bb=tmp - j
                  aa=1D0-bb
                  
                  dum(1:3,l,j,isp,ith)=dum(1:3,l,j,isp,ith)+v*aa
                  dum(1:3,l,j+1,isp,ith)=dum(1:3,l,j+1,isp,ith)+v*bb

                  j=nint((psi-col_vb_pin)*col_vb_inv_dp)
                  dum(4,l,j,isp,ith)=dum(4,l,j,isp,ith)+v(1)
                  dum(5,l,j,isp,ith)=dum(5,l,j,isp,ith)+v(3)
              endif
          endif
      enddo
  enddo

  do ith=2,sml_nthreads
      dum(:,:,:,isp,1)=dum(:,:,:,isp,1)+dum(:,:,:,isp,ith);
  enddo

  dum(3,:,:,isp,1)=dum(3,:,:,isp,1)*c_m_ion  !parallel velocity to MKS
  dum(5,:,:,isp,1)=dum(5,:,:,isp,1)*c_m_ion  !parallel velocity to MKS

  call t_startf("COL_SNAP_RED") 
  call my_mpi_allreduce(dum(:,:,:,isp,1),dum2(:,:,:,isp),5*col_vb_mtheta*(col_vb_m+1))
  call t_stopf("COL_SNAP_RED") 

  col_vb(2,:,:,isp)=2D0/3D0*(dum2(2,:,:,isp)-0.5D0*ptl_mass(isp)*(dum2(3,:,:,isp)/dum2(1,:,:,isp))**2)/dum2(1,:,:,isp)*sml_j2ev !Temperature in (eV)
  col_vb(1,:,:,isp)=dum2(1,:,:,isp)/col_vb_vol ! density
  col_vb(3,:,:,isp)=dum2(5,:,:,isp)/dum2(4,:,:,isp) ! mean flow for moving frame collision

end subroutine col_snapshot
