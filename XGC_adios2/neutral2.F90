
!!***************************************************
!! 2-dimensional neutral particle simulation routine
!! After each evaluation of the plasma profile, stop the plasma routine
!! and evaluate neutral density and temperature profiles for the same physical
!! duration.
!! Assume a hypothetical boundary at Psi_hypo=1.03 or z=zmin.  At Psi >= Psi_hypo
!! or z<= zmin, assume Maxwellian neutral distribution function at room 
!! temperature.
!!**************************************************

subroutine neutral_step(mass)
  use sml_module
  use eq_module,only: eq_x_psi,eq_axis_r,eq_axis_z
  use neu_module
  use sml_module
  use ptl_module
  use perf_monitor
  use random_xgc
  use Ecuyer_random
  implicit none
  real (kind=8), intent(in) :: mass
  real (kind=8) :: rstart,zstart,slope,r0,z0,r,z,vr,vz,vangle,tn,v,pol,vparan,vparai
  integer :: ipol,k,ipart,ptest,istep
  real (kind=8) :: f_ion, f_cx, f_el,rate_ion,rate_cx,rate_el,drate
  real (kind=8) :: denpe,denpi,teev,tiev,enev,xmax,sigma0,jacob,tekev
  real (kind=8) :: weight, den, del,theta,psi_interpol
!  integer , parameter :: mpol=1000
! neutral-ion Elastic scattering crosssection
!  print *, 'neutral_step'

!  eq_axis_z=eq_axis_z*0.
  if(sml_mype==0) print *, 'Neutral profile renewal...'
  sigma0=5d-19/sqrt(2D0)  !m**2

  neu_grid_den=0D0
  neu_grid_temp=0D0
  neu_grid_flow=0D0
  neu_grid_wsum=0D0
  neu_grid_esum=0D0
  neu_grid_vsum=0D0
call check_point('before init_seeds')
call init_seeds(1234*sml_mype,2345*sml_mype+6789,4321*sml_mype+10)
!  slope=sml_pi*0.5D0
  r0=eq_axis_r 
  z0=eq_axis_z
if(sml_mype==0)print *, 'ranx()=',ranx()
  ! Make the neutral time equal to the plasma time interval
  ! Initiation position
!	if(sml_mype==0)print *, 'neu_mpol==',neu_mpol,istep
  do ipol=1,neu_mpol

      !print *, ipol
     pol=2.d0*sml_pi*(real(ipol-1)+ real(sml_mype)/real(sml_totalpe))/real(neu_mpol)  ! mype is for phase difference
!	if(sml_mype==0)print *, 'neu_num==',neu_num,pol
     ! start simulation at each initiation position
     do ipart=1, neu_num
        ! neut_num is the particle number per processor
        call startrz(pol,rstart,zstart,jacob,slope)
        ! Get (R,Z) of psi=1.03 at the starting position.
        r=rstart
        z=zstart
        tn=neu_temp0 ! unit ev
!if(sml_mype==0)print *,'r,z,pol',r,z,jacob,slope	
        call maxwell_neut(r0,z0,r,z,slope,tn,vr,vz,vangle)    
        v=sqrt(vr**2+vz**2)
        vparan=0D0 ! initial parallel velocity of test neutral particle 

        ! get theta value and neutral density
        theta=dacos( (r-eq_axis_r)/dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2) )
        if(z<eq_axis_z) then
           theta=sml_2pi-theta
        endif
        del=min(abs(theta-neu_peak_theta),min( abs(theta-neu_peak_theta+sml_2pi), abs(theta-neu_peak_theta-sml_2pi)))
        den=(1D0 +(neu_delta_n-1.d0)*exp( -(del/neu_delta_theta)**2 ) ) ! density normalize to neu_n0

        ! weight calculation -----------------------------------
        weight=jacob*sml_2pi/real(neu_mpol)*  &  ! area of segment
             den* &  ! neutral density in norm unit
             v*abs(sin(slope-vangle))*neu_dt  & !
             /real(neu_num*sml_totalpe) !---------------------
!if(sml_mype==0)print *,'vr,vz,theta',vr,vz,den,weight,neu_dt,v,neu_num,sml_totalpe,ranx(),ipart	
        ! Record neutral density and temperature at neutral birth position 
        call neut_fluid(r,z,vr,vz,weight,vparan) 

        do istep=1, neu_istep_max
           ! move forward.
           r= r+vr*neu_dt
           z= z+vz*neu_dt

           ! If particle is outside the hypothetical boundary, stop simulation.
           ! If ptest=-2, neutral entered divertor.
           ! If ptest=-1, neutral is on or outside Psi = 1.03.
           call neut_pos_test(r,z,ptest)
!           if(sml_mype==0) write(143,*) r,z
           if(ptest<0) then
              !	   if(sml_mype==0)print *,'r,z=',r,z,istep
              exit
           endif
           call plasma(r,z,denpe,denpi,teev,tiev,vparai,istep)  ! te, ti-- ev, denp -- mks
	   
!           if(sml_mype==0) print *,'denp,te', denpi,teev,tiev
           !------------------------------------------------
           ! Get local collision probability of the neutral particle
	   !
           ! Probability of ionization per sec 
           f_ion=denpe*0.8D-8*sqrt(Teev)*(exp(-13.56/Teev)/(1D0+0.01D0*Teev))*1D-6  
           ! denp is in MKS here.
           f_ion=f_ion! s^-1 --> conversion to normalized unit
           !
           ! probability of charge exchange per sec 
           enev=tiev
           f_cx= denpi * 1.1D-8* (enev)**0.3 /sqrt(ptl_mass_au) *1D-6   
           ! 1.1 K_i(eV)^0.3 * 10^-8 / Sqrt(Mi(AU))  cm s^-1   -> 1D-6 m s^-1
           f_cx=f_cx  ! s^-1 --> conversion to normalized unit. 
           !
           ! probablity of elastic collision
           f_el = denpi*sigma0*v*0.5D0*sml_pi  !check !!!!
           !
           ! -----------------------------------------------------
           ! Rate parameters
           drate= neu_dt*f_ion
!           rate_ion=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
           rate_ion=1D0-dexp(-drate)
!	   rate_ion=drate
           drate= neu_dt*f_cx
!           rate_cx=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
           rate_cx=1D0-dexp(-drate) 
!	   rate_cx=drate
           drate= neu_dt*f_el
!           rate_el=drate-drate**2/2.d0+drate**3/6.d0-drate**4/24.d0
           rate_el=1D0-dexp(-drate)
!	   rate_el=drate
              
           !write(123,*) rate_ion,rate_cx,rate_el
!           rate_cx=0D0; rate_el=0D0; rate_ion=0D0

           ! Ionization of the neutral particle.  
           !  The neutral disappears upon ionization.
!           if(rate_ion .ge. ranx().and. psi_interpol(r,z,0,0) < NC_x_psi) exit  !2003/08/12 - ionize inside only
           if(rate_ion .ge. ranx()) exit  !2003/08/12 - ionize everywhere

           ! Charge Exchange of the neutral particle.
           ! Vangle is randomized and te energy jumps to local Ti.
           if(rate_cx .ge. ranx()) then
              vangle=vangle+2.d0*sml_pi*ranx()
              v=dsqrt(2.d0*tiev*sml_ev2j/mass)
              vr=v*dcos(vangle)
              vz=v*dsin(vangle)
              vparan=vparai
           endif

           ! Elastic collision of the neutral particle.
           ! Vangle is randomized keeping the same energy.
           if((rate_el .ge. ranx()).and.(neu_elastic_col_on==1)) then
              vangle=vangle+2.d0*sml_pi*ranx()
              vr=v*dcos(vangle)
              vz=v*dsin(vangle)
           endif

           call neut_fluid(r,z,vr,vz,weight,vparan)  ! Record neutral density and temperature
        enddo
     enddo
!     if(sml_mype==0)print *, 'collision rate=',f_ion,f_cx,f_el,teev,enev,denpi,denpe
  enddo

  call t_startf("NEUT_FLUID_FINAL")
  call neut_fluid_final
  call t_stopf("NEUT_FLUID_FINAL")
  if(sml_mype==0) print *, 'Neutral profile modified.'
end subroutine

subroutine maxwell_neut(r0,z0,r,z,slope,tn,vr,vz,vangle)
  !***********************************************************
  ! SLOPE is the local slope of the imaginary wall at Psi = 1.03
  ! Created on 5/17/2003
  !**********************************************************
  use sml_module
  use ptl_module
  use random_xgc
  implicit none
  real(kind=8) , intent(in) :: r,z,r0,z0,slope, tn
  real(kind=8) , intent(out) :: vr, vz, vangle
  real (kind=8) :: xdum,xmax, v,tn2
  ! r-z velocity
  xmax=1.d0-dexp(-7.d0)      
  xdum=xmax*ranx()
  tn2=tn*sml_ev2j/ptl_mass(1)
  v=dsqrt(-2.d0*dlog(1.d0-xdum)*tn2)
!if(sml_mype==0)print *, 'energy=',tn,0.5*ptl_mass(1)*v**2*sml_j2ev
  vangle=2.d0*sml_pi*ranx()
!  vangle=slope+0.5D0*sml_pi  !debug   
  vr=v*dcos(vangle)
  vz=v*dsin(vangle)
!  if( (r.gt.r0.and. (slope .ge. vangle .or. vangle .ge. (slope+sml_pi))) .or. &
!  (r.le.r0 .and. slope .le. vangle .and. vangle .le. (slope+sml_pi))) then
!     vangle=vangle+sml_pi
!     vr=v*dcos(vangle)
!     vz=v*dsin(vangle)
!  endif

end subroutine maxwell_neut

subroutine neut_fluid(r,z,vr,vz,weight,vparan)
  use neu_module
  use eq_module
  use sml_module
  use ptl_module,only : ptl_mass
  implicit none
  real (kind=8),intent(in) :: r,z,vr,vz,weight,vparan
  real (kind=8), external :: psi_interpol,neu_psi_interpol
  real (kind=8) :: psi,psi1,theta,rs,a1,a2,energy,pflow
  integer :: ipsi,itheta,itheta_p1

!  eq_axis_z=eq_axis_z*0.
  psi1=psi_interpol(r,z,0,0)
  rs=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  theta=acos( (r-eq_axis_r)/sqrt(rs) )
  if(z<eq_axis_z)  theta = sml_2pi -theta
  if(neu_grid_mode>0) then
     psi=neu_psi_interpol(r,z,psi1)
  else
     psi=psi1
  endif

!  if( psi  > neu_grid_min_psi .and. psi<= neu_grid_max_psi*1.0001D0 .and. z>eq_x_z) then
  if( psi  > neu_grid_min_psi .and. psi<= neu_grid_max_psi*1.0001D0 .and. ((z>eq_x_z) &
     .or.(neu_grid_mode>0))) then
     
     ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
     a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi + 1
     ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
     
     itheta= int(theta/neu_grid_dtheta) + 1 
     itheta= min(neu_grid_mtheta, max(1, itheta))
     a2= theta/neu_grid_dtheta - itheta +1
     itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1

     
     ! sum for neutral density
     neu_grid_wsum(ipsi,itheta)=neu_grid_wsum(ipsi,itheta) + weight*(1D0-a1)*(1D0-a2)
     neu_grid_wsum(ipsi+1,itheta)=neu_grid_wsum(ipsi+1,itheta) + weight*a1*(1D0-a2)
     neu_grid_wsum(ipsi,itheta_p1)=neu_grid_wsum(ipsi,itheta_p1) + weight*(1D0-a1)*a2
     neu_grid_wsum(ipsi+1,itheta_p1)=neu_grid_wsum(ipsi+1,itheta_p1) + weight*a1*a2

     energy=0.5D0*(vr**2 + vz**2)*weight*ptl_mass(1)
     
     ! sum for neutral energy
     neu_grid_esum(ipsi,itheta)=neu_grid_esum(ipsi,itheta) + energy*(1D0-a1)*(1D0-a2)
     neu_grid_esum(ipsi+1,itheta)=neu_grid_esum(ipsi+1,itheta) + energy*a1*(1D0-a2)
     neu_grid_esum(ipsi,itheta_p1)=neu_grid_esum(ipsi,itheta_p1) + energy*(1D0-a1)*a2
     neu_grid_esum(ipsi+1,itheta_p1)=neu_grid_esum(ipsi+1,itheta_p1) + energy*a1*a2  ! bug fix 2003/06/07

     pflow=vparan*weight

     ! sum for neutral parallel flow
     neu_grid_vsum(ipsi,itheta)=neu_grid_vsum(ipsi,itheta) + pflow*(1D0-a1)*(1D0-a2)
     neu_grid_vsum(ipsi+1,itheta)=neu_grid_vsum(ipsi+1,itheta) + pflow*a1*(1D0-a2)
     neu_grid_vsum(ipsi,itheta_p1)=neu_grid_vsum(ipsi,itheta_p1) + pflow*(1D0-a1)*a2
     neu_grid_vsum(ipsi+1,itheta_p1)=neu_grid_vsum(ipsi+1,itheta_p1) + pflow*a1*a2

  endif
  
!  write(141,*) vparan,pflow
  
end subroutine neut_fluid


subroutine neut_fluid_final
  use neu_module
  use sml_module
  use eq_module
  use perf_monitor
  implicit none
  real (kind=8) :: wsum(neu_grid_mpsi,neu_grid_mtheta), esum(neu_grid_mpsi, neu_grid_mtheta)
  real (kind=8) :: vsum(neu_grid_mpsi,neu_grid_mtheta)
  real (kind=8) :: del,theta,den
  integer :: i,j
  real (kind=8) , external :: tempi_ev
  save wsum, esum

  call t_startf("NEUT_FLUID_FINAL_RED")
  call  my_mpi_allreduce(neu_grid_wsum, wsum, neu_grid_mtheta*neu_grid_mpsi)
  call  my_mpi_allreduce(neu_grid_esum, esum, neu_grid_mtheta*neu_grid_mpsi)
  call  my_mpi_allreduce(neu_grid_vsum, vsum, neu_grid_mtheta*neu_grid_mpsi)
  call t_stopf("NEUT_FLUID_FINAL_RED")
!   if(sml_mype==0)print *, 'neu_grid_vsum=',neu_grid_vsum
!   if(sml_mype==0)print *, 'neu_grid_wsum=',vsum
!   if(sml_mype==0)print *, 'neu_grid_esum=',esum
!   if(sml_mype==0)print *, 'neu_grid_vol=',neu_grid_vol   
  neu_grid_den=wsum/neu_grid_vol
  neu_grid_temp=esum/wsum
!  neu_grid_den(neu_grid_mpsi,:)=2D0*neu_grid_den(neu_grid_mpsi-1,:)-neu_grid_den(neu_grid_mpsi-2,:)
!  neu_grid_den(neu_grid_mpsi,:)=neu_grid_den(neu_grid_mpsi-1,:)
  neu_grid_flow=vsum/wsum
  neu_nr=0D0
  do j=1, neu_grid_mtheta
     do i=1, neu_grid_mpsi
        if(neu_grid_den(i,j)<1D-10) then
!           neu_grid_temp(i,j)=tempi_ev(neu_grid_dpsi*real(i-1)+neu_grid_min_psi,0D0,0)*sml_ev2j 
           neu_grid_temp(i,j)=neu_temp0*sml_ev2j 
           neu_grid_flow(i,j)=0D0
        endif
     enddo
     neu_nr=neu_nr+neu_grid_den(neu_grid_mpsi,j) 
  enddo

  neu_nr=neu_nr/real(neu_grid_mtheta)
!   if(sml_mype==0)print *, 'neu_temp0=',neu_temp0*sml_ev2j
!  do j=1, neu_grid_mtheta
!     theta=neu_grid_dtheta*real(j-1)
!     del=min(abs(theta-neu_theta_x),min( abs(theta-neu_theta_x+sml_2pi), abs(theta-neu_theta_x-sml_2pi)))
!     den=(1D0 +(neu_delta_n-1.d0)*exp( -(del/neu_delta_theta)**2 ) ) ! density normalized to neu_n0
!     neu_grid_den(:,j)=neu_grid_den(:,j)*den/neu_grid_den(neu_grid_mpsi,j)
!  enddo
!   if(sml_mype==0)print *, 'sml_ev2j=',neu_temp0,neu_grid_temp
  if(sml_mype==300) then
     do j=1,neu_grid_mtheta 
        do i=1,neu_grid_mpsi
           write(145,8000) (neu_grid_dpsi*real(i-1)+neu_grid_min_psi)/eq_x_psi,neu_grid_dtheta*real(j-1),&
              neu_grid_den(i,j),neu_grid_temp(i,j)*sml_j2ev,neu_grid_flow(i,j),&
              neu_grid_vol(i,j),wsum(i,j)
        enddo
        write(145,*) ' '
     enddo
  endif
  close(145)

8000 format(7(e19.13,' '))
end subroutine


subroutine plasma(r,z,denpe,denpi,teev,tiev,vparai,istep)
  use sml_module
  use eq_module
  use neu_module,only : neu_istep_max
  use ptl_module, only : ptl_charge,ptl_charge_eu 
  use col_module
  implicit none
  real (kind=8) ,intent(in):: r,z
  real (kind=8) , intent(out) :: denpe,denpi,teev,tiev,vparai
  real (kind=8) , external :: psi_interpol,den_ion,den_elec,flow_species
  real (kind=8) :: psi,xd,aa,bb,tmp,theta

  integer :: istep, j, l, isp
  isp=1 !for ion species

  psi=psi_interpol(r,z,0,0)
  theta=dacos( (r-eq_axis_r)/dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2) )

!  call background_ion_profile(theta,r,z,psi,denpi,tiev,vparai)

  denpi=eq_ftn(psi,r,z,eq_den)
  tiev=eq_ftn(psi,r,z,eq_tempi)

  if(sml_electron_on) then
     denpe=eq_ftn(psi,r,z,eq_den) ! collision.f90   MKS  m^-3
  else
     denpe=ptl_charge_eu*denpi ! collision.f90   MKS  m^-3
  endif

  teev=eq_ftn(psi,r,z,eq_tempe)  ! collision.f90 -- ev

  teev=max(teev,1D-9)
  tiev=max(tiev,1D-9)
  denpe=max(denpe,1D0)
  denpi=max(denpi,1D0)
!  vparai=flow_species(psi,z,1)
  vparai=0D0 
!  denp0,denp1,te0,te1,ti0,ti1,width,rs

  ! Plasma density and temperature in a test case
  ! ne=ni assumed

end subroutine plasma

subroutine  neut_pos_test(r,z,ptest)

  use sml_module
  use eq_module, only : eq_axis_r, eq_axis_z 
  implicit none
  real (kind=8) :: r,z
  integer :: ptest
  real (kind=8) :: theta, rw,zw,tmp1,tmp2,rs
  
!  eq_axis_z=eq_axis_z*0.
  rs=(r-eq_axis_r)**2 + (z-eq_axis_z)**2
  theta=dacos((r-eq_axis_r)/dsqrt(rs))
  if(z<eq_axis_z) then
     theta= sml_2pi - theta
  endif
  
  call startrz(theta,rw,zw,tmp1,tmp2)
  if( rs*0.9999D0 <= (rw-eq_axis_r) **2 +(zw-eq_axis_z)**2 ) then
     ptest=1
  else
     ptest=-1
  endif
end subroutine neut_pos_test


subroutine startrz(theta,r,z,jacob,slope)
  use neu_module
  use sml_module
  implicit none
  real (kind=8) ,intent(in) :: theta
  real (kind=8) ,intent(out) :: r,z,jacob,slope
!  integer ,intent(out) :: flag
  integer :: i,j
  real (kind=8) :: aa,dr,dz
!  flag=1
!  ! check if theta is in valid range
!  if( neu_theta1 < theta .AND.  theta <neu_theta2 ) then
!     flag=0
!     return
!  endif

  ! find proper index of pol
  i= int(theta/neu_dtheta) + 1
  i=min(neu_mtheta,max(1,i))
  aa= theta/neu_dtheta -real(i-1)

  j=i+1
  if(i+1 > neu_mtheta)  j=1 ! cyclic

  r=(1D0-aa)*neu_r0(i) + aa*neu_r0(j)
  z=(1D0-aa)*neu_z0(i) + aa*neu_z0(j)
  jacob=(1D0-aa)*neu_jacob(i) + aa*neu_jacob(j)
!  jacob=1D0
  dr=neu_r0(j)-neu_r0(i)
  dz=neu_z0(j)-neu_z0(i)
  slope= acos( dr/sqrt(dz**2+dr**2) )
  if(dz<0) then
     slope= sml_2pi - slope 
  endif

end subroutine startrz


subroutine neutral2_setup
  use neu_module
  use eq_module
  use sml_module
  use lim_module
  implicit none
  real (kind=8):: max_r,min_r,max_z,min_z,theta,zd,rd,right,left,mid,r_tmp,z_tmp,&
       max_psi,max_psi_r,max_psi_z,factor,psi,psi1,a2,tantheta
  real (kind=8) ,external :: psi_interpol,neu_psi_interpol
  integer :: i,j,ip,im,ifactor,itheta,itheta_p1

!  eq_axis_z=eq_axis_z*0.
  ! find set of points equally spaced in poloidal angle on the separatrix and limiter surfaces

  call calbdpoints(neu_sep_mtheta,neu_sep_r_file,neu_sep_z_file,neu_sep_mtheta_file-1,neu_sep_r,neu_sep_z)
  call calbdpoints(neu_sep_mtheta,lim_org_r,lim_org_z,lim_mdata-1,neu_lim_r,neu_lim_z) 

  ! evaluate separation distances between the separatrix and psi=neu_psi_edge points at theta=neu_theta_edge(k)
  call get_x_r_setup
  if(sml_mype==0)write(111,*)neu_sep_r,neu_sep_z,neu_lim_r,neu_lim_z
  ! for each angle theta, find psi=neu_psi0 (1.03?) position
  ! using binary search
  ! left is magnetic axis and right is ?
  max_r=sml_bd_max_r
  min_r=sml_bd_min_r
  max_z=sml_bd_max_z
  min_z=0.5*(eq_x_z+sml_bd_min_z)
!  min_z=max(min_z, sml_bd_min_z)
  min_z=sml_bd_min_z
   if(sml_mype==0) print *, 'min_r=', min_r, min_z, eq_axis_r,neu_mtheta

  do i=1, neu_mtheta
     theta=real(i-1)*neu_dtheta
!     print *,theta/sml_2pi
     left=0D0
     if(theta<sml_pi) then
        zd=max_z-eq_axis_z
     else
        zd=eq_axis_z-min_z
     endif
     if( theta < 0.5D0*sml_pi .or. theta > 1.5D0*sml_pi ) then
        rd=max_r-eq_axis_r
     else
        rd=eq_axis_r-min_r
     endif
 
     if(neu_grid_mode==0) then   ! original method of constructing neutral birth surface  
!        right=sqrt( (rd*cos(theta))**2 + (zd*sin(theta))**2 )
        tantheta=tan(theta)
        right=sqrt( (1D0 + tantheta**2) / (1D0/rd**2 + tantheta**2/zd**2) )
!        if(sml_mype==0)print *,'R1',right,sml_bd_max_r-1D0,rd
        max_psi=0D0
        do ifactor=50, 100
           factor = ifactor/100D0        
           !linear search
           r_tmp=factor*right*cos(theta)+eq_axis_r
           z_tmp=factor*right*sin(theta)+eq_axis_z
           psi=psi_interpol(r_tmp,z_tmp,0,0)
           if(z_tmp<eq_x_z) then  ! xd_z cut-off
              psi=neu_psi_edge*(1D0 + (z_tmp-eq_x_z)**2/eq_axis_r**2)
!           elseif((z_tmp>eq_x2_z).and.(eq_x2_on==1)) then  ! xd_z cut-off
!              psi=neu_psi_edge*(1D0 + (z_tmp-eq_x2_z)**2/eq_axis_r**2)
           endif 
           if(max_psi<psi) then
              max_psi_r=r_tmp
              max_psi_z=z_tmp
              max_psi=psi
           endif
        enddo
     else    ! improved method of constructing neutral birth surface  
        itheta=int(theta/neu_sep_dtheta+1D-10) + 1
        itheta=min(neu_sep_mtheta, max(1, itheta))
        a2=theta/neu_sep_dtheta - real(itheta -1D0)
        itheta_p1=mod(itheta,neu_sep_mtheta) + 1

        max_psi_r=a2*neu_lim_r(itheta_p1)+(1D0-a2)*neu_lim_r(itheta)
        max_psi_z=a2*neu_lim_z(itheta_p1)+(1D0-a2)*neu_lim_z(itheta)
        max_psi=psi_interpol(max_psi_r,max_psi_z,0,0)  
     endif
!     print *,'R2',r_tmp
     if((max_psi<neu_psi_edge).and.(neu_grid_mode==0)) then
        neu_r0(i)=max_psi_r
        neu_z0(i)=max_psi_z
        if(sml_mype==0) print *, 'max_psi', max_psi_r, max_psi_z
     else
           
        !binary search
        right=sqrt((max_psi_r-eq_axis_r)**2+(max_psi_z-eq_axis_z)**2)
        left=0.3D0*right
        
        do j=1, 30
           mid=0.5D0*(right+left)
           r_tmp=mid*cos(theta)+eq_axis_r
           z_tmp=mid*sin(theta)+eq_axis_z
           psi1=psi_interpol(r_tmp,z_tmp,0,0)
           if(neu_grid_mode>0) then   ! improved method of constructing neutral birth surface  
              psi=neu_psi_interpol(r_tmp,z_tmp,psi1)
!              if((z_tmp>eq_x2_z).and.(eq_x2_on==1)) psi=neu_psi_edge*(1D0+(z_tmp-eq_x2_z)**2)  ! xd_z cut-off  
           else    ! original method of constructing neutral birth surfac
              psi=psi1 
              if(z_tmp<eq_x_z) then  ! xd_z cut-off
                 psi=neu_psi_edge*(1D0+(z_tmp-eq_x_z)**2/eq_axis_r**2)
!              elseif((z_tmp>eq_x2_z).and.(eq_x2_on==1)) then  ! xd_z cut-off
!                 psi=neu_psi_edge*(1D0+(z_tmp-eq_x2_z)**2/eq_axis_r**2)
              endif  
           endif
           if(psi<neu_psi_edge) then
              left=mid
           else
              right=mid
           endif
        enddo
        neu_r0(i)=mid*cos(theta)+eq_axis_r
        neu_z0(i)=mid*sin(theta)+eq_axis_z
        
     endif

  enddo

  do i=1, neu_mtheta
     !jacobian ? --> area (ignore psi direction)
     ip=mod(i+1  -1,neu_mtheta) +1
     im=mod(i-1  -1+neu_mtheta,neu_mtheta) +1     
     neu_jacob(i)=sml_2pi*neu_r0(i) *0.5D0/neu_dtheta*( &
!     neu_jacob(i)=sml_2pi*0.5D0/neu_dtheta*( &
          sqrt( (neu_r0(ip)-neu_r0(i))**2 + (neu_z0(ip)-neu_z0(i))**2 ) + &
          sqrt( (neu_r0(i)-neu_r0(im))**2 + (neu_z0(i)-neu_z0(im))**2 )  )
  enddo
     
  if(sml_mype==0) then
     open(190,file='fort.neu_birthpoints',status='replace')
     do i=1, neu_mtheta
        write(190,*) neu_r0(i),neu_z0(i),neu_jacob(i)
     enddo
     close(190)
  endif

  call get_volume_neu   ! monte-carlo volume calculation

end subroutine neutral2_setup
  
subroutine get_volume_neu
  use ptl_module
  use sml_module
  use eq_module
  use neu_module
  use random_xgc
  implicit none
  integer valid, total,itheta,ipsi,itheta_p1
  real (kind=8) :: rdim, zdim, roffset, zoffset, tvolume, vvolume, &
       r,z,psi,psi1,phi,dum(neu_grid_mpsi,neu_grid_mtheta),dum1(1),dum2(1),theta
  real (kind=8) , external :: psi_interpol, neu_psi_interpol
  real (kind=8) :: a1,a2

!  eq_axis_z=eq_axis_z*0.
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  valid=0
  total=0
  if(sml_mype==0)print *,'rdim==',rdim,sml_bd_min_r,zdim,zoffset,neu_grid_max_psi,eq_x_z,neu_grid_min_psi,neu_grid_dtheta,eq_x_psi
  neu_grid_vol=1D-50   
  do while(valid<neu_monte_num)
     r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) 
     z=zdim*ranx() + zoffset
!     r=rdim*ranx() + roffset  ! flat   --- debug
     psi1=psi_interpol(r,z,0,0)     
     if(neu_grid_mode>0) then
        psi=neu_psi_interpol(r,z,psi1)
     else
        psi=psi1
     endif  
     theta=dacos((r-eq_axis_r)/dsqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
     if(z < eq_axis_z) then
        theta= sml_2pi - theta
     endif
     total=total+1
     if(neu_grid_min_psi < psi .AND. psi < neu_grid_max_psi .and. ((z>eq_x_z) &
        .or.(neu_grid_mode>0))) then 
        valid=valid+1

        ipsi= int((psi-neu_grid_min_psi)/neu_grid_dpsi) + 1 
        ipsi= min(neu_grid_mpsi-1, max(1, ipsi))
        a1= (psi-neu_grid_min_psi)/neu_grid_dpsi -  ipsi +1
        
        itheta= int(theta/neu_grid_dtheta) + 1 
        itheta= min(neu_grid_mtheta, max(1, itheta))
        a2= theta/neu_grid_dtheta - itheta +1
        itheta_p1=mod(itheta +1 - 1,neu_grid_mtheta) + 1
        
        neu_grid_vol(ipsi,itheta)=neu_grid_vol(ipsi,itheta) + (1D0-a1)*(1D0-a2)
        neu_grid_vol(ipsi+1,itheta)=neu_grid_vol(ipsi+1,itheta) + a1*(1D0-a2)
        neu_grid_vol(ipsi,itheta_p1)=neu_grid_vol(ipsi,itheta_p1) + (1D0-a1)*a2
        neu_grid_vol(ipsi+1,itheta_p1)=neu_grid_vol(ipsi+1,itheta_p1) + a1*a2

     endif
  enddo

  ! neu_grid_vol sum
  call my_mpi_allreduce(neu_grid_vol,dum,neu_grid_mpsi*neu_grid_mtheta)
  neu_grid_vol=dum
  ! total sum
  dum1(1)=real(total)
  call my_mpi_allreduce(dum1,dum2,1) !dum2(1) is sum of total


!  sml_marker_den= real(ptl_num)*dum2(1)/ (  rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) * sml_monte_num)  ! normalized unit - NC_norm_r**3 is normalization constant

  ! volume of each flux shell -- shell volume = (ptl number of each shell) / (sum of total) * volume
    neu_grid_vol=neu_grid_vol/dum2(1)*( rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) ) ! bug fix 2002/07/04

end subroutine 

!****************************************************************************
! function for evaluating effective psi for neutral simulation
! this effective psi is the same as poloidal flux inside separatrix, while
! outside the separatrix, it is defined in terms of the radial fraction of
! the distance between separatrix and limiter wall
!****************************************************************************

real (kind=8) function neu_psi_interpol(r,z,psi)
  use eq_module,only: eq_x_psi,eq_x_z,eq_axis_r,eq_axis_z
  use neu_module
  use sml_module,only: sml_2pi,sml_pi
  implicit none
  real(kind=8) , intent(in) :: r,z,psi
  real(kind=8) :: theta,r_sep,a2,z_sep,rad_dist,theta_dist
  integer :: itheta, itheta_p1

!  eq_axis_z=eq_axis_z*0.
  neu_psi_interpol=psi
  theta=atan2(z-eq_axis_z,r-eq_axis_r); if(theta<0D0) theta=theta+sml_2pi

  if((theta>neu_theta_edge(1)).and.(theta<neu_theta_edge(2)).and. &
     ((psi>eq_x_psi).or.((psi<eq_x_psi).and.(z<eq_x_z)))) then
     itheta=int(theta/neu_sep_dtheta) + 1
     itheta=min(neu_sep_mtheta, max(1, itheta))
     a2= theta/neu_sep_dtheta - itheta +1
     itheta_p1=mod(itheta,neu_sep_mtheta) + 1

     r_sep=a2*neu_sep_r(itheta_p1)+(1D0-a2)*neu_sep_r(itheta)
     z_sep=a2*neu_sep_z(itheta_p1)+(1D0-a2)*neu_sep_z(itheta)

     theta_dist=neu_theta_edge(2)-neu_theta_edge(1)
     rad_dist=(theta-neu_theta_edge(1))/theta_dist*neu_rad_edge(2)+(neu_theta_edge(2)-theta)/theta_dist*neu_rad_edge(1)
     neu_psi_interpol=eq_x_psi+dsqrt((r-r_sep)**2+(z-z_sep)**2)/rad_dist*(neu_psi_edge-eq_x_psi)
  endif

!  if(eq_x2_on==1) then
!     if((theta>neu_theta_edge2(1)).and.(theta<neu_theta_edge2(2)).and. &
!        ((psi>eq_x_psi).or.((psi<eq_x2_psi).and.(z>eq_x2_z)))) then
!        itheta=int(theta/neu_sep_dtheta) + 1
!        itheta=min(neu_sep_mtheta, max(1, itheta))
!        a2= theta/neu_sep_dtheta - itheta +1
!        itheta_p1=mod(itheta,neu_sep_mtheta) + 1

!        r_sep=a2*neu_sep_r(itheta_p1)+(1D0-a2)*neu_sep_r(itheta)
!        z_sep=a2*neu_sep_z(itheta_p1)+(1D0-a2)*neu_sep_z(itheta)

!        theta_dist=neu_theta_edge2(2)-neu_theta_edge2(1)
!        rad_dist=(theta-neu_theta_edge2(1))/theta_dist*neu_rad_edge2(2)+(neu_theta_edge2(2)-theta)/theta_dist*neu_rad_edge2(1)
!        neu_psi_interpol=eq_x_psi+dsqrt((r-r_sep)**2+(z-z_sep)**2)/rad_dist*(neu_psi_edge-eq_x_psi)
!     endif
!  endif
end function neu_psi_interpol

!*********************************************************************************************************
! evaluate radial separation between the separatrix and psi=neu_psi_edge points at theta=neu_theta_edge(k)
!*********************************************************************************************************

subroutine get_x_r_setup
  use eq_module
  use neu_module
  implicit none
  real (kind=8) :: theta,r,z,r_sep,z_sep,r_lim,z_lim,a2,mid,left,right,psi,r_tmp,z_tmp
  integer :: j,k,itheta,itheta_p1
  real(kind=8), external :: psi_interpol

!  eq_axis_z=eq_axis_z*0.
  do k=1, 2
     theta=neu_theta_edge(k)
     itheta=int(theta/neu_sep_dtheta+1D-10) + 1
     itheta=min(neu_sep_mtheta, max(1, itheta))
     a2=theta/neu_sep_dtheta - real(itheta -1D0)
     itheta_p1=mod(itheta,neu_sep_mtheta) + 1
     call check_point('inside set_x_r_setup')
     r_sep=a2*neu_sep_r(itheta_p1)+(1D0-a2)*neu_sep_r(itheta)
     z_sep=a2*neu_sep_z(itheta_p1)+(1D0-a2)*neu_sep_z(itheta)
     r_lim=a2*neu_lim_r(itheta_p1)+(1D0-a2)*neu_lim_r(itheta)
     z_lim=a2*neu_lim_z(itheta_p1)+(1D0-a2)*neu_lim_z(itheta)

     ! binary search
     right=sqrt((r_lim-eq_axis_r)**2+(z_lim-eq_axis_z)**2)
     left=sqrt((r_sep-eq_axis_r)**2+(z_sep-eq_axis_z)**2)
     do j=1, 30
        mid=0.5D0*(right+left)
        r_tmp=mid*cos(theta)+eq_axis_r
        z_tmp=mid*sin(theta)
        psi=psi_interpol(r_tmp,z_tmp,0,0)
        if(psi<neu_psi_edge) then
           left=mid
        else
           right=mid
        endif
     enddo
     r=mid*cos(theta)+eq_axis_r
     z=mid*sin(theta)
     neu_rad_edge(k)=dsqrt((r-r_sep)**2+(z-z_sep)**2)
  enddo
 
!  if(eq_x2_on==1) then
!     do k=1, 2
!        theta=neu_theta_edge2(k)
!        itheta=int(theta/neu_sep_dtheta+1D-10) + 1
!        itheta=min(neu_sep_mtheta, max(1, itheta))
!        a2=theta/neu_sep_dtheta - real(itheta -1D0)
!        itheta_p1=mod(itheta,neu_sep_mtheta) + 1

!        r_sep=a2*neu_sep_r(itheta_p1)+(1D0-a2)*neu_sep_r(itheta)
!        z_sep=a2*neu_sep_z(itheta_p1)+(1D0-a2)*neu_sep_z(itheta)
!        r_lim=a2*neu_lim_r(itheta_p1)+(1D0-a2)*neu_lim_r(itheta)
!        z_lim=a2*neu_lim_z(itheta_p1)+(1D0-a2)*neu_lim_z(itheta)

        ! binary search
!        right=sqrt((r_lim-eq_axis_r)**2+(z_lim-eq_axis_z)**2)
!        left=sqrt((r_sep-eq_axis_r)**2+(z_sep-eq_axis_z)**2)
!        do j=1, 30
!           mid=0.5D0*(right+left)
!           r_tmp=mid*cos(theta)+eq_axis_r
!           z_tmp=mid*sin(theta)
!           psi=psi_interpol(r_tmp,z_tmp,0,0)
!           if(psi<neu_psi_edge) then
!              left=mid
!           else
!              right=mid
!           endif
!        enddo
!        r=mid*cos(theta)+eq_axis_r
!        z=mid*sin(theta)
!        neu_rad_edge2(k)=dsqrt((r-r_sep)**2+(z-z_sep)**2)
!     enddo
!  endif
end subroutine get_x_r_setup
