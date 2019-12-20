subroutine efield(grid,psn,sp,i,fld,time)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use eq_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  integer, intent(in) :: i
  type(species_type) :: sp
  real (kind=8) :: time

  
  if(sml_00_efield .or. sml_turb_efield) then
     !Gyrokinetic E
     call efield_gk(grid,psn,sp,i,fld)
  else
       fld%Er=0D0
       fld%Ez=0D0
       fld%Ephi=0D0
       fld%Epot=0D0
       fld%Er00=0D0
       fld%Ez00=0D0
       fld%ddpotdt=0D0
  endif

end subroutine efield


subroutine efield_gk(grid,psn,sp,i,fld) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer, intent(in) :: i
  type(fld_type), intent(inout) :: fld

  integer :: itr,ip,node,irho, iphi
  real (kind=8) :: wphi(0:1), wp, rho, rhon, wrho(2)
  real (kind=8) :: pot, E(3),E00(2), B , ddpotdt !, D(3)
  real (kind=8) :: E_phi_ff(3,0:1)
  iphi=floor(fld%phi/grid%delta_phi)
  !wphi(1)=(fld%phi/grid%delta_phi)  - grid%iphi_offset
  wphi(1)=(fld%phi/grid%delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)



  pot=0D0
  E=0D0
  E00=0D0
  !D=0D0 !debug
  ddpotdt=0D0

  !get E
  itr=sp%tr_save(i)
  if(itr>0) then        
     do ip=1, 3
        node=grid%nd(ip,itr)
        wp=sp%p_save(ip,i)
        
        ! find gyro-radius
        if(sp%type==1) then
           rho=sp%rhoi(i)  ! for ion

           !find rho index

           rhon=min(rho,grid%rhomax)/grid%drho
           irho=min(floor(rhon),grid%nrho-1)
           wrho(2)=rhon - real(irho)
           wrho(1)=1D0-wrho(2)
           
           pot = pot + wp*wphi(0)*wrho(1)*psn%pot_rho_ff(0,irho  ,node) 
           pot = pot + wp*wphi(1)*wrho(1)*psn%pot_rho_ff(1,irho  ,node) 
           pot = pot + wp*wphi(0)*wrho(2)*psn%pot_rho_ff(0,irho+1,node) 
           pot = pot + wp*wphi(1)*wrho(2)*psn%pot_rho_ff(1,irho+1,node)
           ! following 4 lines are not working correclty with -fast option - maybe compiler bug
           E   = E   + wp*wphi(0)*wrho(1)*psn%E_rho_ff(:,0,irho  ,node)
           E   = E   + wp*wphi(1)*wrho(1)*psn%E_rho_ff(:,1,irho  ,node)
           E   = E   + wp*wphi(0)*wrho(2)*psn%E_rho_ff(:,0,irho+1,node)
           E   = E   + wp*wphi(1)*wrho(2)*psn%E_rho_ff(:,1,irho+1,node)
        else
           ! for electron -- rho=0 case, but iphi=/=0.           
           pot = pot + wp*wphi(0)*psn%pot_phi_ff(0,node,iphi) 
           pot = pot + wp*wphi(1)*psn%pot_phi_ff(1,node,iphi)

#ifdef USE_CALC_GRADIENT
           call calc_E_phi_ff(grid,psn,node,iphi,E_phi_ff )
           E = E + wp*wphi(0)*E_phi_ff(:,0)
           E = E + wp*wphi(1)*E_phi_ff(:,1)
#else
           E   = E   + wp*wphi(0)*psn%E_phi_ff(:,0,node,iphi)
           E   = E   + wp*wphi(1)*psn%E_phi_ff(:,1,node,iphi)
#endif


           E00 = E00 + wp*wphi(0)*psn%E00_ff(:,0,node)
           E00 = E00 + wp*wphi(1)*psn%E00_ff(:,1,node)

           ddpotdt = ddpotdt + wp*wphi(0)*psn%ddpotdt_phi(node,0,iphi)
           ddpotdt = ddpotdt + wp*wphi(1)*psn%ddpotdt_phi(node,1,iphi)
        endif


!!$        E(1)   = E(1)   + wp*wphi(0)*wrho(1)*psn%E_rho_ff(1,0,irho  ,node)
!!$        E(2)   = E(2)   + wp*wphi(0)*wrho(1)*psn%E_rho_ff(2,0,irho  ,node)
!!$        E(3)   = E(3)   + wp*wphi(0)*wrho(1)*psn%E_rho_ff(3,0,irho  ,node)
!!$
!!$        E(1)   = E(1)   + wp*wphi(1)*wrho(1)*psn%E_rho_ff(1,1,irho  ,node)
!!$        E(2)   = E(2)   + wp*wphi(1)*wrho(1)*psn%E_rho_ff(2,1,irho  ,node)
!!$        E(3)   = E(3)   + wp*wphi(1)*wrho(1)*psn%E_rho_ff(3,1,irho  ,node)
!!$
!!$        E(1)   = E(1)   + wp*wphi(0)*wrho(2)*psn%E_rho_ff(1,0,irho+1,node)
!!$        E(2)   = E(2)   + wp*wphi(0)*wrho(2)*psn%E_rho_ff(2,0,irho+1,node)
!!$        E(3)   = E(3)   + wp*wphi(0)*wrho(2)*psn%E_rho_ff(3,0,irho+1,node)
!!$
!!$        E(1)   = E(1)   + wp*wphi(1)*wrho(2)*psn%E_rho_ff(1,1,irho+1,node)
!!$        E(2)   = E(2)   + wp*wphi(1)*wrho(2)*psn%E_rho_ff(2,1,irho+1,node)
!!$        E(3)   = E(3)   + wp*wphi(1)*wrho(2)*psn%E_rho_ff(3,1,irho+1,node)



!        pot = pot + wp*wphi(0)*wrho(1)*(psn%pot0(node)+psn%dpot(node,0))
!        pot = pot + wp*wphi(1)*wrho(1)*(psn%pot0(node)+psn%dpot(node,1))
!        pot = pot + wp*wphi(0)*wrho(2)*(psn%pot0(node)+psn%dpot(node,0))
!        pot = pot + wp*wphi(1)*wrho(2)*(psn%pot0(node)+psn%dpot(node,1))
!        D(1:2)=D(1:2)  + wp*wphi(0)*wrho(1)*psn%E_perp_node(:,node,0)
!        D(1:2)=D(1:2)  + wp*wphi(1)*wrho(1)*psn%E_perp_node(:,node,1)
!        D(1:2)=D(1:2)  + wp*wphi(0)*wrho(2)*psn%E_perp_node(:,node,0)
!        D(1:2)=D(1:2)  + wp*wphi(1)*wrho(2)*psn%E_perp_node(:,node,1)
!        D(3)=D(3)  + wp*wphi(0)*wrho(1)*psn%E_para(node,0)
!        D(3)=D(3)  + wp*wphi(1)*wrho(1)*psn%E_para(node,1)
!        D(3)=D(3)  + wp*wphi(0)*wrho(2)*psn%E_para(node,0)
!        D(3)=D(3)  + wp*wphi(1)*wrho(2)*psn%E_para(node,1)


     enddo
  else
     print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
  endif
    
  !E(3) was para E-field and becomes Ephi
  if(sml_turb_efield) then
     B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     E(3)=(E(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi

     !debug
     !D(3)=(D(3)*B- D(1)*fld%Br - D(2)*fld%Bz)/fld%Bphi
  else
     E(3)=0D0
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
  fld%Epot=pot
  fld%Er00=E00(1)
  fld%Ez00=E00(2) 
  fld%ddpotdt=ddpotdt
end subroutine efield_gk

subroutine efield_init(grid,psn)
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  
  call get_potential_grad(grid,psn)

end subroutine efield_init




