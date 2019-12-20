attributes(device) &
subroutine efield_gk_gpu(i,fld,itr,p) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module_gpu, only : sml_turb_efield, sml_mype, sml_turb_efield
  use grid_class_gpu, only : grid_delta_phi, grid_nd, grid_rhomax, grid_drho, grid_nrho                             

  use psn_class_gpu, only :                &
     psn_pot_rho_ff,                       &
     psn_E_rho_ff,                         &
     psn_ddpotdt,                          &
     psn_E00_ff,                           &
     psn_pot_phi_ff,                       &
     psn_E_phi_ff0,                         &
     psn_E_phi_ff1,                         &     
     psn_ddpotdt_phi

  use fld_module, only : fld_type
  use ptl_module_gpu, only :  type_gpu, rhoi_gpu
  use precision_mod_gpu
  implicit none
  integer, intent(in) :: i,itr
  type(fld_type), intent(inout) :: fld

  integer :: ip,node,irho, iphi
  real (kind=work_p) :: wphi(0:1), wp, rho, rhon, wrho(2), p(3)
  real (kind=work_p) :: pot, E(3),E00(2), B , ddpotdt !, D(3)

#ifdef USE_CALC_GRADIENT
  real (kind=work_p) :: E_phi_ff
#endif


  iphi=floor(fld%phi/grid_delta_phi)
!  wphi(1)=(fld%phi/grid%delta_phi)  - grid%iphi_offset
  wphi(1)=(fld%phi/grid_delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)

  pot=0D0
  E=0D0
  E00=0D0
  !D=0D0 !debug
  ddpotdt=0D0

  !get E
 
  if(itr>0) then        
     do ip=1, 3
        node=grid_nd(ip,itr)
        wp=p(ip)
        
        ! find gyro-radius
!        if(type_gpu==1) then
!           rho=rhoi_gpu(i)  ! for ion

           !find rho index

!           rhon=min(rho,grid_rhomax)/grid_drho
!           irho=min(floor(rhon),grid_nrho-1)
!           wrho(2)=rhon - real(irho)
!           wrho(1)=1D0-wrho(2)
           
!           pot = pot + wp*wphi(0)*wrho(1)*psn_pot_rho_ff(0,irho  ,node) 
!           pot = pot + wp*wphi(1)*wrho(1)*psn_pot_rho_ff(1,irho  ,node) 
!           pot = pot + wp*wphi(0)*wrho(2)*psn_pot_rho_ff(0,irho+1,node) 
!           pot = pot + wp*wphi(1)*wrho(2)*psn_pot_rho_ff(1,irho+1,node)
!           ! following 4 lines are not working correclty with -fast option - maybe compiler bug
!           E   = E   + wp*wphi(0)*wrho(1)*psn_E_rho_ff(:,0,irho  ,node)
!           E   = E   + wp*wphi(1)*wrho(1)*psn_E_rho_ff(:,1,irho  ,node)
!           E   = E   + wp*wphi(0)*wrho(2)*psn_E_rho_ff(:,0,irho+1,node)
!           E   = E   + wp*wphi(1)*wrho(2)*psn_E_rho_ff(:,1,irho+1,node)
!        else
           ! for electron -- rho=0 case, optimized.           

           pot = pot + wp*wphi(0)*psn_pot_phi_ff(0,node,iphi) 
           pot = pot + wp*wphi(1)*psn_pot_phi_ff(1,node,iphi)


#ifdef USE_CALC_GRADIENT
           call calc_E_phi_ff_gpu(node,iphi,E_phi_ff )
           E = E + wp*wphi(0)*E_phi_ff(:,0)
           E = E + wp*wphi(1)*E_phi_ff(:,1)
#else
           E   = E   + wp*wphi(0)*psn_E_phi_ff0(:,node,iphi)
           E   = E   + wp*wphi(1)*psn_E_phi_ff1(:,node,iphi)
#endif

           E00 = E00 + wp*wphi(0)*psn_E00_ff(:,0,node)
           E00 = E00 + wp*wphi(1)*psn_E00_ff(:,1,node)

           ddpotdt = ddpotdt + wp*wphi(0)*psn_ddpotdt_phi(node,0,iphi)
           ddpotdt = ddpotdt + wp*wphi(1)*psn_ddpotdt_phi(node,1,iphi)
!        endif


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
!      print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
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


end subroutine efield_gk_gpu
