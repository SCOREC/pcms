attributes(device) &
subroutine sheath_calculation_gpu(iptl,epc,type_gpu,itrout,pout,ith,tb)

#ifdef USE_GPU_EMU
#define atomicAdd atomicAdd_d
#endif

  use grid_class_gpu
#ifdef USE_TEXTURE
  use psn_class_gpu, only :                           &
     psn_nwall,                                       &
     psn_wall_nodes,                                  &
     psn_sheath_lost,                                 &
     psn_pot0 => t_psn_pot0
#else
  use psn_class_gpu, only :                           &
     psn_nwall,                                       &
     psn_wall_nodes,                                  &
     psn_sheath_lost,                                 &
     psn_pot0
#endif

  use sml_module_gpu
  use ptl_module_gpu
!  use diag_module_gpu
  use precision_mod_gpu
  use neu_module_gpu, only : neu_weight_sum_lost_gpu

  implicit none

  integer, intent(in) :: type_gpu, iptl, epc, ith
  type(tbuf), intent(inout) :: tb
  integer :: itrout
  integer :: i,l, node,itr
  real (kind=work_p) :: pout(3)
  real (kind=work_p) :: p(3),psave(3)
  integer, parameter :: rgn_wall=100
  real (kind=work_p), parameter :: minus_val=-1D50
  real (kind=work_p) :: rho,b,en_para, x(2),phi, phi_mid, xff(2)
  real (kind=work_p) :: time_now, xn(2), dist_sqr, dist_min
  integer :: node_min
  real (kind=work_p)  :: b_interpol_gpu

  real (kind=work_p) :: new_phase(ptl_nphase), dummy
  integer :: widx, ip
  real (kind=work_p)  :: w1_change, en_perp
  real (kind=work_p)  :: tmp_ph(ptl_nphase), tmp_ct(ptl_nconst)
  ! find nearest wall point

  new_phase = phase0_gpu(iptl,:)
  x = new_phase(1:2)
  phi=new_phase(3)
  phi_mid=(floor(phi/grid_delta_phi) + 0.5_work_p) * grid_delta_phi

  ! get field following posision at 1/2 angle
  call field_following_pos2_gpu(x,phi,phi_mid,xff)

  ! find position of previous time step

     call search_tr2_gpu(xff,itr,p)
 
 
  psave=p

  ! if old position is also outside of grid --> remove particle and return
  if(itr<0) then
     call remove_particle_gpu(iptl,-1,tb)
     itrout=itr
     pout=p
     return
  endif

  ! search three nodes of the triangle and check if it is wall nodes
  do i=1, 3
!     l=maxloc(p,1)

     l = 3
     if ((p(2) >= p(1)).and.(p(2) >= p(3))) l = 2
     if ((p(1) >= p(2)).and.(p(1) >= p(3))) l = 1

     node = grid_nd(l,itr)

     if(grid_rgn(node)==rgn_wall) then
        exit
     else
        p(l)=minus_val
     endif
  enddo

  !if the triangle does not have a wall node
  ! find nearest one using primitive scan
  if(grid_rgn(node)/=rgn_wall) then
     dist_min = 1D99
     do i=1, psn_nwall ! for all wall nodes
        ! check distance
        xn=grid_x(:,psn_wall_nodes(i))
        dist_sqr=(xn(1)-xff(1))**2 + (xn(2)-xff(2))**2
        ! check minimum
        if(dist_min > dist_sqr) then
           dist_min=dist_sqr
           node_min = i
        endif
     enddo
     node=psn_wall_nodes(node_min)
  endif


  !
  !check potential and energy

  rho=new_phase(4)
  b = b_interpol_gpu(x(1),x(2),0.0_work_p)

  en_para=ptl_c2_2m(type_gpu) * (rho*b)**2

  ! reflection
  new_phase(pirho)=-new_phase(pirho)

  widx=psn_node_to_wall(node)
  if(en_para < psn_sheath_pot(widx) *sml_ev2j .and. type_gpu==0) then
     ! do nothing -- just relection due to sheath potential     
  else
     w1_change = 1.0_work_p - new_phase(piw2) + new_phase(piw1)

     ! sum-up to wall node

!     psn_sheath_lost(widx)=psn_sheath_lost(widx) + (1D0 +  new_phase(piw2))*sp%charge*sp%ptl(iptl)%ct(piw0)
     if(sml_ipc==2 .and. epc==1) dummy = atomicAdd(psn_sheath_lost(widx), (w1_change)*ptl_charge(type_gpu)*tb%ct(piw0))
     ! --- reflect with weight change
     new_phase(piw1)=-1.0_work_p+new_phase(piw2)
     ! w2 does not change

     ! for neutral
     if(sml_neutral) then
        if(sml_ipc==2) then
          if(epc==1) then
            ip = mod( iptl, size(neu_weight_sum_lost_gpu)) + 1 
            dummy = atomicAdd( neu_weight_sum_lost_gpu(ip), w1_change*tb%ct(piw0))
!           neu_weight_sum_lost(ith) = neu_weight_sum_lost(ith) + (w1_change)*sp%ptl(iptl)%ct(piw0)
          endif
        endif
     endif

     ! for heat load diagnosie
     if(sml_ipc==2 .and. diag_heat_on .and. epc==1) then
!    if(ipc_gpu==2)then 
        en_perp = tb%ct(pim)*b
        ! arguments:
        ! 1. weight_change, 2. potential, 3. en_para, 4. en_perp, 5. ct, 6. phase(old), 7. phase(new), 8. stype
        ! new_phase is old position
        tmp_ph = tb%ph(:)
        tmp_ct = tb%ct(:)
!        print *, "call diag_heat_port_gpu"
        call diag_heat_port_gpu(w1_change,psn_sheath_pot(widx),en_para, en_perp, tmp_ct, new_phase, tmp_ph,grid_delta_phi, type_gpu,ith)
     endif

  endif

  tb%ph(:)=new_phase

  if(ptl_gid_gpu(iptl)>0) then
     itrout=itr
     pout=psave
  endif

end subroutine sheath_calculation_gpu

