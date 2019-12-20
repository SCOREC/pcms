#define XGC1 1
!***************************************************************************
!*** Poisson: physics
!***************************************************************************
subroutine read_add_pot0(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: num, i

  open(unit=18,file=sml_add_pot0_file,action='read')
  read(18,*) num
  if(num/=grid%nnode) then
     if(sml_mype==0) print *, &
          'Error : Grid number missmatch in read_add_pot0',num, grid%nnode
     stop
  endif
  do i=1, num
     read(18,*) psn%add_pot0(i)
  enddo

  close(18)
end subroutine read_add_pot0

subroutine get_add_pot0_psi(psi_in,pot,grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8), intent(in) :: psi_in
  real (kind=8), intent(out) :: pot
  real (kind=8) :: psi,x(2), p(3)
  integer :: itr, nodes(3), j
  logical, parameter :: USE_SEARCH_TR2 = .true.
  real (kind=8) , external :: get_mid_r

  ! exit with 0
  if(sml_add_pot0<=0) then
     pot=0D0
     return
  endif
  ! set psi
  psi=min(max(psi_in,sml_inpsi),sml_outpsi)

  ! get r,z=0 from psi
#ifdef XGC1
  x(1)=get_mid_R(psi)
#else
  x(1)=get_mid_R(psi,1)
#endif
  x(2)=eq_axis_z

  ! find triangle

  ! find position and save it.
  if (USE_SEARCH_TR2) then
#ifdef XGC1
     call search_tr2(grid,x,itr,p)
#else
     call search_tr2(grid,x,psi,itr,p)
#endif
  else
     call search_tr(grid,x,itr,p)
  endif


  ! get potential value
  pot=0D0
  if( itr > 0 ) then

     nodes=grid%nd(:,itr)

     do j=1, 3
        pot= pot + p(j) * psn%add_pot0(nodes(j))
     enddo
  else
     print *, 'Warning : itr <0 in get_add_pot0_psi', psi, psi_in
     call err_count
  endif


end subroutine get_add_pot0_psi


subroutine zero_out_total_charge(grid,psn,den_org,den_zero)
  use grid_class
  use psn_class
  implicit none
!  integer :: flag
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: den_org(grid%nnode), den_zero(grid%nnode), total,area
  integer :: ibd, obd

!!$  ibd=psn%pbd0%in%end+1
!!$  obd=psn%pbd0%out1%start-1
!!$
!!$  ! Volume?? Area??
!!$  total=sum( den_org(ibd:obd)*grid%node_area(ibd:obd) )
!!$  area=sum(grid%node_area(ibd:obd))
!!$
!!$  den_zero=den_org - total/area

end subroutine zero_out_total_charge


subroutine set_decaying_boundary(den,grid)
  use grid_class
  use sml_module
  implicit none
!  integer :: flag
  type(grid_type) :: grid
  real (kind=8) :: den(grid%nnode)
  real (kind=8) :: sml_f0_psi_c, sml_f0_1_psi_w, alpha
  integer :: i

  sml_f0_psi_c=0.5*(sml_inpsi+sml_outpsi)
  sml_f0_1_psi_w=1D0/( 0.35*(sml_outpsi-sml_inpsi) )

  do i=1, grid%nnode
    alpha= exp( - ((grid%psi(i) - sml_f0_psi_c)*sml_f0_1_psi_w )**6 )
    den(i)=den(i)* alpha
  enddo

end subroutine


!! return gyro radius with a give mu and position
real (kind=8) function gyro_radius(x,mu)
  use ptl_module
  use eq_module
  implicit none
  real (kind=8),intent(in) :: x(2),mu
  real (kind=8) :: B
  real (kind=8),external :: b_interpol


  B=b_interpol(x(1),x(2),0D0)
#ifdef XGC_FLAT_B
  b=eq_axis_b
#endif
  gyro_radius=sqrt(mu/B/ptl_c2_2m(1))


end function gyro_radius

subroutine gyro_radius_vec(nx,x,mu,gyro_radius)
  use ptl_module
  implicit none
  integer, intent(in) :: nx
  real (kind=8),intent(in) :: x(2,nx),mu
  real (kind=8), intent(inout) :: gyro_radius(nx)

  real (kind=8) :: B(nx), phi(nx)

  interface
     subroutine b_interpol_vec(n,r,z,phi,b_interpol)
       implicit none
       integer, intent(in) :: n
       real(kind=8), dimension(n), intent(in) :: r,z,phi
       real(kind=8), dimension(n), intent(inout) :: b_interpol
     end subroutine b_interpol_vec
  end interface


  phi = 0.0d0
!  call b_interpol_vec(nx, x(1,1:nx), x(2,1:nx), phi, B )
  !  B(1:nx)=b_interpol_vec(x(1,1:nx),x(2,1:nx),0D0)
  gyro_radius=sqrt(mu/(B*ptl_c2_2m(1)))


end subroutine  gyro_radius_vec



!!return gyro radius of a thermal ion at a given position
real (kind=8) function gyro_radius2(x)
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  real (kind=8),intent(in) :: x(2)
  real (kind=8) :: B,psi,tev
  real (kind=8),external :: b_interpol,psi_interpol
!  real (kind=8),external :: init_tempi_ev

  B=b_interpol(x(1),x(2),0D0)
  psi=psi_interpol(x(1),x(2),0,0)
  tev=eq_ftn(psi,0D0,0D0,eq_tempi)

  gyro_radius2=sqrt(Tev*sml_ev2j/ptl_c2_2m(1)/2D0)/B


end function gyro_radius2

real (kind=8) function gyro2_tev(x)
  use sml_module
  use ptl_module
  implicit none
  real (kind=8),intent(in) :: x(2)
  real (kind=8) :: B,psi,tev
  real (kind=8),external :: b_interpol,psi_interpol
!  real (kind=8),external :: init_tempi_ev

  B=b_interpol(x(1),x(2),0D0)
  psi=psi_interpol(x(1),x(2),0,0)
!  tev=init_tempi_ev(psi)

!  gyro2_tev=sml_ev2j/ptl_c2_2m(1)/2D0/B**2
  gyro2_tev=sml_ev2j/(ptl_c2_2m(1)*2D0*B**2)

end function gyro2_tev

!! Write out matrix element for poisson equation
!subroutine output_matrix(grid,psn)
!  use psn_class
!  use grid_class
!  implicit none
!  integer :: i,j
!  type(psn_type) :: psn
!  type(grid_type) :: grid
!
!  open(unit=201,file='poisson.matrix',status='replace')
!  do i=1, grid%nnode
!     do j=1, psn%Amat%nelement(i)
!        write(201,*) i,psn%Amat%eindex(j,i),psn%Amat%value(j,i),psn%Amat%value(j,i)/grid%node_area(i)
!     enddo
!  enddo
!  close(201)
!
!end subroutine output_matrix


#ifndef TRIANGLE_GYRO_AVG
#ifdef OPTIM_GYRO_AVG_MAT
subroutine init_gyro_avg_mat(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  integer :: iphi, irho
  real (8) :: rho, phi0, phi
  integer :: i, larmor
  real (8) :: x(2), xff(2), x_ring(2), x_rff(2)
  integer, parameter :: n_gyro=32 !8
  real (8) :: dx_unit(2,n_gyro), dx_unit2(2)
  real (8) :: angle, dpdr, dpdz, dp, cosa, sina
  integer :: itr, node, ip
  real (8) :: p(3), weight
  !
  real (8), external :: psi_interpol
  integer :: loop_action,start_loop

  start_loop=2

  if(start_loop.eq.2)then
    call new_mat(psn%gyro_avg_mat,grid%nnode,n_gyro*3)
  elseif(start_loop.eq.1)then
    call new_mat(psn%gyro_avg_mat,grid%nnode,0)
  endif

#ifdef XGC1
  iphi=mod(sml_plane_mype,2)
  irho=sml_plane_mype/2  + 1

  phi0=0.5D0*grid%delta_phi
  phi=real(iphi)*grid%delta_phi  ! 0 for 0 , 1 for delta_phi
#else
  iphi=0
  phi=0D0
  irho=sml_plane_mype + 1
#endif

  rho=grid%drho*real(irho)

  do larmor=1, n_gyro
     angle=sml_2pi/real(n_gyro)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  do loop_action=start_loop,2

  do i=1, grid%nnode
     !find real position -- follow field line
#ifdef XGC1
     xff=grid%x(:,i)
     call field_following_pos2(xff,phi0,phi,x)
#else
     x=grid%x(:,i)
#endif
     ! find gyro-ring position
     dpdr=psi_interpol(x(1),x(2),1,0)  !### get together 3 of them (psi, dpdr, dpdz)
     dpdz=psi_interpol(x(1),x(2),0,1)
     dp=sqrt(dpdr**2 + dpdz**2)
     cosa=dpdr/dp
     sina=dpdz/dp


     do larmor=1, n_gyro

        !find 8-point with rhoi
        dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
        dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa

        x_ring = x + rho* dx_unit2(:)

#ifdef XGC1
        ! back transform to ff position
        call field_following_pos2(x_ring,phi,phi0,x_rff)
        ! find triangle
        call search_tr2(grid,x_rff,itr,p)
#else
        x_rff=x_ring
        call search_tr2(grid,x_rff,-1D10,itr,p)
#endif

        if(itr>0) then
           do ip=1,3
              weight=p(ip)/real(n_gyro)
              node=grid%nd(ip,itr)
              if(loop_action.eq.2)then
                call set_value(psn%gyro_avg_mat,i,node,weight,1)
              elseif(loop_action.eq.1)then
                call set_value_width(psn%gyro_avg_mat,i,node)
              endif
           enddo
        endif
     enddo
     ! column normalization is required for boundary nodes
  enddo

  if(loop_action.eq.1) call new_value(psn%gyro_avg_mat) ! Init the value array with varying width 

  enddo

end subroutine init_gyro_avg_mat

#else
subroutine init_gyro_avg_mat(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  use ptl_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  integer :: iphi, irho
  real (8) :: rho, phi0, phi
  integer :: i, larmor
  real (8) :: x(2), xff(2), x_ring(2), x_rff(2)
  integer, parameter :: n_gyro=32 !8
  real (8) :: dx_unit(2,n_gyro), dx_unit2(2)
  real (8) :: angle, dpdr, dpdz, dp, cosa, sina
  integer :: itr, node, ip
  real (8) :: p(3), weight
  !
  real (8), external :: psi_interpol
  call new_mat(psn%gyro_avg_mat,grid%nnode,n_gyro*3)

#ifdef XGC1
  iphi=mod(sml_plane_mype,2)
  irho=sml_plane_mype/2  + 1

  phi0=0.5D0*grid%delta_phi
  phi=real(iphi)*grid%delta_phi  ! 0 for 0 , 1 for delta_phi
#else
  iphi=0
  phi=0D0
  irho=sml_plane_mype + 1
#endif

  rho=grid%drho*real(irho)

  do larmor=1, n_gyro
     angle=sml_2pi/real(n_gyro)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  do i=1, grid%nnode
     !find real position -- follow field line
#ifdef XGC1
     xff=grid%x(:,i)
     call field_following_pos2(xff,phi0,phi,x)
#else
     x=grid%x(:,i)
#endif


     ! find gyro-ring position
     dpdr=psi_interpol(x(1),x(2),1,0)  !### get together 3 of them (psi, dpdr, dpdz)
     dpdz=psi_interpol(x(1),x(2),0,1)
     dp=sqrt(dpdr**2 + dpdz**2)
     cosa=dpdr/dp
     sina=dpdz/dp


     do larmor=1, n_gyro

        !find 8-point with rhoi
        dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
        dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa

        x_ring = x + rho* dx_unit2(:)

#ifdef XGC1
        ! back transform to ff position
        call field_following_pos2(x_ring,phi,phi0,x_rff)
        ! find triangle
        call search_tr2(grid,x_rff,itr,p)
#else
        x_rff=x_ring
        call search_tr2(grid,x_rff,-1D10,itr,p)
#endif

        if(itr>0) then
           do ip=1,3
              weight=p(ip)/real(n_gyro)
              node=grid%nd(ip,itr)
              call set_value(psn%gyro_avg_mat,i,node,weight,1)
           enddo
        endif
     enddo

     ! column normalization is required for boundary nodes
  enddo

end subroutine init_gyro_avg_mat
#endif
#else

subroutine init_gyro_avg_mat(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  integer :: iphi, irho
  real (8) :: rho, phi0, phi
  integer, parameter :: n_gyro=32 !32, nsub=3
  integer :: larmor, i, nd(3), j, k, itmp
  real (8) :: dx_unit(2,n_gyro), dx_unit2(2)
  real (8) :: dx1(2), dx2(2), area, c1, c2, c3
  real (8) :: x(2), xff(2), x_ring(2), x_rff(2)
  real (8) :: angle, dpdr, dpdz, dp, cosa, sina
  integer :: itr, node, ip
  real (8) :: p(3), weight
  !
  real (8), external :: psi_interpol
  !

  ! matrix size??  factor 2 at the end is approximate factor
  call new_mat(psn%gyro_avg_mat,grid%nnode,n_gyro*3*nsub*8)

#ifdef XGC1
  iphi=mod(sml_plane_mype,2)
  irho=sml_plane_mype/2  + 1

  phi0=0.5D0*grid%delta_phi
  phi=real(iphi)*grid%delta_phi  ! 0 for 0 , 1 for delta_phi
#else
  iphi=0
  phi=0D0
  irho=sml_plane_mype + 1
#endif

  rho=grid%drho*real(irho)


  do larmor=1, n_gyro
     angle=sml_2pi/real(n_gyro)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo


  do i=1, grid%ntriangle ! for all triangle

     nd=grid%nd(:,i)
     dx1= ( grid%x(:,nd(1)) - grid%x(:,nd(3)) )/real(nsub,8)
     dx2= ( grid%x(:,nd(2)) - grid%x(:,nd(3)) )/real(nsub,8)
     area= 0.5D0*abs( dx1(1)*dx2(2) - dx2(1)*dx1(2) )

     do j=1, nsub       ! for all subtriangle
        do k=1, 2*j-1 ! for all subtriangle
           itmp=  (k-1)/2   ! make integer if not integer
           c1=(j-1)- itmp + real(mod(k,2)*2-1,8)/3D0
           itmp= k/2
           c2=itmp + real(mod(k,2)*2-1,8)/3D0
           c3=real(nsub,8) - c1 - c2

           xff=grid%x(:,nd(3)) + c1*dx1 + c2*dx2 ! center of subtriangle
           !volume=area*xff(1)    ! volume of subtriangle (2pi is missing)


#ifdef XGC1
           ! coordinate transform
           call field_following_pos2(xff,phi0,phi,x)
#else
           x=xff
#endif

           ! find gyro-ring position
           dpdr=psi_interpol(x(1),x(2),1,0)  !### get together 3 of them (psi, dpdr, dpdz)
           dpdz=psi_interpol(x(1),x(2),0,1)
           dp=sqrt(dpdr**2 + dpdz**2)
           cosa=dpdr/dp
           sina=dpdz/dp


           do larmor=1, n_gyro

              !find 8-point with rhoi
              dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
              dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa

              x_ring = x + rho* dx_unit2(:)

#ifdef XGC1
              ! back transform to ff position
              call field_following_pos2(x_ring,phi,phi0,x_rff)
              ! find triangle
              call search_tr2(grid,x_rff,itr,p)
#else
              x_rff=x_ring
              call search_tr2(grid,x_rff,-1D10,itr,p)
#endif


              if(itr>0) then
                 do ip=1,3
                    weight=p(ip)/real(n_gyro)
                    node=grid%nd(ip,itr)
                    call set_value(psn%gyro_avg_mat,nd(1),node,weight*c1*area,1)
                    call set_value(psn%gyro_avg_mat,nd(2),node,weight*c2*area,1)
                    call set_value(psn%gyro_avg_mat,nd(3),node,weight*c3*area,1)
                 enddo
              endif

           enddo
        enddo
     enddo
  enddo

  !normalize average operation
  call normalize_const(psn%gyro_avg_mat)

contains
  subroutine normalize_const(mat)
    implicit none
    type(mat_type) :: mat
    integer :: i,j
    real (8) :: sum

    do i=1, mat%n
       sum=1D-99
       do j=1, mat%nelement(i)
          sum=sum + mat%value(j,i)
       enddo

       mat%value(:,i)=mat%value(:,i)/sum
    enddo
  end subroutine normalize_const
end subroutine init_gyro_avg_mat

#endif

subroutine apply_wall_boundary_condition(grid,psn,den)
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: den(grid%nnode)
  integer :: i

  do i=1, psn%nwall
     den(psn%wall_nodes(i))=psn%sheath_pot(i)
  enddo

end subroutine apply_wall_boundary_condition

subroutine sheath_pot_init(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i
  real (kind=8) :: x(2), Te, psi
  real (kind=8) :: psi_interpol
! apply sheat potential according to temperature


   do i=1, psn%nwall
      x=grid%x(:,psn%wall_nodes(i))
      psi=psi_interpol(x(1),x(2),0,0)
      ! get initial temperature
      Te=eq_ftn(psi,x(1),x(2),eq_tempe)
      ! apply alpha *  T
      psn%sheath_pot(i)=sml_sheath_init_pot_factor*Te
   enddo


end subroutine

subroutine simple00(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn

  real (8) :: area(grid%npsi00), dadp(grid%npsi00), dpdp00(grid%npsi00)
  integer :: npsi, i, bd1, bd2, bd3, bd4
  real (8) :: psi, in_bd_psi, out_bd_psi


  npsi = grid%npsi00

  ! prepare coefficient
  dadp = psn%vol00 / (sml_2pi*eq_axis_r*grid%dpsi00)  ! not real area -- approx sense.

  area(1)=dadp(1)  ! area defined in the middle point

  do i=2, npsi
     area(i)=area(i-1) + dadp(i)
  enddo

  ! prepare boundary
  bd1=1   !potential boundary - where potential is zero
  bd2=1
  bd3=npsi
  bd4=npsi
  do i=1, npsi
     psi=grid%psi00min + real(i-1)*grid%dpsi00

     in_bd_psi  = sml_inpsi  + sml_bd_ext_delta1*eq_x_psi
     out_bd_psi = sml_outpsi - sml_bd_ext_delta3*eq_x_psi
!     if(sml_rgn1_pot0_only) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta3), out_bd_psi)
     if(.true.) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta3), out_bd_psi) !enforce region1 only

     if(sml_zero_inner_bd/=1) in_bd_psi=-99D0  ! some minus number
     if( psi <= in_bd_psi  .and. psi + grid%dpsi00 >   in_bd_psi ) bd1=i
     if( psi < out_bd_psi  .and. psi + grid%dpsi00 >= out_bd_psi ) bd3=i+1

     in_bd_psi  = sml_inpsi  + sml_bd_ext_delta2*eq_x_psi
     out_bd_psi = sml_outpsi - sml_bd_ext_delta4*eq_x_psi
!     if(sml_rgn1_pot0_only) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta4), out_bd_psi)
     if(.true.) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta4), out_bd_psi) !enforce region1 only

     if(sml_zero_inner_bd/=1) in_bd_psi=-99D0  ! some minus number
     if( psi <= in_bd_psi  .and. psi + grid%dpsi00 >   in_bd_psi ) bd2=i
     if( psi < out_bd_psi  .and. psi + grid%dpsi00 >= out_bd_psi ) bd4=i+1
  enddo



  dpdp00=0D0  ! defined in the middle point
  do i=bd2 + 1, bd4  ! bd3?
     dpdp00(i) = dpdp00(i-1) + psn%cden00_1d(i)*dadp(i)/sml_2pi * grid%dpsi00
  enddo

  dpdp00(1:npsi-1)=dpdp00(1:npsi-1)*0.5D0*(dadp(1:npsi-1) + dadp(2:npsi)) / ( 2D0*area(1:npsi-1)*psn%mn_eb2(1:npsi-1) )


  psn%pot00_1d=0D0
  do i=bd1 + 1, bd3
     psn%pot00_1d(i) = psn%pot00_1d(i-1) + dpdp00(i-1)
  enddo

  psn%pot00_1d(1:bd3)=psn%pot00_1d(1:bd3) - psn%pot00_1d(bd3)

end subroutine simple00

subroutine init_simple00(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i
  real (8) :: psi, den

  allocate(psn%mn_eb2(grid%npsi00))

  do i=1, grid%npsi00

     psi=grid%psi00min + (real(i)-0.5D0)*grid%dpsi00
     den=eq_ftn(psi,eq_axis_r,eq_axis_z,eq_den)

     psn%mn_eb2(i)=-ptl_mass(1)/sml_e_charge*den/(eq_axis_b*eq_axis_b)
  enddo
end subroutine init_simple00

