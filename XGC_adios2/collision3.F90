subroutine vpic_col(grid, spall)
  use ptl_module
  use sml_module
  use col_module
  use perf_monitor
  use grid_class
  implicit none
  type(grid_type) :: grid
  type(species_type), target :: spall(0:ptl_nsp_max)
  type(species_type), pointer :: spi, spe
  integer :: ierror, i, j, k, vpic_ierr
  real (kind=8) :: col_3_massi, col_3_masse
  real (8),external :: psi_interpol
  type(col_3_type) :: col_3_spi, col_3_spe

  real (kind=8) :: vpic_tstep_gamma(4, 0:col_3_nset_allocated-1) ! 1: ion, 2:electron 

  interface
    subroutine vpic_p2m_bl2(grid, sp, col_3_sp)
    use grid_class
    use ptl_module
    use sml_module
    use eq_module
    use col_module
    use perf_monitor
    implicit none
    type(grid_type) :: grid
    type(species_type) :: sp
    type(col_3_type) :: col_3_sp
    end subroutine vpic_p2m_bl2
    subroutine vpic_m2p_bl2(sp, col_3_sp)
    use ptl_module
    use sml_module
    use eq_module
    use col_module
    use perf_monitor
    implicit none
    type(species_type) :: sp
    type(col_3_type) :: col_3_sp
    end subroutine vpic_m2p_bl2
  end interface

  spi=>spall(1)
  if(sml_electron_on) spe=>spall(0)
  

  if(sml_mype .eq. 0) print *, 'enter NL-VPIC'
  col_3_massi = spi%mass                         !Simulation mass
  if(sml_electron_on) col_3_masse = spe%mass                         !Simulation mass

  ! I need to check if there is time-consuming double execution script
  call vpic_p2m_bl2(grid, spi, col_3_spi)
  if(sml_electron_on) call vpic_p2m_bl2(grid, spe, col_3_spe)

  call  col_3_lambda_gamma(col_3_spi, col_3_spe, col_3_massi, col_3_masse, vpic_tstep_gamma)

  if(sml_electron_on) then
      do k=0, col_3_nset_allocated-1
          if(sml_mype .eq. 0) print *, 'col_3_core :  in', k
          vpic_ierr=0
          if(col_3_spi%num(k) .lt. col_3_min_popul .or. col_3_spe%num(k) .lt. col_3_min_popul ) then
              if(sml_mype .eq. 0) print *, 'the number of sampling particle is too small - col_3 skipped', sml_mype, sml_plane_mype, k, col_3_spi%num(k), col_3_spe%num(k), col_3_min_popul
              vpic_ierr=-11 
          else
                  call col_3_core_m(spi, col_3_spi%lz(k), col_3_spi%deltay(k), col_3_spi%deltax(k), col_3_spi%col_3_vmap(0:col_3_nvr-1,0:col_3_nvz-1,k), col_3_spi%vol(0:col_3_nvr-1,k), &
                                    spe, col_3_spe%lz(k), col_3_spe%deltay(k), col_3_spe%deltax(k), col_3_spe%col_3_vmap(0:col_3_nvr-1,0:col_3_nvz-1,k), col_3_spe%vol(0:col_3_nvr-1,k), &
                                    vpic_tstep_gamma(:,k), vpic_ierr)
          endif
          if(sml_mype .eq. 0) print *, 'col_3_core : out', k
      enddo
  else
      do k=0, col_3_nset_allocated-1
          if(sml_mype .eq. 0) print *, 'col_3_core :  in', k
          vpic_ierr=0
          if(col_3_spi%num(k) .lt. col_3_min_popul ) then
              if(sml_mype .eq. 0) print *, 'the number of sampling particle is too small - col_3 skipped', sml_mype, sml_plane_mype, k, col_3_spi%num(k)
              vpic_ierr=-11 
          else
                  call col_3_core_s(spi, col_3_spi%lz(k), col_3_spi%deltay(k), col_3_spi%deltax(k), col_3_spi%col_3_vmap(0:col_3_nvr-1,0:col_3_nvz-1,k), col_3_spi%vol(0:col_3_nvr-1,k), &
                                    vpic_tstep_gamma(1,k), vpic_ierr)
          endif
          if(sml_mype .eq. 0) print *, 'col_3_core : out', k
      enddo
  endif

  call vpic_m2p_bl2(spi, col_3_spi)
  if(sml_electron_on) call vpic_m2p_bl2(spe, col_3_spe)

  !multispecies conservation routine is required.
  !call consv_end(sp, col_3_sp)

  call col_3_sp_deallocation(col_3_spi)
  if(sml_electron_on) call col_3_sp_deallocation(col_3_spe)

end subroutine vpic_col

subroutine col_3_core_s(sp, Dlx, mesh_dr, mesh_dz, dist_n, vol, vpic_gamma, vpic_ierr)
    use sml_module, only : sml_e_charge, sml_pi, sml_mype, sml_plane_mype
    use ptl_module 
    use col_module
    implicit none
    type(species_type) :: sp
    real (kind=8) :: Dlx, mesh_dz, mesh_dr  !original code : -FPL_rx->Dlx
    real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_n, dist_iter   ! local 
    real (kind=8), dimension(1:col_3_nvr) :: vol      ! local
    real (kind=8), dimension (col_3_ntotal_v) :: dist_col, dist_n_col
    real (kind=8) :: vpic_gamma
    integer :: vpic_ierr
    integer :: mesh_Nr, mesh_Nz, mesh_N
    integer :: vpic_inner_iter_max, mesh_Nrm1, mesh_Nzm1, mesh_Nrzm1
    real (kind=8) :: vpic_tstep, vpic_tstep_gamma 
    real (kind=8), dimension (:,:,:), allocatable :: EDs
    real (kind=8), dimension (:,:), allocatable :: f_half, dfdr, dfdz
    real (kind=8), dimension (:,:,:), allocatable :: Ms 
    real (kind=8) :: numeric_T, vpic_mass, tmpr1, tmpr2, tmpr3
    integer :: index_ip, index_jp, itmp, index_I, index_J
    integer :: iter_inter
    real (kind=8) :: vpic_dn,vpic_dw,vpic_dw_prev, vpic_dfc, vpic_exit_prev_num, vpic_exit_prev_en, vpic_exit_num, vpic_exit_en
    ! Super LU Variables
    real (kind=8), dimension(LU_nnz) :: LU_values 
    !-------------------------------------------------------------Super LU
    integer :: i, j, m, ierr
    type(col_3_core_type) :: cs

    if(sml_mype .eq. 0) print *, 'col_3 single species'

    mesh_Nr = col_3_nvr
    mesh_Nz = col_3_nvz
    mesh_N = col_3_ntotal_v

    mesh_Nrm1 = mesh_Nr-1
    mesh_Nzm1 = mesh_Nz-1
    mesh_Nrzm1 = mesh_Nrm1*mesh_Nzm1

    vpic_inner_iter_max = 100
    vpic_tstep = col_3_dt
    vpic_tstep_gamma = vpic_tstep*vpic_gamma 

    vpic_exit_prev_num = 1D0
    vpic_exit_prev_en = 1D0

    vpic_ierr = 0
    vpic_mass = sp%mass

    call col_3_core_s_init(sp, cs, dist_n, mesh_dr, mesh_dz, dlx, vol)
    !call check_point('col_3_core_s_init')

    allocate(Ms(4,mesh_Nrzm1,mesh_Nrzm1)) ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra 
    Ms = 0D0
    call col_3_angle_avg_s(cs, Ms)

    dist_iter = dist_n
    dist_n_col = reshape(dist_n, (/mesh_Nr*mesh_Nz/))         !FOR COLUMN for b in every linear algebra

    do iter_inter=1, vpic_inner_iter_max
    !NOTE THAT "dist_iter" is the PDF iterated and being updated in this LOOP

       dist_col = dist_n_col                                 !dist_col SHOULD BE always "f_n"

       !f_(I+1/2,J+1/2) = f_half, dfdr, dfdz
       allocate( f_half(mesh_Nrm1, mesh_Nzm1), dfdr(mesh_Nrm1, mesh_Nzm1), dfdz(mesh_Nrm1, mesh_Nzm1) )
       call col_3_f_df(cs%delta_r(1,:), cs%delta_z(1,:), cs%mesh_dr, cs%mesh_dz, dist_iter, f_half, dfdr, dfdz)
       !call check_point('col_3_f_df')

       ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
       allocate(EDs(6,mesh_Nrm1, mesh_Nzm1))
       call col_3_E_and_D_s(cs, cs, f_half, dfdr, dfdz, Ms, EDs)  
       !call check_point('col_3_E_and_D')
       deallocate( f_half, dfdr, dfdz )

       call col_3_LU_matrix(1, cs, EDs, iter_inter, LU_values)
       !call check_point('col_3_LU_matrix')
       LU_values = LU_values*vpic_tstep_gamma
       
       deallocate(EDs)

       call col_3_picard_step(iter_inter, LU_values, dist_col, dist_iter)

       ! DIST_ITER has new updated PDF.

       vpic_dn = 0D0
       vpic_dw = 0D0
       vpic_dw_prev = 0D0
       do index_I=1,mesh_Nz
           tmpr1 = cs%mesh_z(index_I)
           tmpr1 = tmpr1*tmpr1     ! mesh_z^2
           do index_J=1,mesh_Nr
               tmpr2 = cs%mesh_r(index_J)
               tmpr2 = tmpr2*tmpr2     ! mesh_r^2
               tmpr3 = tmpr1+tmpr2     ! mesh_z^2+mesh_r^2
               ! The error difference between guide program(MATLAB) and this stems from dfdtc and dfc
               vpic_dfc = (dist_iter(index_J,index_I) - dist_n(index_J, index_I))*vol(index_J)
               vpic_dn = vpic_dn + vpic_dfc
               vpic_dw = vpic_dw +  vpic_dfc * tmpr3
               vpic_dw_prev = vpic_dw_prev + dist_n(index_J,index_I)*tmpr3*vol(index_J)
           enddo
       enddo
!#ifdef COL_3_CORE_MSG
       if(sml_mype .eq. 0) then
           print *, 'dn/n = ', vpic_dn/cs%dens, 'dw/w = ', vpic_dw/vpic_dw_prev 
       endif
!#endif       
       ! exit process : we might need to tag number for each exit method
       vpic_exit_en = dabs(vpic_dw/vpic_dw_prev)
       vpic_exit_num = dabs(vpic_dn/cs%dens)
       if( vpic_exit_en .le. 1D-10 .and. vpic_exit_num .le. 1D-10 .and. iter_inter .ne. 1 ) then
           ! Accuracy is enough to go
           exit
       else if( (1D0 - vpic_exit_prev_en/vpic_exit_en) .lt. 1D-1 &
                  .and. (1D0 - vpic_exit_prev_num/vpic_exit_num) .lt. 1D-1 &
                  .and. vpic_exit_en .le. 1D-5 .and. vpic_exit_num .le. 1D-5 &
                  .and. iter_inter .ne. 1 ) then
           ! Slow convergence exit
           exit
       else
           ! FOR ITERATION
           vpic_exit_prev_en = vpic_exit_en
           vpic_exit_prev_num = vpic_exit_num
           if( iter_inter .eq. 50 ) print *, sml_mype,': iteration is over 50', vpic_exit_num, vpic_exit_en
       endif

    enddo !iter_inter

    if (iter_inter .eq. vpic_inner_iter_max) then
        print *, 'inner iteration went to maximum number at :', sml_mype, sml_plane_mype,  vpic_dn/cs%dens, vpic_dw/vpic_dw_prev
    endif

    dist_n = dist_iter  !This line has been changed from original source for simplicity

    deallocate(Ms)

    call col_3_core_s_deallocation(cs) 

end subroutine col_3_core_s

subroutine col_3_core_s_deallocation(cs)
    use col_module
    implicit none
    type(col_3_core_type) :: cs

    deallocate(cs%mesh_r, cs%mesh_z, cs%mesh_r_half, cs%mesh_z_half, cs%local_center_volume, cs%vol)
    deallocate(cs%delta_r, cs%delta_z)
end subroutine col_3_core_s_deallocation

subroutine col_3_lambda_gamma(col_3_spi, col_3_spe, col_3_massi, col_3_masse, gammac)
   use sml_module, only : sml_pi, sml_e_charge, sml_electron_on
   use col_module
   implicit none
   type(col_3_type) :: col_3_spi, col_3_spe
   real (kind=8) :: col_3_massi, col_3_masse
   real (kind=8) :: lambda(3, 0:col_3_nset_allocated-1) ! 1: i-i, 2:i-e, 3: e-e
   real (kind=8) :: gammac(4, 0:col_3_nset_allocated-1) ! 1: i-i, 2:i-e, 3: e-i, 4: e-e 
   integer :: i

  ! NRL Plasma Formulary 2011 version, p.34
  !(c) Mixed ion-ion collisions (here, same species ion-ion collisions) 
  lambda(1,:) = 2.3D1 - log(sqrt( 2D0*col_3_spi%den*1D-6/col_3_spi%t_ev)/col_3_spi%t_ev)
  
  ! i-i
  gammac(1,:) = (sml_e_charge**4) * lambda(1,:) / ( col_3_massi*((8.8542D-12)**2) * 8D0 * sml_pi )
   
  if( sml_electron_on ) then 
      !(b) Electron-ion (and ion-electron) collisions 
      do i=0, col_3_nset_allocated-1
          if(col_3_spi%t_ev(i) * col_3_masse/col_3_massi .lt. col_3_spe%t_ev(i)) then
              if(col_3_spe%t_ev(i) .lt. 1D1) then
                  lambda(2,i) = 2.3D1 - log(sqrt(col_3_spe%den(i)*1D-6/(col_3_spe%t_ev(i)**3)))
              else
                  lambda(2,i) = 2.4D1 - log(sqrt(col_3_spe%den(i)*1D-6)/col_3_spe%t_ev(i))
              endif
          else
              lambda(2,i) = 3D1 - log(sqrt(col_3_spi%den(i)*1D-6/(col_3_spi%t_ev(i)**3)))
          endif
      enddo

      !(a) Thermal electron-electron collisions
      lambda(3,:) = 2.35D1 - log(sqrt(col_3_spe%den*1D-6)*(col_3_spe%t_ev**(-1.25D0))) - sqrt(1D-5+((log(col_3_spe%t_ev) -2D0)**2)/16D0)

      ! i-e
      gammac(2,:) = (sml_e_charge**4) * lambda(2,:) / ( col_3_massi*((8.8542D-12)**2) * 8D0 * sml_pi ) 
      ! e-i
      gammac(3,:) = (sml_e_charge**4) * lambda(2,:) / ( col_3_masse*((8.8542D-12)**2) * 8D0 * sml_pi ) 
      ! e-e  
      gammac(4,:) = (sml_e_charge**4) * lambda(3,:) / ( col_3_masse*((8.8542D-12)**2) * 8D0 * sml_pi ) 
  endif


end subroutine col_3_lambda_gamma

subroutine col_3_sp_deallocation(col_3_sp)
  use col_module
  implicit none
  type(col_3_type) :: col_3_sp

   deallocate(col_3_sp%lz, col_3_sp%lr, col_3_sp%Vel_avg_z, col_3_sp%deltax, col_3_sp%deltay)
   deallocate(col_3_sp%den, col_3_sp%t_ev)
   deallocate(col_3_sp%Evolume_inV_proc, col_3_sp%vol)
   deallocate(col_3_sp%col_3_vmap, col_3_sp%col_3_vmap_prev)
   deallocate(col_3_sp%ipsth, col_3_sp%ivel_r, col_3_sp%ivel_z)
   deallocate(col_3_sp%vpic_ptl_v_parallel, col_3_sp%vpic_ptl_v_perp)
   deallocate(col_3_sp%p2m_fc)
   deallocate(col_3_sp%org_weight_sum, col_3_sp%org_v2_avg)
   deallocate(col_3_sp%negative_weight_sum, col_3_sp%negative_v2_sum, col_3_sp%negative_vp_sum)
   deallocate(col_3_sp%negative_tag)
   deallocate(col_3_sp%num)

end subroutine col_3_sp_deallocation

subroutine vpic_p2m_bl2(grid, sp, col_3_sp)
  use grid_class
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  type(species_type) :: sp
  type(col_3_type) :: col_3_sp

  real (kind=8) :: r, z, b, psi, rho, r_cyl, theta, mu, weight
  real (kind=8) :: c_m_sp, vpic_mass
  real (kind=8) :: dum_w(0:col_3_nset_allocated-1)
  real (kind=8) :: v_para_min(0:col_3_nset_allocated-1)
  real (kind=8) :: v_para_max(0:col_3_nset_allocated-1)
  real (kind=8) :: dz, dr, dr2, dzpi2dr2
  real (kind=8) :: domain_l, ddum, ddum0, ddum1, C_par, C_perp
  integer :: Nmesh_R_local,Nmesh_R_local_m1,Nmesh_R_local_m2, posX, posY
  integer :: i, j, k, indx, ierror, itarget, rtarget, ztarget
  integer :: idumA(0:col_3_nset_allocated-1)
  real (8),external :: psi_interpol, b_interpol
  type(ptl_type) :: ptli
  integer, external :: col_3_search_indx_subg

   !allocation
   !lz : max(abs(v_para)), lr : max(v_perp), Vel_avg_z : avg v_para
   allocate(col_3_sp%lz(0:col_3_nset_allocated-1), col_3_sp%lr(0:col_3_nset_allocated-1), &
            col_3_sp%Vel_avg_z(0:col_3_nset_allocated-1),&
            col_3_sp%deltax(0:col_3_nset_allocated-1), col_3_sp%deltay(0:col_3_nset_allocated-1),&
            col_3_sp%den(0:col_3_nset_allocated-1), col_3_sp%t_ev(0:col_3_nset_allocated-1))
   !Evolume_inV_proc : v space volume?
   allocate(col_3_sp%Evolume_inV_proc(0:col_3_nvr-1,0:col_3_nset_allocated-1), &
            col_3_sp%vol(0:col_3_nvr-1,0:col_3_nset_allocated-1))
   !vmap : vspace arr
   allocate(col_3_sp%col_3_vmap(0:col_3_nvr-1,0:col_3_nvz-1,0:col_3_nset_allocated-1),&
            col_3_sp%col_3_vmap_prev(0:col_3_nvr-1,0:col_3_nvz-1,0:col_3_nset_allocated-1))

   !ipsth : cell number, ivel_r : r location in vspace, ivel_z : z location in vspace
   allocate(col_3_sp%ipsth(sp%num), col_3_sp%ivel_r(sp%num), col_3_sp%ivel_z(sp%num))
   !vpic_ptl_v_parallel : v_para in moving frame, vpic_ptl_v_perp : v_perp
   allocate(col_3_sp%vpic_ptl_v_parallel(sp%num), col_3_sp%vpic_ptl_v_perp(sp%num))
   allocate(col_3_sp%p2m_fc(4,sp%num))
   allocate(col_3_sp%org_weight_sum(0:col_3_nset_allocated-1), &
            col_3_sp%org_v2_avg(0:col_3_nset_allocated-1))
   allocate(col_3_sp%num(0:col_3_nset_allocated-1))

  col_3_sp%col_3_vmap=0D0
  c_m_sp=ptl_c_m(sp%type)
  vpic_mass = sp%mass
  idumA=0
  col_3_sp%ipsth=-1
  col_3_sp%Vel_avg_z=0D0
  dum_w=1D-99
  v_para_min=1D50
  v_para_max=-1D50
  col_3_sp%lr=0D0
  col_3_sp%den=0D0
  col_3_sp%t_ev=0D0
  col_3_sp%org_weight_sum=0D0
  col_3_sp%org_v2_avg=0D0

  !1. evaluate lz, lr, Vel_avg_z, ipsth, vpic_ptl_v_parallel, vpic_ptl_v_perp 

  do i=1, sp%num
     if(sp%ptl(i)%gid > 0) then
        ptli = sp%ptl(i) ! the role of ptli would be better to be replaced as a pointer

        r=ptli%ph(1)
        z=ptli%ph(2)
        b=b_interpol(r,z,0D0)
        psi=psi_interpol(r,z,0,0)
        if(psi .gt. sml_outpsi .or. psi .lt. sml_inpsi) cycle

        rho=ptli%ph(pirho)
        mu =ptli%ct(pim)

        if (ptl_deltaf_sp(sp%type)) then
           !below line : f = f0(0) + delta f + df0(t)
           weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
        else
           weight = ptli%ct(piw0)
        endif

        r_cyl = sqrt((r-eq_axis_r)*(r-eq_axis_r)+(z-eq_axis_z)*(z-eq_axis_z))
        theta= acos((r-eq_axis_r)/r_cyl)
        if(z .lt. eq_axis_z) then
            theta= sml_2pi-theta
        endif

        indx = col_3_search_indx_subg(r,z,psi,theta,grid,2)
        if(indx .eq. -1) cycle  !This is limited case 
!ipsth setup
        do j=0, col_3_nset_allocated-1
           if(indx .eq. col_3_tag(j)) then
              col_3_sp%ipsth(i)=j
              exit
           endif
        enddo

        if (col_3_sp%ipsth(i) < 0) then
           print *,'ipsth find error', col_3_sp%ipsth(i), indx, col_3_nset_allocated-1, col_3_tag
        endif

        idumA(col_3_sp%ipsth(i))=idumA(col_3_sp%ipsth(i))+1

!v_para and v_perp setup
        col_3_sp%vpic_ptl_v_parallel(i)=c_m_sp*rho*b !in lab frame
        col_3_sp%vpic_ptl_v_perp(i)=sqrt(mu*b*2D0/vpic_mass)

        col_3_sp%Vel_avg_z(col_3_sp%ipsth(i))=col_3_sp%Vel_avg_z(col_3_sp%ipsth(i))+&
                                     col_3_sp%vpic_ptl_v_parallel(i)*weight
        dum_w(col_3_sp%ipsth(i))=dum_w(col_3_sp%ipsth(i))+weight

        v_para_min(col_3_sp%ipsth(i))=min(v_para_min(col_3_sp%ipsth(i)), col_3_sp%vpic_ptl_v_parallel(i))
        v_para_max(col_3_sp%ipsth(i))=max(v_para_max(col_3_sp%ipsth(i)), col_3_sp%vpic_ptl_v_parallel(i))

        col_3_sp%lr(col_3_sp%ipsth(i))=max(col_3_sp%lr(col_3_sp%ipsth(i)), col_3_sp%vpic_ptl_v_perp(i))
     endif
  enddo
  col_3_sp%num = idumA

  do i=0, col_3_nset_allocated-1
      if(col_3_sp%num(i) .lt. col_3_min_popul) then
          print *, 'The number of particles in a certain cell is insufficient for col_3. no collision for this cell', sml_mype, col_3_tag(i), col_3_sp%num(i), col_3_min_popul
          !ipsth for corresponding particles to be -1
      endif
  enddo
  
  !if(minval(idumA)==0) then
  !   print *, 'the number of particles for VPIC collision is not enough.'
  !   print *, 'Empty bin for collision exists.'
  !   stop
  !endif

!Vel_avg_z setup
  col_3_sp%Vel_avg_z=col_3_sp%Vel_avg_z/dum_w

!lr and lz setup
  col_3_sp%lr=col_3_sp%lr*1.001D0
  do j=0, col_3_nset_allocated-1
     col_3_sp%lz(j)=max(abs(v_para_min(j)-col_3_sp%Vel_avg_z(j)), abs(v_para_max(j)-col_3_sp%Vel_avg_z(j)))
  enddo
  col_3_sp%lz=col_3_sp%lz*1.001D0

!deltax : vpara dv size, deltay : vperp dv size
  col_3_sp%deltax=2D0*col_3_sp%lz/real(col_3_nvz-1)
  col_3_sp%deltay=col_3_sp%lr/real(col_3_nvr-1)

! lz setup
 col_3_sp%lz = col_3_sp%Vel_avg_z - col_3_sp%lz

!volume (or vol) : volume in velocity space
  col_3_sp%vol = 0D0
  do i=0,col_3_nset_allocated-1
      Nmesh_R_local = col_3_nvr 
      ddum0 = col_3_sp%deltax(i)
      ddum1 = col_3_sp%deltay(i)

      ddum = sml_2pi* ddum1 * ddum1 * ddum0
      !volume(0,k) = 0.25D0*sml_pi*(tmp_r2*tmp_r2)*tmp_r1
      col_3_sp%vol(0,i) = 0.125D0*ddum
      do j=1,Nmesh_R_local-1
          col_3_sp%vol(j,i) = ddum*j
      enddo
  enddo 

!Evolume_inV_proc setup
  do i=0, col_3_nset_allocated-1
     dz = col_3_sp%deltax(i)
     dr = col_3_sp%deltay(i)
     dr2 = dr*dr
     Nmesh_R_local=col_3_nvr
     Nmesh_R_local_m2 = Nmesh_R_local-2

     dzpi2dr2 = dz*sml_pi*2D0*dr2

     do j=1, Nmesh_R_local_m2
        col_3_sp%Evolume_inV_proc(j,i) = dzpi2dr2*j
     enddo ! j
     col_3_sp%Evolume_inV_proc(0,i) = dzpi2dr2 / 6D0
     col_3_sp%Evolume_inV_proc(Nmesh_R_local-1,i) = dzpi2dr2 *(0.5D0*real(Nmesh_R_local)-2D0/3D0)
  enddo ! i


  do i=1, sp%num
     if(sp%ptl(i)%gid > 0 .and. col_3_sp%ipsth(i) .ne. -1) then
        ptli = sp%ptl(i)
        itarget = col_3_sp%ipsth(i)
        if(col_3_sp%num(itarget) .lt. col_3_min_popul) then
            col_3_sp%ipsth(i) = -1 
            cycle
        endif
     
        dz = col_3_sp%deltax(itarget)
        dr = col_3_sp%deltay(itarget)
        domain_l = col_3_sp%lz(itarget)
        Nmesh_R_local=col_3_nvr
        Nmesh_R_local_m1 = Nmesh_R_local-1

        if (ptl_deltaf_sp(sp%type)) then
           weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
        else
           weight = ptli%ct(piw0)
        endif

!v_para in moving frame
        !col_3_sp%vpic_ptl_v_parallel(i)=col_3_sp%vpic_ptl_v_parallel(i)-col_3_sp%Vel_avg_z(itarget)

        ddum0=col_3_sp%vpic_ptl_v_parallel(i)
        ddum1=col_3_sp%vpic_ptl_v_perp(i)

        col_3_sp%den(itarget)=col_3_sp%den(itarget)+weight
        col_3_sp%t_ev(itarget)=col_3_sp%t_ev(itarget)+&
                               weight*0.5D0*vpic_mass*(ddum0**2+ddum1**2)
        col_3_sp%org_weight_sum(itarget)=col_3_sp%org_weight_sum(itarget)+weight
        col_3_sp%org_v2_avg(itarget)=col_3_sp%org_v2_avg(itarget)+&
                               weight*(ddum0**2+ddum1**2)
        posX = floor((ddum0-domain_l)/dz)
        !DEBUGGING 26JUL11
        if( posX*dz+domain_l .gt. ddum0 ) then
           print *, 'vpic : rounding error correction'
           posX = posX-1
        endif
        posY = floor(ddum1/dr)
        if( posY*dr .gt. ddum1 ) then
           print *, 'vpic : rounding error correction'
           posY = posY-1
        endif
        if(posX .gt. ubound(col_3_sp%col_3_vmap,2) .or. posX .lt. lbound(col_3_sp%col_3_vmap,2) &
           .or. posY .gt. ubound(col_3_sp%col_3_vmap,1) .or. posY .lt. lbound(col_3_sp%col_3_vmap,1)) then
            print *, 'array bound error in p2m_bl2 of col_3:', sml_mype, i, col_3_sp%ipsth(i), ptli%gid, col_3_sp%num(col_3_sp%ipsth(i)),ddum0, domain_l, dz, ddum1, dr, col_3_sp%Vel_avg_z(itarget)
        endif

!ivel_r, ivel_z setup
        col_3_sp%ivel_r(i)=posY
        col_3_sp%ivel_z(i)=posX

        !------------------------- bilinear interpolation
        ddum0 = (ddum0-domain_l)/dz-posX   !!local normalized v_z in cell
        ddum1 = ddum1/dr- posY                  !!local normalized v_r in cell
              
        C_par = 1D0 - ddum0
        C_perp = 1D0 - ddum1
              
!p2m_fc setup
        col_3_sp%p2m_fc(1,i) = C_par          
        col_3_sp%p2m_fc(2,i) = C_perp
        
        ! Note that vmap is not actual distribution function in that real volume has not been normalized(devided) 
        ! It is hard to evaluate real volume because collision sections are splited based on "psi", not "r"
        ! Normalized distribution function need to be used for FP operation.

        if( weight .lt. 0D0 .or. C_par .lt. 0D0 .or. C_perp .lt. 0D0 .or. C_par .gt. 1D0 .or. C_perp .gt. 1D0 ) then
            print *, 'in : weight error 1:', posX, posY, itarget, Nmesh_R_local
            print *, 'in : weight error 2:', ddum0, dz, domain_l
            print *, weight, C_par, C_perp
            print *, 'error position info', sml_mype, i
            if(C_par .gt. 1D0) then
                C_par = 1D0
                col_3_sp%p2m_fc(1,i) = C_par
                print *, 'weight error due to C_par has been corrected to', C_par
                print *, 'this is due to function acurracy of floor'
            else
                stop
            endif
        endif
        ! SW
        col_3_sp%col_3_vmap(posY, posX, itarget) = col_3_sp%col_3_vmap(posY, posX, itarget) +&
                                                   weight*C_par*C_perp
        ! NW
        col_3_sp%col_3_vmap(posY+1, posX, itarget) = col_3_sp%col_3_vmap(posY+1, posX, itarget) +&
                                                     weight*C_par*(1D0-C_perp)
        ! SE
        col_3_sp%col_3_vmap(posY, posX+1, itarget) = col_3_sp%col_3_vmap(posY, posX+1, itarget) +&
                                                     weight*(1D0-C_par)*C_perp
        ! NE 
        col_3_sp%col_3_vmap(posY+1, posX+1,itarget) = col_3_sp%col_3_vmap(posY+1, posX+1,itarget) +&
                                                      weight*(1D0-C_par)*(1D0-C_perp)
     endif
  enddo
  col_3_sp%org_v2_avg = col_3_sp%org_v2_avg /col_3_sp%org_weight_sum !this line was missed - es 

  !density, temperature calculation.
  col_3_sp%t_ev=col_3_sp%t_ev/col_3_sp%den/sml_e_charge*2D0/3D0
  do k=0, col_3_nset_allocated-1
     i=col_3_tag(k)
     col_3_sp%den(k)=col_3_sp%den(k) / col_3_vol_r(i)
  enddo

  ! 8.1.2 Contribution of each particles
  !     Using saved fraction of weight, 
  !       save [(reduced p.d.f.) / ((saved fraction of weight) * (ptcl_weight))]
  do i=1, sp%num
     if(sp%ptl(i)%gid > 0 .and. col_3_sp%ipsth(i) .ne. -1) then
        ptli = sp%ptl(i)

        ddum0 = col_3_sp%p2m_fc(1,i) ! parallel
        ddum1 = col_3_sp%p2m_fc(2,i) ! perp
        if (ptl_deltaf_sp(sp%type)) then
           !below line : f = f0(0) + delta f + df0(t)
           weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
        else
           weight = sp%ptl(i)%ct(piw0)
        endif
        itarget = col_3_sp%ipsth(i)
        rtarget=col_3_sp%ivel_r(i)
        ztarget=col_3_sp%ivel_z(i)

        col_3_sp%p2m_fc(1,i) = weight * ddum0 * ddum1 /&
                               col_3_sp%col_3_vmap(rtarget,ztarget,itarget)
        col_3_sp%p2m_fc(2,i) = weight * (1D0-ddum0) * ddum1 /&
                               col_3_sp%col_3_vmap(rtarget,ztarget+1,itarget)
        col_3_sp%p2m_fc(3,i) = weight * ddum0 * (1D0-ddum1) /&
                               col_3_sp%col_3_vmap(rtarget+1,ztarget,itarget)
        col_3_sp%p2m_fc(4,i) = weight * (1D0-ddum0) * (1D0-ddum1) /&
                               col_3_sp%col_3_vmap(rtarget+1,ztarget+1,itarget)
     endif
  enddo

  !=======================================================================
  !  Construction of Distribution Function
  !  1 : N -> n   operation : dividing N by volume of configuration space
  !  2 : n -> f   operation : dividing n by volume of velocity space
  !======================================================================
    ! n -> n / V_vel
  do k=0, col_3_nset_allocated-1
     i=col_3_tag(k)

     col_3_sp%col_3_vmap(:,:,k)=col_3_sp%col_3_vmap(:,:,k) / col_3_vol_r(i)
  enddo
  do j=0, col_3_nvz - 1
     col_3_sp%col_3_vmap(:,j,:) = col_3_sp%col_3_vmap(:,j,:) /&
                                  col_3_sp%Evolume_inV_proc(:,:)
  enddo

  col_3_sp%col_3_vmap_prev = col_3_sp%col_3_vmap

#ifdef TEST1
 if (sml_mype==0) then
    open(unit=151,file='vmap_all.txt',status='replace')
    open(unit=152,file='vmap_rz.txt',status='replace')
    open(unit=153,file='vmap_xy.txt',status='replace')
    open(unit=154,file='vel_avg.txt',status='replace')

    do i=0, col_3_nvz-1
       do j=0, col_3_nvr-1
          write(151,2035) col_3_sp%col_3_vmap(j,i,:)
       enddo
    enddo

    do i=0,col_3_nset_allocated-1
       write(152,*) col_3_nvr, col_3_nvz
       write(153,*) col_3_sp%deltax(i), col_3_sp%deltay(i)
       write(154,*) col_3_sp%Vel_avg_z(i)
    enddo

    2035 format(800(e19.13,1x))

    close(151)
    close(152)
    close(153)
    close(154)
    print *, 'nl col profile is printed with TEST1 option'
    stop
 endif


#endif
end subroutine vpic_p2m_bl2

subroutine vpic_m2p_bl2(sp, col_3_sp)
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  use perf_monitor
  implicit none
  type(species_type) :: sp
  type(col_3_type) :: col_3_sp

  integer :: i, j, k
  integer :: itarget, rtarget, ztarget
  real (8) :: r, z, psi, weight
  real (8) :: ddum0, ddum1, ddum2, ddum3
  real (8),external :: psi_interpol
  type(ptl_type) :: ptli
  character (len=100) :: fname
  integer :: debug_i, debug_j, debug_k

  do j=0, col_3_nvz-1
      col_3_sp%col_3_vmap(:,j,:)=col_3_sp%col_3_vmap(:,j,:)*col_3_sp%Evolume_inV_proc
  enddo    
  do k=0, col_3_nset_allocated-1
     i=col_3_tag(k)
     col_3_sp%col_3_vmap(:,:,k)=col_3_sp%col_3_vmap(:,:,k) * col_3_vol_r(i)
  enddo

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! 11. Update particle weight (COMPLETE)
  ! ES 21-Feb-11 :
  ! m2p : Back-Operation of p2m
  !       i.e. Saving fractional contribution of particles to each mesh
  !            Then, apply it back from mesh to particles with changed p.d.f
  !       This requires ALLREDUCE for velocity map to gather p.d.f. 
  !             rather than REDUCE_SCATTER.
  !       REDUCE_SCATTER can be used only with domain decomposition method.
  !       In addition, this operation requires 4xPN memory to save 
  !         fractional geometric contribution of each particles to p.d.f,
  !       OR, requires computation for fractional geometric contribution
  !            of all particles again.
  !       For simplicity, the way to save 4xPN double data is implemented
  !          But this can be optional with other 2 options.
  !          1. PN memory(for perp-because of more computation than parallel)
  !                with parallel calculation
  !          2. Recalculation for all direction
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  allocate(col_3_sp%negative_weight_sum(0:col_3_nset_allocated-1), &
           col_3_sp%negative_v2_sum(0:col_3_nset_allocated-1), &
           col_3_sp%negative_vp_sum(0:col_3_nset_allocated-1))
  allocate(col_3_sp%negative_tag(1000))

  col_3_sp%negative_weight_sum = 1D-99
  col_3_sp%negative_v2_sum = 0
  col_3_sp%negative_vp_sum = 0
  col_3_sp%negative_count = 0
  col_3_sp%negative_tag = 0

  do i=1, sp%num
      if(sp%ptl(i)%gid > 0 .and. col_3_sp%ipsth(i) .ne. -1) then
         ptli = sp%ptl(i)
         itarget = col_3_sp%ipsth(i)  !psi&theta index
         rtarget = col_3_sp%ivel_r(i)   !velocity index
         ztarget = col_3_sp%ivel_z(i)   !velocity index

         if( col_3_sp%col_3_vmap(rtarget,ztarget,itarget) .lt. 0D0 &
             .or. col_3_sp%col_3_vmap(rtarget,ztarget+1,itarget) .lt. 0D0 &
             .or. col_3_sp%col_3_vmap(rtarget+1,ztarget,itarget) .lt. 0D0 &
             .or. col_3_sp%col_3_vmap(rtarget+1,ztarget+1,itarget) .lt. 0D0 ) then
             col_3_sp%negative_count = col_3_sp%negative_count+1
             if(col_3_sp%negative_count .gt. 1000) then
                 print *, sml_mype, sml_plane_mype, ': negative count is over 1000, too much negative!'
                 do debug_i=0,col_3_nset_allocated-1 
                    write(fname,'("col_3_negative_",i5.5,"_",i5.5,".txt")') sml_mype, col_3_tag(debug_i)
                    open(unit=777,file=fname,status='replace')
                    do debug_j=0,col_3_nvz-1
                       do debug_k=0,col_3_nvr-1
                          write(777,*) col_3_sp%col_3_vmap(debug_k,debug_j,debug_i) 
                       enddo
                    enddo  
                    close(777)
                 enddo
                 stop
             endif
             col_3_sp%negative_tag(col_3_sp%negative_count) = i

             if (ptl_deltaf_sp(sp%type)) then
                !below line : f = f0(0) + delta f + df0(t)
                weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
             else
                weight = ptli%ct(piw0)
             endif

             ! Shifted Parallel Velocity(ddum0) and Perpendicular Velocity(ddum1)
             ddum0=col_3_sp%vpic_ptl_v_parallel(i) 
             ddum1=col_3_sp%vpic_ptl_v_perp(i)

             col_3_sp%negative_weight_sum(itarget)= col_3_sp%negative_weight_sum(itarget)+weight
             col_3_sp%negative_v2_sum(itarget)    = col_3_sp%negative_v2_sum(itarget)+&
                                                    weight*(ddum0**2+ddum1**2) ! 2 kin_en sum
             col_3_sp%negative_vp_sum(itarget)    = col_3_sp%negative_vp_sum(itarget)+&
                                                    weight*ddum0
         else
             ddum0 = col_3_sp%col_3_vmap(rtarget,ztarget,itarget)    * col_3_sp%p2m_fc(1,i)
             ddum1 = col_3_sp%col_3_vmap(rtarget,ztarget+1,itarget)  * col_3_sp%p2m_fc(2,i)
             ddum2 = col_3_sp%col_3_vmap(rtarget+1,ztarget,itarget)  * col_3_sp%p2m_fc(3,i)
             ddum3 = col_3_sp%col_3_vmap(rtarget+1,ztarget+1,itarget)* col_3_sp%p2m_fc(4,i)

             weight = ddum0 + ddum1 + ddum2 + ddum3

             if (ptl_deltaf_sp(sp%type)) then
                !!below line : f = f0(0) + delta f + df0(t)
                !p2m : weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
                sp%ptl(i)%ph(piw1) = weight/ptli%ct(piw0) - 1D0 + ptli%ph(piw2)  !Note that ptli is not pointer
             else
                sp%ptl(i)%ct(piw0) = weight
             endif
         endif
      endif
  enddo
   

end subroutine vpic_m2p_bl2

subroutine consv_end(sp, col_3_sp)
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(species_type) :: sp
  type(col_3_type) :: col_3_sp

  real (kind=8) :: r, z, phi, b, psi, rho, mu, theta, weight
  real (kind=8) :: ddum0, ddum1

  real (kind=8), allocatable, dimension(:) :: colV_aft_weight_sum, colV_aft_v2_avg, colV_aft_vp_avg
  real (kind=8), allocatable, dimension(:) :: colV_org_weight_sum_eff, colV_org_v2_avg_eff, colV_org_vp_avg_eff
  real (kind=8), allocatable, dimension(:) :: dumr
  real (kind=8), allocatable, dimension(:) :: factor_w, factor_v2, factor_vp
  integer :: i, itarget, ierr, negative_ct
  real (kind=8) :: c_m_sp, vpic_mass
  real (kind=8),external :: b_interpol, psi_interpol
  type(ptl_type) :: ptli

  ALLOCATE(colV_aft_weight_sum(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(colV_aft_v2_avg(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(colV_aft_vp_avg(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(colV_org_weight_sum_eff(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(colV_org_v2_avg_eff(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(colV_org_vp_avg_eff(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(dumr(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(factor_w(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(factor_vp(0:col_3_nset_allocated-1), STAT=ierr)
  ALLOCATE(factor_v2(0:col_3_nset_allocated-1), STAT=ierr)

  colV_aft_weight_sum = 0D0
  colV_aft_v2_avg = 0D0
  colV_aft_vp_avg = 0D0

  c_m_sp=ptl_c_m(sp%type)
  vpic_mass = sp%mass


  negative_ct = 1
  do i=1, sp%num
     if(sp%ptl(i)%gid > 0 .and. col_3_sp%ipsth(i) .ne. -1) then
        ptli = sp%ptl(i)

        r=ptli%ph(1)
        z=ptli%ph(2)
        psi=psi_interpol(r,z,0,0)

        if(sml_inpsi<psi .AND. psi<sml_outpsi) then
           if(col_3_sp%negative_tag(negative_ct) .ne. i) then

              !-------------------------r, z, psi
              itarget = col_3_sp%ipsth(i)
      
              if (ptl_deltaf_sp(sp%type)) then
                 !below line : f = f0(0) + delta f + df0(t)
                 weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
              else
                 weight = ptli%ct(piw0)
              endif

              ! Shifted Parallel Velocity(ddum0) and Perpendicular Velocity(ddum1)
              ddum0=col_3_sp%vpic_ptl_v_parallel(i) 
              ddum1=col_3_sp%vpic_ptl_v_perp(i)
                    
              ! Save moments in shifted velocity domain for conservation
              colV_aft_weight_sum(itarget)=colV_aft_weight_sum(itarget)+weight
              colV_aft_v2_avg(itarget)    =colV_aft_v2_avg(itarget)+ weight*(ddum0**2+ddum1**2) ! 2 kin_en sum
              colV_aft_vp_avg(itarget)    =colV_aft_vp_avg(itarget)+ weight*ddum0 ! v_|| sum
           else
              negative_ct = negative_ct+1
           endif
        endif
     endif
  enddo

  colV_aft_v2_avg=colV_aft_v2_avg/colV_aft_weight_sum
  colV_aft_vp_avg=colV_aft_vp_avg/colV_aft_weight_sum

  ! -------------------------------------OMITTING PARTICLES around NEGATIVE CELL
  ! colV_aft_weight_sum ? containing ommited particles???
  ! colV_org_weight_sum = colV_org_weight_sum - negative_weight_sum
  ! colV_org_vp_avg = colV_org_vp_avg*(1D0+negative_weight_sum/colV_org_weight_sum) - negative_vp_avg*negative_weight_sum/colV_org_weight_sum
  ! colV_org_v2_avg = colV_org_v2_avg*(1D0+negative_weight_sum/colV_org_weight_sum) - negative_v2_avg*negative_weight_sum/colV_org_weight_sum
  ! 28JUL11

  !colV_org_vp_avg_eff = (colV_org_vp_avg*colV_org_weight_sum - negative_vp_sum)/(colV_org_weight_sum-negative_weight_sum)  !for optimization 2012-02-21
  colV_org_vp_avg_eff = - col_3_sp%negative_vp_sum/(col_3_sp%org_weight_sum-col_3_sp%negative_weight_sum)
  colV_org_v2_avg_eff = (col_3_sp%org_v2_avg*col_3_sp%org_weight_sum - col_3_sp%negative_v2_sum)/&
                        (col_3_sp%org_weight_sum-col_3_sp%negative_weight_sum)
  colV_org_weight_sum_eff = col_3_sp%org_weight_sum - col_3_sp%negative_weight_sum

  ! -------------------------------------CORRECTION METHOD
  ! v2_"org"_avg = SUM( alpha*w*gamma*( (vp-beta)^2 + vperp^2 ) / SUM(alpha*w)
  ! w' = alpha * w
  ! vp' = sqrt(gamma)*(vp - beta)
  ! vperp' = sqrt(gamma)*vperp
  ! factor_w = alpha
  ! factor_v2 = sqrt(gamma)
  ! --------------------------------------
  negative_ct = 1
  factor_w  = colV_org_weight_sum_eff / colV_aft_weight_sum  ! Note that factor_w should be MULTIPLIED
  ! Considering this system is in moving frame, "colV_aft_vp_avg"(=factor_vp) should be SUBTRACTED
  ! Note that because of velocity correction, the amount of correction to be corrected is changed!
  factor_v2 = sqrt((colV_org_v2_avg_eff - colV_org_vp_avg_eff**2)/(colV_aft_v2_avg - colV_aft_vp_avg**2))
  factor_vp = colV_aft_vp_avg - colV_org_vp_avg_eff/factor_v2
#ifdef COL_3_DEBUG  
  call mpi_barrier(sml_comm, ierr)
  if(sml_mype == 0) print *, 'Factor calculation has been finished'
  call mpi_barrier(sml_comm, ierr)
#endif
  ! CORRECTIONS need to be proceeded in following order : w -> vp -> v2
  do i=1, sp%num
     if(sp%ptl(i)%gid > 0 .and. col_3_sp%ipsth(i) .ne. -1) then
        ptli = sp%ptl(i)

        r=ptli%ph(1)
        z=ptli%ph(2)
        psi=psi_interpol(r,z,0,0)

        if(sml_inpsi<psi .AND. psi<sml_outpsi) then
           if( (col_3_sp%negative_count .eq. 0) .or. (col_3_sp%negative_tag(negative_ct) .ne. i)) then
              b=b_interpol(r,z,0D0)
              itarget = col_3_sp%ipsth(i)
      
              if (ptl_deltaf_sp(sp%type)) then
                 !below line : f = f0(0) + delta f + df0(t)
                 weight = ptli%ct(piw0)*(1D0+ptli%ph(piw1)-ptli%ph(piw2))
              else
                 weight = ptli%ct(piw0)
              endif

              ! Shifted Parallel Velocity(ddum0) and Perpendicular Velocity(ddum1)
              ddum0=col_3_sp%vpic_ptl_v_parallel(i) 
              ddum1=col_3_sp%vpic_ptl_v_perp(i)
        
              weight = weight * factor_w(itarget)
              ddum0 = (ddum0 - factor_vp(itarget))*factor_v2(itarget)
              ddum1 = ddum1*factor_v2(itarget)
      
              ! Corrections back to particles
              if (ptl_deltaf_sp(sp%type)) then
                 sp%ptl(i)%ph(piw1) = weight/ptli%ct(piw0) - 1D0 + ptli%ph(piw2)
                 sp%ptl(i)%ph(pirho) = (ddum0+col_3_sp%Vel_avg_z(itarget))/(b*c_m_sp)
                 sp%ptl(i)%ct(pim) = (ddum1**2)*0.5D0/b*vpic_mass
              else
                 sp%ptl(i)%ct(piw0) = weight
                 sp%ptl(i)%ph(pirho) = (ddum0+col_3_sp%Vel_avg_z(itarget))/(b*c_m_sp)
                 sp%ptl(i)%ct(pim) = (ddum1**2)*0.5D0/b*vpic_mass
              endif
      
#ifdef VPIC_DEBUG_CONSV
         ! Correction Check
          colV_corr_weight_sum(itarget)=colV_corr_weight_sum(itarget)+weight
          colV_corr_v2_avg(itarget)    =colV_corr_v2_avg(itarget)+ weight*(ddum0**2+ddum1**2) ! 2 kin_en sum
          colV_corr_vp_avg(itarget)    =colV_corr_vp_avg(itarget)+ weight*ddum0 ! v_|| sum
#endif
           else 
              negative_ct = negative_ct+1
           endif
        endif
     endif
  enddo 

  if(negative_ct-1 .ne. col_3_sp%negative_count) then
     print *, 'particle has not been omitted', negative_ct-1, '.ne.', col_3_sp%negative_count
     print *, 'negative_count / particle total number =', col_3_sp%negative_count, sp%num
  endif
#ifdef COL_3_DEBUG  
  call mpi_barrier(sml_comm, ierr)
  if(sml_mype == 0) print *, 'Particle correction has been finished'
  call mpi_barrier(sml_comm, ierr)
#endif
  deallocate(factor_w, factor_vp, factor_v2, dumr)
  deallocate(colV_aft_weight_sum, colV_aft_v2_avg, colV_aft_vp_avg)
  deallocate(colV_org_weight_sum_eff, colV_org_v2_avg_eff, colV_org_vp_avg_eff)

end subroutine consv_end

subroutine col_3_decomp(grid,sp)
  use grid_class
  use ptl_module
  use sml_module
  use eq_module
  use col_module
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(species_type) :: sp
  integer :: i, j, s, m, indx, ierror, n_send_total,n_recv_total, dindx, sep, htotal, findx, fhindx
  integer, allocatable :: dest_count(:), dest_count_cumul(:), dest_count_r(:), &
                          scounts_ptl(:), rcounts_ptl(:), sdispls_ptl(:), rdispls_ptl(:)
  integer, allocatable :: hole_list(:), rgn_indx(:)
  real (8) :: r,z,r_cyl,psi,theta
  real (8),external :: psi_interpol
  logical :: snd_logic
  real (8), allocatable :: grid_psi_norm(:), grid_theta(:)
  integer, allocatable  :: grid_ir(:), grid_itheta(:), grid_ptl_num(:), myyard(:)
  character(len=7) :: mype_num
  integer, external :: col_3_search_indx_subg, col_3_search_indx_subg_decomp 
  
  integer (8):: idum0, idum1, idum2, idum3 
  integer (8) :: ptl_popul(sml_totalpe)

  type col_3_decomp_ptl
  sequence
  real (8) :: ph(ptl_nphase) !6
  real (8) :: ct(ptl_nconst) !3
  real (8) :: phase0(ptl_nphase) !6
  integer :: gid
  end type col_3_decomp_ptl
  
  type(col_3_decomp_ptl), dimension(:), allocatable :: ptl_send, ptl_recv

  ! we assume that toroidal decomposition has been already done.

#ifdef COL_3D_DEBUG
     idum3=sp%num
     call mpi_allreduce(idum3,idum1,1,mpi_integer8,mpi_sum,sml_comm,ierror)
     call mpi_allreduce(idum3,idum2,1,mpi_integer8,mpi_max,sml_comm,ierror)
     call mpi_allreduce(idum3,idum0,1,mpi_integer8,mpi_min,sml_comm,ierror)
     if(sml_mype==0) print *, 'col_3_decomp enter : ', idum1, idum2, idum0

 call check_point('col_3_decomp : begins')
#endif
  allocate(dest_count(0:col_3_total_subg_num-1), dest_count_cumul(0:col_3_total_subg_num-1), &
           dest_count_r(0:col_3_nset_allocated*sml_pe_per_plane-1), rgn_indx(sp%num)) !local

  ! 1. counting particles to be moved
  dest_count = 0
  dest_count_cumul = 0

#ifdef COL_3D_DEBUG
  allocate(myyard(0:col_3_nset_allocated-1))
  myyard = 0
#endif
  rgn_indx = -1 !This means if rgn_indx is -1, it is located in outside of domain(psi>sml_outpsi)
  do i=1, sp%num
     r=sp%ptl(i)%ph(1)
     z=sp%ptl(i)%ph(2)
     psi=psi_interpol(r,z,0,0)

     if(psi .gt. sml_outpsi .or. psi .lt. sml_inpsi) cycle

     r_cyl = sqrt((r-eq_axis_r)*(r-eq_axis_r)+(z-eq_axis_z)*(z-eq_axis_z))
     theta= acos((r-eq_axis_r)/r_cyl)
     if(z .lt. eq_axis_z) then
         theta= sml_2pi-theta
     endif

     indx = col_3_search_indx_subg(r,z,psi,theta,grid,2)

#ifdef COL_3D_DEBUG
     !if(indx .ge. col_3_total_subg_num .or. indx .lt. 0) print *, sml_mype, sml_plane_mype, floor((psi-sml_inpsi)/col_3_dpsi_core), floor(theta/col_3_dtheta_core), indx, col_3_total_subg_num, psi, sml_outpsi
     call assert(indx .lt. col_3_total_subg_num .and. indx .ge. 0,'col_3 index problem',ierror)
#endif

     rgn_indx(i) = indx
     if(indx .eq. -1) cycle
     
     dest_count(indx) = dest_count(indx) + 1
     !filter out paraticles in my yard
     do j=0, col_3_nset_allocated-1
         if(indx .eq. col_3_tag(j)) then
              dest_count(indx) = dest_count(indx)-1
#ifdef COL_3D_DEBUG
              myyard(j) = myyard(j)+1  !myyard does not count -1 region
#endif
              exit
         endif
     enddo
  enddo
  do i=0,col_3_total_subg_num-2
      dest_count_cumul(i+1)=dest_count_cumul(i)+dest_count(i)
  enddo
 
  ! 2. mpi_alltoallv countings
  dest_count_r = 0
  call mpi_alltoallv(dest_count, col_3d_scounts, col_3d_sdispls, mpi_integer, &
                     dest_count_r, col_3d_rcounts, col_3d_rdispls, mpi_integer, &
                     sml_plane_comm, ierror) 
  n_recv_total = sum(dest_count_r)
  n_send_total = sum(dest_count)

!if(sml_intpl_mype .eq. 0) then
!  open(unit=1000+sml_mype,status='replace')
!  write(1000+sml_mype,*) dest_count
!  write(1000+sml_mype,*) '-------------------'
!  write(1000+sml_mype,*) dest_count_r
!  close(1000+sml_mype)
!endif

#ifdef COL_3D_DEBUG
  allocate(grid_ptl_num(0:col_3_nset_allocated-1))

!  grid_ptl_num = 0
!  do i=0, col_3_nset_allocated-1
!      !grid_ptl_num
!      do j=0, sml_plane_totalpe-1
!          grid_ptl_num(i) = grid_ptl_num(i) + dest_count_r(j*col_3_nset_allocated+i)
!      enddo
!  enddo
!  grid_ptl_num = grid_ptl_num + myyard
!  
!  !result, i.e. after
!  write(mype_num,"(I4.4)") sml_mype 
!  open(unit=789,file='col3_decomp_'//trim(mype_num)//'.txt',status='replace')
!  do i=0, col_3_nset_allocated-1
!     write(789,987) col_3_tag(i), grid_ptl_num(i), col_3_vol_r(col_3_tag(i))
!  enddo
!  close(789)
!987 format(2(I8,1x),1(e19.13,1x))

  deallocate(grid_ptl_num, myyard)
#endif

  if(sp%num+n_recv_total-n_send_total > sp%maxnum) then
      write(*,899) sml_mype, sml_plane_mype, sp%num+n_recv_total-n_send_total, sp%maxnum
899 format('not enough ptl array size in collision3 due to col_3_decomp :', 2(I8,1x),':',2(I8,1x))
      allocate(grid_psi_norm(0:col_3_nset_allocated-1), grid_theta(0:col_3_nset_allocated-1), &
               grid_ir(0:col_3_nset_allocated-1), grid_itheta(0:col_3_nset_allocated-1))
      do i=0, col_3_nset_allocated-1
          if(col_3_tag(i) .lt. col_3_npsi*col_3_ntheta) then
              grid_ir(i) = col_3_tag(i)/col_3_ntheta            !Index starts from 0-base
              grid_itheta(i) = mod(col_3_tag(i),col_3_ntheta)   !Note that this index starts from col_3_theta_offset
              grid_psi_norm(i) = grid_ir(i)*col_3_dpsi_core/eq_x_psi
              grid_theta(i) = mod(grid_itheta(i)*col_3_dtheta_core-col_3_theta_offset+sml_2pi, sml_2pi) !not offset base, but outboard midplane base
          else
              grid_ir(i) = (col_3_tag(i)-col_3_npsi*col_3_ntheta)/col_3_ntheta_sol
              grid_itheta(i) = mod(col_3_tag(i)-col_3_npsi*col_3_ntheta,col_3_ntheta_sol)   !Note that this index starts from col_3_theta_offset
              grid_psi_norm(i) = grid_ir(i)*col_3_dpsi_sol/eq_x_psi + 1D0
              grid_theta(i) = mod(grid_itheta(i)*col_3_dtheta_sol-col_3_theta_offset+sml_2pi, sml_2pi) !not offset base, but outboard midplane base
              grid_ir(i) = grid_ir(i) + col_3_npsi
          endif
!          write(*,900)  sml_mype, sml_plane_mype, i,grid_ir(i), grid_itheta(i), grid_psi_norm(i),grid_theta(i) 
      enddo
900 format('not enough ptl array size in collision3 due to col_3_decomp :',3(I8,1x),':',2(I8,1x),2(e19.13,1x)) 
      deallocate(grid_ir, grid_itheta, grid_psi_norm, grid_theta)
      call MPI_ABORT(sml_comm, 1, ierror)
  endif

  allocate(hole_list(0:n_send_total-1))
  hole_list = -1 !initialize
  !call check_point('col_3_decomp : 1st mpi_alltoallv is okay')


  ! 3. Packing particles and preaper to send and recieve
  allocate(ptl_send(0:n_send_total-1), ptl_recv(0:n_recv_total-1))
  allocate(scounts_ptl(0:sml_pe_per_plane-1), rcounts_ptl(0:sml_pe_per_plane-1), &
           sdispls_ptl(0:sml_pe_per_plane-1), rdispls_ptl(0:sml_pe_per_plane-1))
  ! scounts, rcounts
  scounts_ptl(:) = 0
  rcounts_ptl(:) = 0
  do i=0, sml_pe_per_plane-2
      do j=col_3d_sdispls(i), col_3d_sdispls(i+1)-1
          scounts_ptl(i) = scounts_ptl(i)+dest_count(j)
      enddo
      do j=col_3d_rdispls(i), col_3d_rdispls(i+1)-1
          rcounts_ptl(i) = rcounts_ptl(i)+dest_count_r(j)
      enddo
  enddo
  scounts_ptl(sml_pe_per_plane-1) = sum(dest_count(col_3_total_subg_num-col_3_nset:col_3_total_subg_num-1))
  rcounts_ptl(sml_pe_per_plane-1) = sum(dest_count_r(col_3_nset_allocated*(sml_pe_per_plane-1):col_3_nset_allocated*sml_pe_per_plane-1))
 
!if(sml_intpl_mype .eq. 3) then
!  open(unit=1100+sml_plane_mype,status='replace')
!  write(1100+sml_plane_mype,*) scounts_ptl
!  write(1100+sml_plane_mype,*) '-------------------'
!  write(1100+sml_plane_mype,*) rcounts_ptl
!  close(1100+sml_plane_mype)
!endif
  
  ! sdispls
  sdispls_ptl(:) = 0
  rdispls_ptl(:) = 0
  do i=0, sml_pe_per_plane-2
     sdispls_ptl(i+1) = sdispls_ptl(i) + scounts_ptl(i)
     rdispls_ptl(i+1) = rdispls_ptl(i) + rcounts_ptl(i)
  enddo

  htotal = 0 !index for holes
  do i=1, sp%num
   
     !in order to filter out my particles (just for sure)
     snd_logic = .true.  
     indx = rgn_indx(i)
     if(indx .eq. -1) cycle   !psi>sml_outpsi 
     do j=0, col_3_nset_allocated-1
         if( indx .eq. col_3_tag(j) ) then
             snd_logic = .false.
             exit
         endif
     enddo
     
     if(snd_logic) then 
         dindx=dest_count_cumul(indx)
         ptl_send(dindx)%ph = sp%ptl(i)%ph
         ptl_send(dindx)%ct = sp%ptl(i)%ct
         ptl_send(dindx)%phase0 = sp%phase0(:,i)
         ptl_send(dindx)%gid = sp%ptl(i)%gid

         dest_count_cumul(indx) = dindx+1
         hole_list(htotal) = i
         htotal = htotal+1
     endif
  enddo

  ! 4. MPI_ALLtoALLv
  call mpi_alltoallv(ptl_send, scounts_ptl, sdispls_ptl, mpi_ptl_type, ptl_recv, rcounts_ptl, rdispls_ptl, mpi_ptl_type, sml_plane_comm,ierror) 

!if(sml_intpl_mype .eq. 0) then
!  !real (8) :: ph(ptl_nphase) !6
!  !real (8) :: ct(ptl_nconst) !3
!  !real (8) :: phase0(ptl_nphase) !6
!  open(unit=1200+sml_plane_mype,position='append')
!  do i=0, n_send_total-1
!      write(1200+sml_plane_mype,1201) i, ptl_send(i)%ph, ptl_send(i)%ct, ptl_send(i)%phase0, ptl_send(i)%gid
!  enddo
!  close(1200+sml_plane_mype)
!1201 format(I8,1x, 15(e19.13,1x), I8)
!
!  open(unit=1300+sml_plane_mype,position='append')
!  do i=0, n_recv_total-1
!      write(1300+sml_plane_mype,1201) i, ptl_recv(i)%ph, ptl_recv(i)%ct, ptl_recv(i)%phase0, ptl_recv(i)%gid
!  enddo
!  close(1300+sml_plane_mype)
!endif


  ! Up to now, decomposition has been done to collision cell level. However, since 1 CPU deal with many collision cells, we don't need this.
  ! Therefore, for fast code, we can modify above routine to reduce MPI communication. 
  !!! 5. Unpacking recieved particles - hole treatment
  findx = 0

  do i=0, n_recv_total-1
     if(i .lt. htotal) then
         fhindx = hole_list(findx)
         hole_list(findx) = -1  !this implies "filled"
         findx = findx+1
     else
         fhindx = sp%num+1+i-findx 
     endif
     sp%ptl(fhindx)%ph = ptl_recv(i)%ph
     sp%ptl(fhindx)%ct = ptl_recv(i)%ct
     sp%phase0(:,fhindx) = ptl_recv(i)%phase0
     sp%ptl(fhindx)%gid = ptl_recv(i)%gid
  enddo
  !call check_point('check filling-hole')

  ! Below is not best algorithm as you can see.
  ! There should be better altorithm than this!!! Please simplify this part if you can! -ES
  if(n_recv_total .lt. htotal) then
      s=sp%num
      m=htotal-1
      do i=0, htotal-n_recv_total-1
          !sp%num-i => hole_list
          fhindx = hole_list(findx)
          if(s .eq. fhindx) exit   !This is a case where hole exists at the end of ptl array. So exit.

          ! s: beginning points for existing particles
          ! j: particles to be put into holes
          !     fhindx(hole) <=> j
          j=s
          do
#ifdef COL_3D_DEBUG
             if( j .lt. 0 .or. hole_list(m) .eq. -1) then
                 print *, 'hole filling : infinite loops in col_3_decomp : logical error', sml_mype, j, m, findx, fhindx, s
                 call MPI_ABORT(sml_comm, 1, ierror)
             endif
#endif
             if(hole_list(m) .ne. j .or. j .le. fhindx) then
                 s = j-1  !next starting point becomes right below the j
                 exit
             else
                 m=m-1
                 j=j-1
             endif
          enddo
#ifdef COL_3D_DEBUG
          if( j .eq. sp%num-htotal-2) then
              print *, 'hole filling logic error'
              call MPI_ABORT(sml_comm, 1, ierror)
          endif
#endif
          if(j .le. fhindx) exit

          hole_list(findx) = j  !j becomes holes now
          findx = findx+1

          sp%ptl(fhindx)%ph = sp%ptl(j)%ph
          sp%ptl(fhindx)%ct = sp%ptl(j)%ct
          sp%phase0(:,fhindx) = sp%phase0(:,j)
          sp%ptl(fhindx)%gid = sp%ptl(j)%gid
      enddo
  endif

  sp%num = sp%num - htotal + n_recv_total  !htotal = n_send_total
  !if(sml_mype .eq. 0) print *, sml_mype, 'ptl num', sp%num, n_send_total, n_recv_total, sp%type
#ifdef COL_3D_DEBUG
  call assert( htotal .eq. n_send_total, 'col_3 logical error at fill-hole', ierror)
#endif


#ifdef COL_3D_DEBUG
  ! 6. check if all paraticles are in our regions
  do i=1,sp%num
     r=sp%ptl(i)%ph(1)
     z=sp%ptl(i)%ph(2)
     psi=psi_interpol(r,z,0,0)
     if(psi .gt. sml_outpsi .or. psi .lt. sml_inpsi) cycle
     r_cyl = sqrt((r-eq_axis_r)*(r-eq_axis_r)+(z-eq_axis_z)*(z-eq_axis_z))
     theta= acos((r-eq_axis_r)/r_cyl)
     if(z .lt. eq_axis_z) then
         theta= sml_2pi-theta
     endif

     indx = col_3_search_indx_subg(r,z,psi,theta,grid,2)

     snd_logic = .true.
     
     do j=0, col_3_nset_allocated-1
         if( indx .eq. col_3_tag(j) ) then
             snd_logic = .false.
             exit
         endif
     enddo

     if(snd_logic) then
         print *, 'filling error', sml_mype, sml_plane_mype, i, indx, psi, r, z
     endif

  enddo
#endif

  deallocate(hole_list,rgn_indx) 
  deallocate(ptl_send, ptl_recv)
  deallocate(scounts_ptl, rcounts_ptl, sdispls_ptl, rdispls_ptl)
  deallocate(dest_count, dest_count_cumul, dest_count_r)
  
   
#ifdef COL_3D_DEBUG
  call check_point('col_3_decomp : moving is completed')

     idum3=sp%num
     call mpi_allreduce(idum3,idum1,1,mpi_integer8,mpi_sum,sml_comm,ierror)
     call mpi_allreduce(idum3,idum2,1,mpi_integer8,mpi_max,sml_comm,ierror)
     call mpi_allreduce(idum3,idum0,1,mpi_integer8,mpi_min,sml_comm,ierror)
     if(sml_mype==0) print *, 'col_3_decomp exit : ', idum1, idum2, idum0
#endif

     idum3=sp%num
     call mpi_allgather(idum3, 1, mpi_integer8, ptl_popul, 1, mpi_integer8, sml_comm, ierror)
     idum2 = maxloc(ptl_popul,1)
     if(sml_mype .eq. idum2-1) then
         open(unit=702,status='replace')
         write(702,*) col_3_tag
         close(702)
     endif
end subroutine col_3_decomp

!Note that original input FPL_rx has been changed to Dlx
!2013-02-20 : '_m' implies use for multi-species
!2013-03-01 : Merging all operations 
subroutine col_3_core_m(spi, Dlxi, mesh_dri, mesh_dzi, dist_ni, voli, &
                        spe, Dlxe, mesh_dre, mesh_dze, dist_ne, vole, &
                        vpic_gamma, vpic_ierr)
    use sml_module, only : sml_e_charge, sml_pi, sml_mype, sml_plane_mype
    use ptl_module 
    use col_module
    implicit none
    type(species_type) :: spi, spe 
    real (kind=8) :: Dlxi, mesh_dzi, mesh_dri  !original code : -FPL_rx->Dlx
    real (kind=8) :: Dlxe, mesh_dze, mesh_dre  !original code : -FPL_rx->Dlx
    real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_ni, dist_iteri   ! local 
    real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_ne, dist_itere   ! local 
    real (kind=8), dimension(1:col_3_nvr) :: voli,vole      ! local
    real (kind=8), dimension (col_3_ntotal_v) :: dist_coli, dist_n_coli
    real (kind=8), dimension (col_3_ntotal_v) :: dist_cole, dist_n_cole
    real (kind=8) :: vpic_gamma(4)  ! 1: i-i, 2: i-e, 3: e-i, 4: e-e
    integer :: vpic_ierr


    integer :: mesh_Nr, mesh_Nz
    integer :: vpic_inner_iter_max, mesh_Nrm1, mesh_Nzm1, mesh_Nrzm1
    real (kind=8) :: vpic_tstep, vpic_tstep_gamma(4) 
    real (kind=8), dimension (:,:,:,:), allocatable :: EDi, EDe
    real (kind=8), dimension (:,:,:), allocatable ::fi_half, fe_half, dfidr, dfidz, dfedr, dfedz 
    real (kind=8), dimension (:,:,:), allocatable :: M_i, M_e 
    real (kind=8), dimension (:,:,:), allocatable :: M_ie, M_ei 
    real (kind=8) :: numeric_T, vpic_mass, tmpr1, tmpr2, tmpr3
    integer :: index_ip, index_jp, itmp, index_I, index_J
    integer :: iter_inter
    real (kind=8) :: vpic_exit_dni, vpic_exit_dwi, vpic_exit_n_previ, vpic_exit_w_previ, &
                     vpic_exit_dne, vpic_exit_dwe, vpic_exit_n_preve, vpic_exit_w_preve, &
                     vpic_exit_dw_w, vpic_exit_dw_w_prev, vpic_exit_dwi_wi, vpic_exit_dwe_we, vpic_exit_dni_ni, vpic_exit_dne_ne, &
                     vpic_exit_dw_rel, vpic_exit_dni_rel, vpic_exit_dne_rel 
    ! Super LU Variables
    real (kind=8), dimension(LU_nnz) :: LU_valuesi, LU_valuese, LU_values_tmp 
    !-------------------------------------------------------------Super LU
    integer :: i, j, m, ierr
    type(col_3_core_type) :: ci, ce

    mesh_Nr = col_3_nvr
    mesh_Nz = col_3_nvz

    mesh_Nrm1 = mesh_Nr-1
    mesh_Nzm1 = mesh_Nz-1
    mesh_Nrzm1 = mesh_Nrm1*mesh_Nzm1

    vpic_inner_iter_max = 100
    vpic_tstep = col_3_dt
    vpic_tstep_gamma = vpic_tstep*vpic_gamma 

    vpic_ierr = 0

    call col_3_core_init(spi, ci, dist_ni, mesh_dri, mesh_dzi, dlxi, voli, &
                         spe, ce, dist_ne, mesh_dre, mesh_dze, dlxe, vole )

    allocate(M_ie(5,mesh_Nrzm1,mesh_Nrzm1),M_ei(5,mesh_Nrzm1,mesh_Nrzm1)) ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za 
    allocate(M_i(4,mesh_Nrzm1,mesh_Nrzm1), M_e(4,mesh_Nrzm1,mesh_Nrzm1) ) ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra 
    M_ie = 0D0
    M_ei = 0D0
    M_i = 0D0
    M_e = 0D0
    call col_3_angle_avg_s(ci, M_i)
    call col_3_angle_avg_s(ce, M_e)
    call col_3_angle_avg_1st(ci, ce, M_ie)
    call col_3_angle_avg_2nd(ci, ce, M_ie, M_ei)
           
    dist_iteri = dist_ni   !dist_iter <-implicitly updated distribution  (2D)
    dist_n_coli = reshape(dist_ni, (/mesh_Nr*mesh_Nz/))  !dist_n_col <- saved nth column distribution (1D) 
    dist_itere = dist_ne   !dist_iter <-implicitly updated distribution  (2D)
    dist_n_cole = reshape(dist_ne, (/mesh_Nr*mesh_Nz/))  !dist_n_col <- saved nth column distribution (1D) 

    !NOTE THAT "dist_iter" is the PDF iterated and being updated in this LOOP
    do iter_inter=1, vpic_inner_iter_max

       dist_coli = dist_n_coli                                 !dist_col SHOULD BE always "f_n"
       dist_cole = dist_n_cole                                 !dist_col SHOULD BE always "f_n"



       !f_(I+1/2,J+1/2) = f_half, dfdr, dfdz
       !To be deallocated after col_3_E_and_D_m in this subroutine
       allocate( fi_half(2,mesh_Nrm1, mesh_Nzm1), dfidr(2, mesh_Nrm1, mesh_Nzm1), dfidz(2,mesh_Nrm1, mesh_Nzm1), &
                 fe_half(2,mesh_Nrm1, mesh_Nzm1), dfedr(2, mesh_Nrm1, mesh_Nzm1), dfedz(2,mesh_Nrm1, mesh_Nzm1) )
       call col_3_f_df(ci%delta_r(1,:), ci%delta_z(1,:), ci%mesh_dr, ci%mesh_dz, dist_iteri, fi_half(1,:,:), dfidr(1,:,:), dfidz(1,:,:))  !i-i
       call col_3_f_df(ci%delta_r(2,:), ci%delta_z(2,:), ci%mesh_dr, ci%mesh_dz, dist_iteri, fi_half(2,:,:), dfidr(2,:,:), dfidz(2,:,:))  !i-e
       call col_3_f_df(ce%delta_r(1,:), ce%delta_z(1,:), ce%mesh_dr, ce%mesh_dz, dist_itere, fe_half(1,:,:), dfedr(1,:,:), dfedz(1,:,:))  !e-e
       call col_3_f_df(ce%delta_r(2,:), ce%delta_z(2,:), ce%mesh_dr, ce%mesh_dz, dist_itere, fe_half(2,:,:), dfedr(2,:,:), dfedz(2,:,:))  !e-i


       ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
       !To be deallocated after col_3_LU_matrix in this subroutine 
       allocate(EDi(6,mesh_Nrm1, mesh_Nzm1,2), EDe(6,mesh_Nrm1, mesh_Nzm1,2))
       call col_3_E_and_D_s(ci, ci, fi_half(1,:,:), dfidr(1,:,:), dfidz(1,:,:), M_i, EDi(:,:,:,1)) ! i-i 
       call col_3_E_and_D_s(ce, ce, fe_half(1,:,:), dfedr(1,:,:), dfedz(1,:,:), M_e, EDe(:,:,:,1)) ! e-e
       call col_3_E_and_D_m(ci, ce, fe_half(2,:,:), dfedr(2,:,:), dfedz(2,:,:), M_ie, EDe(:,:,:,2)) ! i-e
       call col_3_E_and_D_m(ce, ci, fi_half(2,:,:), dfidr(2,:,:), dfidz(2,:,:), M_ei, EDi(:,:,:,2)) ! e-i

       deallocate( fi_half, dfidr, dfidz, fe_half, dfedr, dfedz )


       ! LU_valuesi, LU_valuese, dist_iteri, dist_itere
       call col_3_LU_matrix(1, ci, EDi(:,:,:,1), iter_inter, LU_values_tmp)  ! i-i
       LU_valuesi = LU_values_tmp*vpic_tstep_gamma(1)
       call col_3_LU_matrix(2, ci, EDe(:,:,:,2), iter_inter, LU_values_tmp)  ! i-e
       LU_valuesi = LU_valuesi + LU_values_tmp*vpic_tstep_gamma(2)
       !LU_valuesi = LU_values_tmp*vpic_tstep_gamma(2)   !For debugging
       call col_3_LU_matrix(2, ce, EDi(:,:,:,2), iter_inter, LU_values_tmp)  ! e-i
       LU_valuese = LU_values_tmp*vpic_tstep_gamma(3)
       call col_3_LU_matrix(1, ce, EDe(:,:,:,1), iter_inter, LU_values_tmp)  ! e-e
       LU_valuese = LU_valuese + LU_values_tmp*vpic_tstep_gamma(4)
       !LU_valuese = LU_values_tmp*vpic_tstep_gamma(4)  !For debugging

       
       deallocate(EDi, EDe)

       call col_3_picard_step(iter_inter, LU_valuesi, dist_coli, dist_iteri)
       call col_3_picard_step(iter_inter, LU_valuese, dist_cole, dist_itere)


       call col_3_convergence_eval(ci, dist_ni, dist_iteri, vpic_exit_dni, vpic_exit_dwi, vpic_exit_n_previ, vpic_exit_w_previ)
       call col_3_convergence_eval(ce, dist_ne, dist_itere, vpic_exit_dne, vpic_exit_dwe, vpic_exit_n_preve, vpic_exit_w_preve)

       vpic_exit_dni_ni  = dabs(vpic_exit_dni/ci%dens)
       vpic_exit_dne_ne  = dabs(vpic_exit_dne/ce%dens)
       vpic_exit_dw_w    = dabs((vpic_exit_dwi+vpic_exit_dwe)/(ci%ens+ce%ens))
       vpic_exit_dwi_wi  = dabs(vpic_exit_dwi/ci%ens) !For debugging
       vpic_exit_dwe_we  = dabs(vpic_exit_dwe/ce%ens) !For debugging
       vpic_exit_dni_rel = dabs(vpic_exit_dni/vpic_exit_n_previ)
       vpic_exit_dne_rel = dabs(vpic_exit_dne/vpic_exit_n_preve)
       vpic_exit_dw_rel  = dabs(vpic_exit_dw_w/vpic_exit_dw_w_prev)

       if(sml_mype .eq. 0) then
           write(*,1001) vpic_exit_dni_ni, vpic_exit_dne_ne, vpic_exit_dw_w
1001 format('(dn/n)_i=', e19.13, 1x, '(dn/n)_e=', e19.13,1x, 'dw/w=', e19.13,1x)  
       endif

       ! DIST_ITER has new updated PDF.
       if( vpic_exit_dni_ni .le. 1D-10 .and. vpic_exit_dne_ne .le. 1D-10 &
           .and. vpic_exit_dw_w .le. 1D-10 .and. iter_inter .ne. 1 ) then
           ! Accuracy is enough to go
           exit
       else if(      vpic_exit_dni_rel .lt. 1D-1 .and. vpic_exit_dne_rel .lt. 1D-1 .and. vpic_exit_dw_rel .lt. 1D-1 &
               .and. vpic_exit_dni_ni  .lt. 1D-5 .and. vpic_exit_dne_ne  .lt. 1D-5 .and. vpic_exit_dw_w .lt. 1D-5 &
               .and. iter_inter .ne. 1 ) then
                    print *, 'col_3 picard iteration - slow convergence exit', &
                    sml_mype, sml_plane_mype, vpic_exit_dni_rel, vpic_exit_dne_rel, vpic_exit_dw_rel, &
                    vpic_exit_dni_ni, vpic_exit_dne_ne, vpic_exit_dw_w 
                    exit
       else
           ! FOR ITERATION
           if( iter_inter .eq. 50 ) print *, sml_mype,': iteration is over 50', & 
                    sml_mype, sml_plane_mype, vpic_exit_dni_rel, vpic_exit_dne_rel, vpic_exit_dw_rel, &
                    vpic_exit_dni_ni, vpic_exit_dne_ne, vpic_exit_dw_w
        
           vpic_exit_w_previ = vpic_exit_dwi
           vpic_exit_w_preve = vpic_exit_dwe
           vpic_exit_dw_w_prev = vpic_exit_dw_w 
       endif


    enddo !iter_inter

    if (iter_inter .eq. vpic_inner_iter_max) then
        print *, 'inner iteration went to maximum number at :', sml_mype, sml_plane_mype, & 
                    sml_mype, sml_plane_mype, vpic_exit_dni_rel, vpic_exit_dne_rel, vpic_exit_dw_rel, &
                    vpic_exit_dni_ni, vpic_exit_dne_ne, vpic_exit_dw_w 
    endif

    dist_ni = dist_iteri 
    dist_ne = dist_itere 

    deallocate(M_i, M_e, M_ie, M_ei)
    call col_3_core_deallocation(ci ,ce) 

end subroutine col_3_core_m

subroutine col_3_core_deallocation(ci, ce)
    use col_module
    implicit none
    type(col_3_core_type) :: ci, ce

    deallocate(ci%mesh_r, ci%mesh_z, ci%mesh_r_half, ci%mesh_z_half, ci%local_center_volume, ci%vol)
    deallocate(ce%mesh_r, ce%mesh_z, ce%mesh_r_half, ce%mesh_z_half, ce%local_center_volume, ce%vol)
end subroutine col_3_core_deallocation

subroutine col_3_angle_avg_1st(c1, c2, M_ie)
    use sml_module, only : sml_pi, sml_mype
    use col_module
    implicit none
    type(col_3_core_type) :: c1, c2
    real (8), dimension(5, (col_3_nvr-1)*(col_3_nvz-1), (col_3_nvr-1)*(col_3_nvz-1)) :: M_ie 
    integer :: index_I, index_ip, index_J, index_jp
    integer :: index_rz, index_rc, index_ac, index_az
    real (8) :: r, z, a, c, dz, dz2, r2, a2, lambda, k, k_eff, kp1_sqrt, km1, EK, EE
    integer :: vpic_ierr
    real (8) :: I1, I2, I3, temp_cons, temp_cons_dz, temp_cons_I1, temp_cons_I2, temp_cons_I3, temp_cons_dz_I2, temp_cons_dz_I3
    real (8) :: tmp_a(5)
    integer :: mesh_Nrm1, mesh_Nzm1

    mesh_Nrm1 = col_3_nvr-1
    mesh_Nzm1 = col_3_nvz-1

    ! new routine
    ! (r,z) : colliding paritcle "c1" - capital index, (a,c) : target particle "c2" - small index
    
    do index_I=1,mesh_Nzm1
        z = c1%mesh_z_half(index_I)

        do index_J=1,mesh_Nrm1
            index_rz = index_J +mesh_Nrm1*(index_I-1)
            r=c1%mesh_r_half(index_J)
            r2 = r*r

            do index_ip=1, mesh_Nzm1
                c = c2%mesh_z_half(index_ip)
                dz = z-c
                dz2 = dz*dz

                do index_jp=1, mesh_Nrm1
                    a=c2%mesh_r_half(index_jp)
                    a2 = a*a
                    index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                    !lambda
                    lambda = r2+a2+dz2     !symmetric (r<->a, z<->c) 
                                           ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                           !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
                                           
                    k = 2D0*r*a/lambda     !symmetric
                    k_eff = 2D0*k/(1D0+k)
                    kp1_sqrt = sqrt(1D0+k)
                    km1 = k-1D0

                    if(k_eff * 1.1D0 .eq. k_eff) then
                        print *, 'Nan appears', sml_mype, k_eff, r, z, a, c, c1%mesh_z_half
                        stop
                    endif

                    call elliptics(k_eff,EK,EE,vpic_ierr)

                    !Calulation of M coeff. (all symmetric)
                    I1 = -4D0*((1D0+k)*EE-EK)/(k*k*kp1_sqrt)
                    I2 = -2D0*EE/(km1*kp1_sqrt)
                    I3 = -2D0*(EE+km1*EK)/(km1*k*kp1_sqrt )
                    !I4 = I2-I1 !For exact Numerical Conservation, and It should be mathematically
                    temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda))
                    temp_cons_dz = temp_cons*dz
                    temp_cons_I1 = temp_cons*I1
                    temp_cons_I2 = temp_cons*I2
                    temp_cons_I3 = temp_cons*I3
                    temp_cons_dz_I2 = temp_cons_dz*I2
                    temp_cons_dz_I3 = temp_cons_dz*I3


                    ! index_rz = index_J +mesh_Nrm1*index_I
                    ! index_ac = index_jp+mesh_Nrm1*index_ip
                    ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za
                    tmp_a(1) = temp_cons_I1*a2+temp_cons_I2*dz2    
                    tmp_a(2) = temp_cons_dz_I3*a-temp_cons_dz_I2*r
                    tmp_a(3) = temp_cons_I2*(r2+a2) - temp_cons_I3*2D0*r*a 
                    tmp_a(4) = temp_cons_I1*r*a + temp_cons_I3*dz2 
                    tmp_a(5) = temp_cons_dz_I2*a-temp_cons_dz_I3*r 

                    ! index order : (target, colliding) - What about M_zc?
                    M_ie(:,index_ac, index_rz) = tmp_a  
                    
 !=============CONSERVATION
 
 
 
 !=========================                   
                     
                enddo !index_jp
            enddo !index_J
        enddo !index_ip
    enddo ! index_I

end subroutine col_3_angle_avg_1st


subroutine col_3_angle_avg_2nd(c1, c2, M1, M2 )
    use sml_module, only : sml_pi
    use col_module
    implicit none
    type(col_3_core_type) :: c1, c2
    real (8), dimension(5, (col_3_nvr-1)*(col_3_nvz-1), (col_3_nvr-1)*(col_3_nvz-1)) :: M1, M2  
    integer :: index_I, index_ip, index_J, index_jp
    integer :: index_rz, index_rc, index_ac, index_az
    real (8) :: r, z, a, c, dz, dz2, r2, a2, lambda, k, k_eff, kp1_sqrt, km1, EK, EE
    integer :: vpic_ierr
    real (8) :: I1, I2, temp_cons
    integer :: mesh_Nrm1, mesh_Nzm1

    mesh_Nrm1 = col_3_nvr-1
    mesh_Nzm1 = col_3_nvz-1
                    
    !M1 : M_ie, M2 : M_ei 
    ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za
    M2(2,:,:) = transpose(M1(5,:,:)) ! M_rz(v,v')=M_za(v',v)
    M2(5,:,:) = transpose(M1(2,:,:)) ! M_za(v,v')=M_rz(v',v)
    M2(3,:,:) = transpose(M1(3,:,:))   !Symmetric
    M2(4,:,:) = transpose(M1(4,:,:))   !Symmetric

    !c2->c1
    do index_I=1,mesh_Nzm1
        z = c2%mesh_z_half(index_I)

        do index_J=1,mesh_Nrm1
            index_rz = index_J +mesh_Nrm1*(index_I-1)
            r=c2%mesh_r_half(index_J)
            r2 = r*r

            do index_ip=1, mesh_Nzm1
                c = c1%mesh_z_half(index_ip)
                dz = z-c
                dz2 = dz*dz

                do index_jp=1, mesh_Nrm1
                    a=c1%mesh_r_half(index_jp)
                    a2 = a*a
                    index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                    !lambda
                    lambda = r2+a2+dz2     !symmetric (r<->a, z<->c) 
                                           ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                           !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
                                           
                    k = 2D0*r*a/lambda     !symmetric
                    k_eff = 2D0*k/(1D0+k)
                    kp1_sqrt = sqrt(1D0+k)
                    km1 = k-1D0

                    call elliptics(k_eff,EK,EE,vpic_ierr)

                    !Calulation of M coeff. (all symmetric)
                    I1 = -4D0*((1D0+k)*EE-EK)/(k*k*kp1_sqrt)
                    I2 = -2D0*EE/(km1*kp1_sqrt)
                    temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda))

                    ! index order : (target, colliding) 
                    M2(1,index_ac, index_rz) = temp_cons*I1*a2+temp_cons*I2*dz2
                enddo !index_jp
            enddo !index_J
        enddo !index_ip
    enddo ! index_I

end subroutine col_3_angle_avg_2nd


subroutine col_3_convergence_eval(cs, dist_n, dist_iter, vpic_dn, vpic_dw, vpic_n_prev, vpic_w_prev)
       use col_module
       implicit none
       type(col_3_core_type) :: cs
       real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_n, dist_iter   ! local 
       real (kind=8) :: vpic_exit_num, vpic_exit_en 

       real (8) :: vpic_dn, vpic_dw, vpic_n_prev, vpic_w_prev, vpic_dfc
       integer :: index_I, index_J
       real (8) :: tmpr1, tmpr2, tmpr3 
       ! DIST_ITER has new updated PDF.

       vpic_dn = 0D0
       vpic_dw = 0D0
       vpic_n_prev = 0D0
       vpic_w_prev = 0D0
       do index_I=1,col_3_nvz
           tmpr1 = cs%mesh_z(index_I)
           tmpr1 = tmpr1*tmpr1     ! mesh_z^2
           do index_J=1,col_3_nvr
               tmpr2 = cs%mesh_r(index_J)
               tmpr2 = tmpr2*tmpr2     ! mesh_r^2
               tmpr3 = tmpr1+tmpr2     ! mesh_z^2+mesh_r^2
               vpic_dfc = (dist_iter(index_J,index_I) - dist_n(index_J, index_I))*cs%vol(index_J)

               vpic_dn = vpic_dn + vpic_dfc
               vpic_dw = vpic_dw + vpic_dfc * tmpr3
               vpic_n_prev = vpic_n_prev + dist_n(index_J, index_I)*cs%vol(index_J) 
               vpic_w_prev = vpic_w_prev + dist_n(index_J,index_I)*tmpr3*cs%vol(index_J)
           enddo
       enddo

       vpic_dw     = vpic_dw*cs%mass
       vpic_w_prev = vpic_w_prev*cs%mass
#ifdef COL_3_CORE_MSG
       if(sml_mype .eq. 0) then
           print *, 'dn/n = ', vpic_dn/cs%dens, 'dw/w = ', vpic_dw/vpic_dw_prev 
       endif
#endif       
end subroutine col_3_convergence_eval

subroutine col_3_picard_step(iter_inter, LU_values, dist_col, dist_iter)
    use sml_module, only : sml_mype
    use col_module
    implicit none
    integer :: iter_inter
    real (kind=8), dimension(LU_nnz) :: LU_values
    real (kind=8), dimension (col_3_ntotal_v) :: dist_col
    real (kind=8), dimension(col_3_nvr,col_3_nvz) :: dist_iter
    integer :: LU_factors(8)
    integer :: LU_info, LU_info2
    integer :: index_I, index_J, mat_pos, mat_pos_d, local_ij, local_i, local_j

    if( iter_inter .eq. 1 ) then
        dist_iter = 0D0
        do index_I=1,col_3_nvz
            do index_J=1, col_3_nvr
               mat_pos = index_J+(index_I-1)*col_3_nvr

               !explicit time marching for 1st guess
               do local_ij=1,9
                  local_i = (local_ij-1)/3 - 1
                  local_j = mod(local_ij-1,3) - 1
                  if( index_J+local_j .gt. 0 .and. index_I+local_i .gt. 0 &
                      .and. index_J+local_j .lt. col_3_nvr+1 .and. index_I+local_i .lt. col_3_nvz+1 ) then
                      !valid mesh

                      mat_pos_d = (index_J+local_j)+((index_I+local_i)-1)*col_3_nvr
                      if( local_ij .eq. 5 ) then
                          ! (I + L(f^n) dt) f^n  : diagonal part
                          dist_iter(index_J, index_I) = dist_iter(index_J, index_I) + (1D0 + LU_values(LU_Cvalues(mat_pos_d))) * dist_col(mat_pos_d)
                          !print *, LU_values(LU_Cvalues(mat_pos))*vpic_tstep*vpic_gamma, vpic_tstep
                      else
                          ! L(f^n) dt f^n : off-diagonal part
                          dist_iter(index_J, index_I) = dist_iter(index_J, index_I) + LU_values(index_map_LU(local_ij,mat_pos)) * dist_col(mat_pos_d) 
                      endif
                  endif
               enddo   !local_ij
            enddo  !index_J
        enddo  !index_I
    else 
        !IMPLICIT TIME MARCHING
        ! NOTE!!!! THAT the below "-" sign has some reasons.
        ! There can be easy mistakes that I(identity) - M doesn't mean
        ! 1-diagonal part. Off-diagonal part will change to -M. Be careful!
        LU_values = -LU_values

        LU_values(LU_cvalues) = 1D0 + LU_values(LU_cvalues)


        !Then, finally solve! Super LU!! =)
        !SUPER LU!!!! - We can directly call superLU wihtout below wrapper function!
        call c_fortran_dgssv(LU_n, LU_nnz, LU_nrhs, LU_values, LU_rowindx, LU_colptr, &
                              dist_col, LU_ldb, LU_factors, LU_info, LU_info2 )



        if(LU_info .ne. 0) then
            write(*,*) 'SuperLU : Info = ',LU_info, 'mype :', sml_mype
            write(*,*) 'It is very possible that you got NAN or INF since matrix component is NAN or INF'

            if(sml_mype .eq. 23) then
                open(unit=987,file='luinfo.txt',status='replace')
                write(987,*) LU_values
                close(987)
            endif
        endif

        dist_iter = reshape(dist_col,(/col_3_nvr,col_3_nvz/))
    endif 

end subroutine col_3_picard_step

subroutine col_3_LU_matrix_ftn(op_mode, cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs_comp, LU_values)
   use col_module
   implicit none
   integer :: op_mode
   type(col_3_core_type) :: cs
   integer :: mat_pos_rel_indx(4)
   integer :: mat_pos, cell_I, cell_J
   real (8) :: coeff1, coeff2, coeff_loc1, coeff_loc2
   real (8) :: EDs_comp(6)
   real (kind=8), dimension(LU_nnz) :: LU_values 
   real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz

    coeff_loc1 = coeff1*cs%mesh_r_half(cell_J)
    coeff_loc2 = coeff2*cs%mesh_r_half(cell_J)
    ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
    M_rr_fOVERdr = EDs_comp(1)/cs%mesh_dr
    M_rz_fOVERdz = EDs_comp(2)/cs%mesh_dz
    M_zz_fOVERdz = EDs_comp(3)/cs%mesh_dz
    M_zr_fOVERdr = EDs_comp(4)/cs%mesh_dr
    tmp_ER_sum   = EDs_comp(5)
    tmp_EZ_sum   = EDs_comp(6)
    delr = cs%delta_r(op_mode, cell_J)
    cdelr = (1D0-delr)
    delz = cs%delta_z(op_mode,cell_I)
    cdelz = (1D0-delz)


    !for (I,J)
    LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) & 
                               + coeff_loc1* (tmp_ER_sum*delr*delz - M_rr_fOVERdr*(-delz) - M_rz_fOVERdz*(-delr)) &
                               + coeff_loc2* (tmp_EZ_sum*delr*delz - M_zr_fOVERdr*(-delz) - M_zz_fOVERdz*(-delr)) 
    !for (I,J+1)
    LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) &
                               + coeff_loc1* (tmp_ER_sum*cdelr*delz - M_rr_fOVERdr*(delz) - M_rz_fOVERdz*(-cdelr)) &
                               + coeff_loc2* (tmp_EZ_sum*cdelr*delz - M_zr_fOVERdr*(delz) - M_zz_fOVERdz*(-cdelr)) 
    !for(I+1,J)
    LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) &
                               + coeff_loc1* (tmp_ER_sum*delr*cdelz - M_rr_fOVERdr*(-cdelz) - M_rz_fOVERdz*(delr)) &
                               + coeff_loc2* (tmp_EZ_sum*delr*cdelz - M_zr_fOVERdr*(-cdelz) - M_zz_fOVERdz*(delr)) 
    !for(I+1,J+1)
    LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) &
                               + coeff_loc1* (tmp_ER_sum*cdelr*cdelz - M_rr_fOVERdr*(cdelz) - M_rz_fOVERdz*(cdelr)) &
                               + coeff_loc2* (tmp_EZ_sum*cdelr*cdelz - M_zr_fOVERdr*(cdelz) - M_zz_fOVERdz*(cdelr))


end subroutine col_3_LU_matrix_ftn

subroutine col_3_LU_matrix(op_mode, cs, EDs, iter_inter ,LU_values)
! EDs : target
   use sml_module, only : sml_mype
   use col_module
   implicit none
   integer :: op_mode
   type(col_3_core_type) :: cs
   real (8), dimension(6,col_3_nvr-1, col_3_nvz-1) :: EDs 
   integer :: iter_inter
   real (kind=8), dimension(LU_nnz) :: LU_values 
   integer :: index_I, index_J, cell_I, cell_J, mesh_Nr, mesh_Nz
   integer :: mat_pos_rel_indx(4)
   real (8) :: coeff1_ab, coeff2_ab, coeff1, coeff2, local_vol
   real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz
   integer :: local_ij, local_i, local_j, mat_pos, mat_pos_d 

   mesh_Nr = col_3_nvr
   mesh_Nz = col_3_nvz

   LU_values(:) = 0D0
   do index_I=1,mesh_Nz
       do index_J=1, mesh_Nr
          mat_pos = index_J+(index_I-1)*col_3_nvr
            ! What we know : the number of cells related to cell (J,I) from LU_colptr
            !                index number of relevant cells from LU_rowindx
            
            ! NOTE that local_vol is for volume in LHS. Therefore, fixed values!
            local_vol = cs%vol(index_J)
            coeff1_ab    = 0.5D0*cs%mesh_dz/local_vol
            coeff2_ab    = 0.5D0*cs%mesh_dr/local_vol

            if((index_I .ne. mesh_Nz) .and. (index_J .ne. mesh_Nr)) then
            ! Existence of (I+1/2, J+1/2)
                cell_I = index_I
                cell_J = index_J
                mat_pos_rel_indx=(/5, 6, 8, 9/) !just reuse array
                coeff1 = -coeff1_ab 
                coeff2 = -coeff2_ab 
                call col_3_LU_matrix_ftn(op_mode, cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
            endif

            if((index_I .ne. 1) .and. (index_J .ne. mesh_Nr)) then
            ! Existence of (I-1/2, J+1/2)
                cell_I = index_I-1
                cell_J = index_J
                mat_pos_rel_indx=(/2, 3, 5, 6/)
                coeff1 = -coeff1_ab
                coeff2 =  coeff2_ab 
                call col_3_LU_matrix_ftn(op_mode, cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
            endif
            
            if((index_I .ne. mesh_Nz) .and. (index_J .ne. 1)) then
            ! Existence of (I+1/2, J-1/2)
                cell_I = index_I
                cell_J = index_J-1
                mat_pos_rel_indx=(/4, 5, 7, 8/)
                coeff1 =  coeff1_ab
                coeff2 = -coeff2_ab 
                call col_3_LU_matrix_ftn(op_mode, cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
            endif

            if( (index_I .ne. 1) .and. (index_J .ne. 1) ) then
            ! Existence of (I-1/2, J-1/2)
                cell_I = index_I-1
                cell_J = index_J-1
                mat_pos_rel_indx=(/1, 2, 4, 5/)
                coeff1 = coeff1_ab
                coeff2 = coeff2_ab
                call col_3_LU_matrix_ftn(op_mode, cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
            endif
              
           enddo !index_J
   enddo !index_I


end subroutine col_3_LU_matrix

subroutine col_3_E_and_D_s(cs1, cs2, f_half, dfdr, dfdz, Ms, EDs)
! cs1 : colliding, cs2 : target
! f_half, dfdr, dfdz => cs2
    use sml_module, only : sml_mype
    use col_module
    implicit none
    type(col_3_core_type) :: cs1, cs2
    real (8), dimension(col_3_nvr-1,col_3_nvz-1) :: f_half, dfdr, dfdz 
    real (8), dimension(4,(col_3_nvr-1)*(col_3_nvz-1), (col_3_nvr-1)*(col_3_nvz-1)) :: Ms
    real (8), dimension(6,col_3_nvr-1, col_3_nvz-1) :: EDs
    integer :: index_I, index_J, index_ip, index_jp, mesh_Nrm1, mesh_Nzm1, index_2dp, index_2D, tmpi
    real (8) :: tmpr(6), tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v

    mesh_Nrm1 = col_3_nvr-1
    mesh_Nzm1 = col_3_nvz-1

    !M_zr => Ms(2,:,:) 
    !M_rc => Ms(2,:,:)
    !M_zc => Ms(3,:,:)

    do index_I=1, mesh_Nzm1
        do index_J=1, mesh_Nrm1

            index_2D = index_J+mesh_Nrm1*(index_I-1)
            tmpr=0D0

            do index_ip = 1, mesh_Nzm1
                do index_jp = 1, mesh_Nrm1

                    index_2dp = index_jp+mesh_Nrm1*(index_ip-1)
                    tmp_vol = cs2%local_center_volume(index_jp)
                    tmp_f_half_v = f_half(index_jp, index_ip) * tmp_vol
                    tmp_dfdr_v = dfdr(index_jp, index_ip) * tmp_vol
                    tmp_dfdz_v = dfdz(index_jp, index_ip) * tmp_vol

                    ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
                    !  Ms = 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za 
                    tmpr(1:3) = tmpr(1:3) + Ms(1:3,index_2dp,index_2D)*tmp_f_half_v
                    tmpr(4)   = tmpr(4)   + Ms(2,index_2dp,index_2D)*tmp_f_half_v
                    tmpr(5)   = tmpr(5)   + Ms(4,index_2dp,index_2D)*tmp_dfdr_v + Ms(2,index_2dp,index_2D)*tmp_dfdz_v
                    tmpr(6)   = tmpr(6)   + Ms(2,index_2D,index_2dp)*tmp_dfdr_v + Ms(3,index_2dp,index_2D)*tmp_dfdz_v  !This works for same species

                    do tmpi=1,6
                       if(tmpr(tmpi) .eq. tmpr(tmpi)*1.1D0 .and. tmpr(tmpi) .eq. tmpr(tmpi)+0.1D0) then
                           print *, 'Nan appears in calculating ED_values', sml_mype, index_I, index_J, index_ip, index_jp, Ms(:,index_2dp, index_2D), tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v, tmp_vol
                           stop
                       endif
                    enddo
                enddo !index_jp
            enddo !index_ip

            ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
            !mass correction
            tmpr(1:4) = tmpr(1:4)/cs1%mass
            tmpr(5:6) = tmpr(5:6)/cs2%mass
            EDs(:,index_J,index_I) = tmpr
        enddo
    enddo

end subroutine col_3_E_and_D_s


subroutine col_3_E_and_D_m(cs1, cs2, f_half, dfdr, dfdz, Ms, EDs)
! cs1 : colliding, cs2 : target
! f_half, dfdr, dfdz => cs2
    use sml_module, only : sml_mype
    use col_module
    implicit none
    type(col_3_core_type) :: cs1, cs2
    real (8), dimension(col_3_nvr-1,col_3_nvz-1) :: f_half, dfdr, dfdz 
    real (8), dimension(5,(col_3_nvr-1)*(col_3_nvz-1), (col_3_nvr-1)*(col_3_nvz-1)) :: Ms
    real (8), dimension(6,col_3_nvr-1, col_3_nvz-1) :: EDs
    integer :: index_I, index_J, index_ip, index_jp, mesh_Nrm1, mesh_Nzm1, index_2dp, index_2D, tmpi
    real (8) :: tmpr(6), tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v

    mesh_Nrm1 = col_3_nvr-1
    mesh_Nzm1 = col_3_nvz-1

    !M_zr => Ms(2,:,:) 
    !M_rc => Ms(2,:,:)
    !M_zc => Ms(3,:,:)

    do index_I=1, mesh_Nzm1
        do index_J=1, mesh_Nrm1

            index_2D = index_J+mesh_Nrm1*(index_I-1)
            tmpr=0D0

            do index_ip = 1, mesh_Nzm1
                do index_jp = 1, mesh_Nrm1

                    index_2dp = index_jp+mesh_Nrm1*(index_ip-1)
                    tmp_vol = cs2%local_center_volume(index_jp)
                    tmp_f_half_v = f_half(index_jp, index_ip) * tmp_vol
                    tmp_dfdr_v = dfdr(index_jp, index_ip) * tmp_vol
                    tmp_dfdz_v = dfdz(index_jp, index_ip) * tmp_vol

                    ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
                    !  Ms = 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za 
                    tmpr(1:3) = tmpr(1:3) + Ms(1:3,index_2dp,index_2D)*tmp_f_half_v
                    tmpr(4)   = tmpr(4)   + Ms(2,index_2dp,index_2D)*tmp_f_half_v
                    tmpr(5)   = tmpr(5)   + Ms(4,index_2dp,index_2D)*tmp_dfdr_v + Ms(2,index_2dp,index_2D)*tmp_dfdz_v
                    tmpr(6)   = tmpr(6)   + Ms(5,index_2dp,index_2D)*tmp_dfdr_v + Ms(3,index_2dp,index_2D)*tmp_dfdz_v  !Multi-species 

                    do tmpi=1,6
                       if(tmpr(tmpi) .eq. tmpr(tmpi)*1.1D0 .and. tmpr(tmpi) .eq. tmpr(tmpi)+0.1D0) then
                           print *, 'Nan appears in calculating ED_values', sml_mype, index_I, index_J, index_ip, index_jp, Ms(:,index_2dp, index_2D), tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v, tmp_vol
                           stop
                       endif
                    enddo
                enddo !index_jp
            enddo !index_ip

            ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez 
            !mass correction
            tmpr(1:4) = tmpr(1:4)/cs1%mass
            tmpr(5:6) = tmpr(5:6)/cs2%mass
            EDs(:,index_J,index_I) = tmpr
        enddo
    enddo

end subroutine col_3_E_and_D_m


subroutine col_3_f_df(delta_r, delta_z, mesh_dr, mesh_dz, f, f_half, dfdr, dfdz)
   use col_module
   implicit none
   real(kind=8) :: delta_z(col_3_nvz-1), delta_r(col_3_nvr-1)
   real(kind=8) :: mesh_dr, mesh_dz
   real(kind=8), dimension(col_3_nvr, col_3_nvz) :: f
   real(kind=8), dimension(col_3_nvr-1, col_3_nvz-1) :: f_half, dfdr, dfdz  
   real(kind=8) :: delz, cdelz, delr, cdelr, tmpr1, tmpr2, tmpr3, tmpr4
   integer :: index_I, index_J, mesh_Nzm1, mesh_Nrm1

   mesh_Nrm1 = col_3_nvr-1
   mesh_Nzm1 = col_3_nvz-1

   do index_I=1,mesh_Nzm1
       delz = delta_z(index_I)
       cdelz = 1D0-delz
       do index_J=1, mesh_Nrm1
           delr = delta_r(index_J)
           cdelr = 1D0-delr

           tmpr1 = f(index_J, index_I)
           tmpr2 = f(index_J+1, index_I)
           tmpr3 = f(index_J,index_I+1)
           tmpr4 = f(index_J+1,index_I+1)

           f_half(index_J, index_I) = tmpr1 * delr*delz &
                                    + tmpr3 * delr*cdelz &
                                    + tmpr2 * cdelr*delz &
                                    + tmpr4 * cdelr*cdelz
           dfdr(index_J, index_I) = (tmpr2 - tmpr1)*delz &
                                    +(tmpr4 - tmpr3)*cdelz
           dfdz(index_J, index_I) = (tmpr3 - tmpr1)*delr &
                                    +(tmpr4 - tmpr2)*cdelr
       enddo
   enddo
   dfdr = dfdr/mesh_dr
   dfdz = dfdz/mesh_dz
end subroutine col_3_f_df

subroutine col_3_angle_avg_s(cs, Ms)
    use sml_module, only : sml_pi, sml_mype
    use col_module
    implicit none
    type(col_3_core_type) :: cs
    real (8), dimension(4, (col_3_nvr-1)*(col_3_nvz-1), (col_3_nvr-1)*(col_3_nvz-1)) :: Ms 
    integer :: index_I, index_ip, index_J, index_jp
    integer :: index_rz, index_rc, index_ac, index_az
    real (8) :: r, z, a, c, dz, dz2, r2, a2, lambda, k, k_eff, kp1_sqrt, km1, EK, EE
    integer :: vpic_ierr
    real (8) :: I1, I2, I3, temp_cons, temp_cons_dz, temp_cons_I1, temp_cons_I2, temp_cons_I3, temp_cons_dz_I2, temp_cons_dz_I3
    real (8) :: tmp_a(6)
    integer :: mesh_Nrm1, mesh_Nzm1

    mesh_Nrm1 = col_3_nvr-1
    mesh_Nzm1 = col_3_nvz-1

    ! below 'do' loops from preconditioner 'DOUBLE_SYMMETRY' 
    !FOR OPTIMIZATION ON THIS ROUTINE, WE NEED TO SEPERATE elliptic calculation.
    do index_I=1,mesh_Nzm1
        z = cs%mesh_z_half(index_I)

        index_ip=1
        do while(index_ip .le. index_I)
            c = cs%mesh_z_half(index_ip)
            dz = z-c
            dz2 = dz*dz

            do index_J=1,mesh_Nrm1
                r=cs%mesh_r_half(index_J)
                r2 = r*r
                index_rz = index_J +mesh_Nrm1*(index_I-1)
                index_rc = index_J +mesh_Nrm1*(index_ip-1)

                index_jp=1
                do while(index_jp .le. index_J)
                    a=cs%mesh_r_half(index_jp)
                    a2 = a*a
                    index_ac = index_jp+mesh_Nrm1*(index_ip-1)
                    index_az = index_jp+mesh_Nrm1*(index_I-1)


                    lambda = r2+a2+dz2     !symmetric (r<->a, z<->c) 
                                           ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z)
                    k = 2D0*r*a/lambda     !symmetric
                    k_eff = 2D0*k/(1D0+k)
                    kp1_sqrt = sqrt(1D0+k)
                    km1 = k-1D0

                    call elliptics(k_eff,EK,EE,vpic_ierr)

                    !Calulation of M coeff. (all symmetric)
                    I1 = -4D0*((1D0+k)*EE-EK)/(k*k*kp1_sqrt)
                    I2 = -2D0*EE/(km1*kp1_sqrt)
                    I3 = -2D0*(EE+km1*EK)/(km1*k*kp1_sqrt )
                    !I4 = I2-I1 !For exact Numerical Conservation, and It shoudl be mathematically
                    temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda))
                    temp_cons_dz = temp_cons*dz
                    temp_cons_I1 = temp_cons*I1
                    temp_cons_I2 = temp_cons*I2
                    temp_cons_I3 = temp_cons*I3
                    temp_cons_dz_I2 = temp_cons_dz*I2
                    temp_cons_dz_I3 = temp_cons_dz*I3


                    !NOTE THAT CONSIDERING SYMMETRY, we already have 
                    ! I1, I2, I3 at (r,z,a,c), (a,z,r,c), (r,c,a,z), and (a,c,r,z) ->(jp,ip,J,I)
                    ! Using this values, we calculate symmetric and asymmetric compoenents of matrix
                    ! index_rz = index_J +mesh_Nrm1*index_I
                    ! index_rc = index_J +mesh_Nrm1*index_ip 
                    ! index_ac = index_jp+mesh_Nrm1*index_ip
                    ! index_az = index_jp+mesh_Nrm1*index_I
                    tmp_a(1) = temp_cons_I1*a2+temp_cons_I2*dz2    
                    tmp_a(2) = temp_cons_I1*r2+temp_cons_I2*dz2
                    tmp_a(3) = temp_cons_dz_I3*a-temp_cons_dz_I2*r
                    tmp_a(4) = temp_cons_dz_I2*a-temp_cons_dz_I3*r
                    tmp_a(5) = temp_cons_I2*(r2+a2) - temp_cons_I3*2D0*r*a 
                    tmp_a(6) = temp_cons_I1*r*a + temp_cons_I3*dz2  

                    !(a,c,r,z)=(r,z,a,c)
                    ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra
                    if(index_ac .eq. index_rz) then
                        Ms(1, index_ac, index_ac) = 0D0
                        Ms(2, index_ac, index_ac) = 0D0
                    else
                        Ms(1, index_ac, index_rz) = tmp_a(1)   
                        Ms(1, index_rz, index_ac) = tmp_a(2)
                        Ms(2, index_ac, index_rz) = tmp_a(3)
                        Ms(2, index_rz, index_ac) = tmp_a(4)

                        Ms(3, index_ac, index_rz) = tmp_a(5)
                        Ms(3, index_rz, index_ac) = tmp_a(5) ! Using Transpose property 
                        Ms(4, index_ac, index_rz) = tmp_a(6)
                        Ms(4, index_rz, index_ac) = tmp_a(6) ! Using Transpose property
                    endif

                    !(a,z,r,c) = (r,c,a,z)
                    if(index_az .eq. index_rc) then
                        Ms(1, index_az, index_az) = 0D0
                        Ms(2, index_az, index_az) = 0D0
                    else
                        Ms(1, index_az, index_rc) = tmp_a(1)
                        Ms(1, index_rc, index_az) = tmp_a(2)
                        Ms(2, index_az, index_rc) = -tmp_a(3)
                        Ms(2, index_rc, index_az) = -tmp_a(4)

                        Ms(3, index_rc, index_az) = tmp_a(5)
                        Ms(3, index_az, index_rc) = tmp_a(5) ! Using Transpose property 
                        Ms(4, index_rc, index_az) = tmp_a(6)
                        Ms(4, index_az, index_rc) = tmp_a(6) ! Using Transpose property 
                    endif

                    index_jp=index_jp+1
                enddo
            enddo

            index_ip = index_ip+1
        enddo
    enddo
end subroutine col_3_angle_avg_s

subroutine col_3_core_init(spi, ci, dist_ni, mesh_dri, mesh_dzi, dlxi, voli, &
                             spe, ce, dist_ne, mesh_dre, mesh_dze, dlxe, vole )
    use sml_module, only : sml_e_charge
    use ptl_module
    use col_module
    implicit none
    type(species_type) :: spi, spe
    type(col_3_core_type) :: ci, ce
    real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_ni, dist_ne   ! local 
    real (8) :: mesh_dri, mesh_dzi, mesh_dre, mesh_dze, dlxi, dlxe
    real (kind=8), dimension(1:col_3_nvr) :: voli, vole      ! local
    real (8) :: densi, dense, numeric_Ti, numeric_Te, eni, ene
    integer :: index_i, mesh_Nz

    mesh_Nz = col_3_nvz

    call col_3_core_sm_init_s(ci, spi%mass, mesh_dri, mesh_dzi, dlxi, voli)
    call col_3_core_sm_init_s(ce, spe%mass, mesh_dre, mesh_dze, dlxe, vole)

    ! Get Eq. Temperature
    densi = 0D0;
    dense = 0D0;
    eni = 0D0;
    ene = 0D0;
    do index_i=1, mesh_Nz
        densi = densi + sum(dist_ni(:,index_i)*voli)
        dense = dense + sum(dist_ne(:,index_i)*vole)
        eni = eni+sum(dist_ni(:,index_i)*(ci%mesh_r**2+(dlxi+mesh_dzi*(index_i-1))**2)*voli) !n_i*v_i^2 
        ene = ene+sum(dist_ne(:,index_i)*(ce%mesh_r**2+(dlxe+mesh_dze*(index_i-1))**2)*vole) !n_e*v_e^2
    end do
    eni = spi%mass*eni
    ene = spe%mass*ene
    numeric_Ti = eni/(3D0*densi*sml_e_charge)
    numeric_Te = ene/(3D0*dense*sml_e_charge)

    call col_3_core_s_init_e(ci, spi%mass, mesh_dri, mesh_dzi, numeric_Ti, densi, eni)
    call col_3_core_s_init_e(ce, spe%mass, mesh_dre, mesh_dze, numeric_Te, dense, ene)

end subroutine col_3_core_init

subroutine col_3_core_s_init(sp, cs, dist_n, mesh_dr, mesh_dz, dlx, vol)
    use sml_module, only : sml_mype, sml_e_charge
    use ptl_module
    use col_module
    implicit none
    type(species_type) :: sp
    type(col_3_core_type) :: cs
    real (kind=8), dimension(col_3_nvr, col_3_nvz) :: dist_n   ! local 
    real (8) :: mesh_dr, mesh_dz, dlx
    real (kind=8), dimension(1:col_3_nvr) :: vol      ! local
    real (8) :: dens, ens, numeric_T
    integer :: index_i, mesh_Nz

    mesh_Nz = col_3_nvz

    call col_3_core_sm_init_s(cs, sp%mass, mesh_dr, mesh_dz, dlx, vol)

    ! Get Eq. Temperature
    dens = 0D0;
    ens = 0D0;
    do index_i=1, mesh_Nz
        dens = dens + sum(dist_n(:,index_i)*vol)
        ens = ens+sum(dist_n(:,index_i)*(cs%mesh_r**2+(dlx+mesh_dz*(index_i-1))**2)*vol) !n_i*v_i^2 
    end do
    ens = ens*sp%mass
    numeric_T = ens/(3D0*dens*sml_e_charge)

    call col_3_core_s_init_e(cs, sp%mass, mesh_dr, mesh_dz, numeric_T, dens, ens)

end subroutine col_3_core_s_init

subroutine col_3_core_sm_init_s( cs, mass, mesh_dr, mesh_dz, dlx, vol )
    use sml_module, only : sml_e_charge
    use col_module
    implicit none
    type(col_3_core_type) :: cs
    real (kind=8) :: mass, mesh_dr, mesh_dz,  dlx
    real (kind=8), dimension(1:col_3_nvr) :: vol      ! local
    integer :: index_i, index_e, index_s
    integer :: mesh_Nr, mesh_Nz, mesh_Nrm1, mesh_Nzm1, mesh_Nrzm1
    real (kind=8) :: ddum0

    ! To be deallocated at `col_3_core_s_deallocation' or `col_3_core_deallocation' subroutine under col_3_core_m subroutine
    allocate(cs%mesh_r(col_3_nvr), cs%mesh_z(col_3_nvz), cs%mesh_r_half(col_3_nvr-1), cs%mesh_z_half(col_3_nvz-1), &
             cs%local_center_volume(col_3_nvr-1), cs%vol(col_3_nvr))

    mesh_Nr = col_3_nvr
    mesh_Nz = col_3_nvz

    mesh_Nrm1 = mesh_Nr-1
    mesh_Nzm1 = mesh_Nz-1
    mesh_Nrzm1 = mesh_Nrm1*mesh_Nzm1

    ! mass, mesh_dr, mesh_dz
    cs%mass        = mass
    cs%mesh_dr     = mesh_dr
    cs%mesh_dz     = mesh_dz

    ! mesh_r, mesh_r_half, mesh_z, mesh_z_half
    do index_i=1, col_3_nvr-1
        cs%mesh_r(index_i) = mesh_dr*(index_i-1)
        cs%mesh_r_half(index_i) = mesh_dr*(index_i-0.5D0)
    end do
    cs%mesh_r(mesh_Nr) = mesh_dr*(mesh_Nrm1)
    do index_i=1, col_3_nvz-1
        cs%mesh_z(index_i) = dlx + mesh_dz*(index_i-1)
        cs%mesh_z_half(index_i) = dlx + mesh_dz*(index_i-0.5D0)
    end do
    cs%mesh_z(mesh_Nz) = dlx + mesh_dz*(mesh_Nzm1)
    
    ! local_center_volume 
    cs%local_center_volume = cs%mesh_r_half*mesh_dr*mesh_dz  !volume centered at a cell

    ! volume
    cs%vol = vol !volume evaluated at NODEs

end subroutine col_3_core_sm_init_s

subroutine col_3_core_s_init_e( cs, mass, mesh_dr, mesh_dz, numeric_Ts, dens , ens) 
    use sml_module, only : sml_e_charge, sml_mype, sml_electron_on
    use col_module
    implicit none
    type(col_3_core_type) :: cs
    real (kind=8) :: mass, mesh_dr, mesh_dz, numeric_Ts, dens, ens
    integer :: i, index_e, index_s
    integer :: mesh_Nr, mesh_Nz, mesh_Nrm1, mesh_Nzm1
    real (kind=8) :: factor_1r(col_3_nvr-1), factor_2r(col_3_nvr-1), factor_3r(col_3_nvr-1), factor_Qr(col_3_nvr-1) 
    real (kind=8) :: factor_1z(col_3_nvz-1), factor_2z(col_3_nvz-1), factor_3z(col_3_nvz-1), factor_Qz(col_3_nvz-1) 

    if(sml_electron_on) then
        allocate( cs%delta_r(2, col_3_nvr-1), cs%delta_z(2, col_3_nvz-1) )   ! 2nd dimension : 1->same species, 2->different species
    else
        allocate( cs%delta_r(1, col_3_nvr-1), cs%delta_z(1, col_3_nvz-1) )   ! 2nd dimension : 1->same species, 2->different species
    endif

    mesh_Nr = col_3_nvr
    mesh_Nz = col_3_nvz

    mesh_Nrm1 = mesh_Nr-1
    mesh_Nzm1 = mesh_Nz-1

    ! numeric_vth2, delta_r
    cs%numeric_vth2 = numeric_Ts*sml_e_charge/mass


    if( .not. sml_electron_on) then
        !Below will be used since ver 0.1 - originally from preconditioner (.not. ORG_COEFF)
        factor_1r = cs%numeric_vth2/(mesh_dr*cs%mesh_r_half)
        factor_2r = 1D0/(exp(1D0/factor_1r)-1D0)
        !factor_3r = 2D0*exp(1D0/factor_1r)-factor_1r  : This has an issue -20130228
        factor_3r = 2D0*factor_2r-factor_1r
        factor_Qr = sqrt(factor_3r*factor_3r-4D0*(factor_2r-factor_1r)*factor_2r)

        cs%delta_r(1,:) = (-factor_3r+ factor_Qr)*0.5D0
        
        ! delta_z
        factor_1z = cs%numeric_vth2/(mesh_dz*cs%mesh_z_half)
        factor_2z = 1D0/(exp(1D0/factor_1z)-1D0)
        !factor_3z = 2D0*exp(1D0/factor_1z)-factor_1z : This has an issue -20130228
        factor_3z = 2D0*factor_2z-factor_1z
        factor_Qz = sqrt(factor_3z*factor_3z-4D0*(factor_2z-factor_1z)*factor_2z)

        if(mod(mesh_Nz,2) .eq. 0) then
            index_e = mesh_Nz/2-1;
            index_s = mesh_Nz/2+1;

            cs%delta_z(1,1:index_e) = (-factor_3z(1:index_e)- factor_Qz(1:index_e))*0.5D0
            cs%delta_z(1,index_s:mesh_Nzm1) = (-factor_3z(index_s:mesh_Nzm1)+factor_Qz(index_s:mesh_Nzm1))*0.5D0
            cs%delta_z(1,mesh_Nz/2) = 0.5D0
        else
            index_e= (mesh_Nz-1)/2
            index_s = (mesh_Nz+1)/2

            cs%delta_z(1,1:index_e) = (-factor_3z(1:index_e) - factor_Qz(1:index_e))*0.5D0
            cs%delta_z(1,index_s:mesh_Nzm1) = (-factor_3z(index_s:mesh_Nzm1)+factor_Qz(index_s:mesh_Nzm1))*0.5D0
        endif

        !do i=1, mesh_Nzm1
        !    if( cs%delta_z(i) .eq. cs%delta_z(i)*1.1D0 .and. cs%delta_z(i) .eq. cs%delta_z(i) + 1D0 ) then
        !        print *, 'Nan for delta_z', factor_1z(i), factor_2z(i), factor_3z(i), factor_Qz(i), cs%delta_z(i), numeric_Teq, mesh_dz, cs%mesh_z_half(i)
        !        stop
        !    endif
        !enddo
    else
       !for multi-spcies
       ! 2013-02-27
       ! This routine should be updated with extended analytic calculation.
       ! Contrast to same species collisions, each species in different species collisions has
       ! their own temperature and momentum. Therefore, if we want to evaluate coefficients of delta for equilibrium state,
       ! we need to calculate final temperature and momentum, i.e. shifted Maxwellian.
       ! Since I don't have enough time to do this, I just postpone this extention, and just use 1/2.
       ! Based on my experiences in past, constant 1/2 doest not lead to bad results.
       cs%delta_r = 0.5D0
       cs%delta_z = 0.5D0
    endif

    cs%numeric_Teq = numeric_Ts
    cs%numeric_T   = numeric_Ts
    cs%dens        = dens
    cs%ens         = ens
 
end subroutine col_3_core_s_init_e

subroutine elliptics(k,ellip_K,ellip_E,vpic_ierr)
    use sml_module, only : sml_pi
    use col_module
    implicit none
    real(kind=8) :: k, ellip_K, ellip_E
    integer :: vpic_ierr, k_index
    real(kind=8) :: ell_B, ell_D, BX_star, B0_star, DX_star, D0_star, X, mc, kc
    integer :: iter, Bnp1, Dnp1
    real (kind=8), dimension(12) :: &
       B0to1 = (/ 0.790401413584395132,  0.102006266220019155, 0.039878395558551461, &
             0.021737136375982167,  0.013960979767622058, 0.009892518822669142, &
             0.007484612400663336,  0.005934625664295474, 0.004874249053581664, &
             0.004114606930310886,  0.003550452989196177, 0.003119229959988475 /), &
       D0to1 = (/ 0.800602040206397048,  0.313994477771767757, 0.205913118705551955, &
             0.157744346538923994,  0.130595077319933092, 0.113308474489758567, &
             0.101454199173630195,  0.092918784207297437, 0.086565380148168087, &
             0.081727984665103014,  0.077990665729107038, 0.075080426851268007 /), &
       B1to2 = (/ 0.801024064452844894,  0.110695344529634015, 0.047348746716993718, &
             0.028484367255041423,  0.020277811444003597, 0.015965005853099119, &
             0.013441320273553635,  0.011871565736951440, 0.010868363672485521, &
             0.010231587232710565,  0.009849585546666211, 0.009656606347153765 /), &
       D1to2 = (/ 0.834232667811735098,  0.360495281619098276, 0.262379664114505869, &
             0.223723944518094276,  0.206447811775681053, 0.199809440876486856, &
             0.199667451603795275,  0.204157558868236842, 0.212387467960572375, &
             0.223948914061499360,  0.238708097425597860, 0.256707203545463756 /)
    real (kind=8), dimension(13) :: &
       B2to3 = (/ 0.812597772919920493,  0.121109617945510113, 0.057293376831239877, &
             0.038509451602167328,  0.030783430301775233, 0.027290564934732527, &
             0.025916369289445199,  0.025847203343361799, 0.026740923539348855, &
             0.028464314554825705,  0.030995446237278954, 0.034384369179940976, &
             0.038738002072493936 /), &
       D2to3 = (/ 0.873152581892675550,  0.420622230667770216, 0.344231061559450379, &
             0.331133021818721762,  0.345277285052808412, 0.377945322150393392, &
             0.427378012464553881,  0.494671744307822406, 0.582685115665646201, &
             0.695799207728083165,  0.840018401472533403, 1.023268503573606061, &
             1.255859085136282496 /), &
       B3to4 = (/ 0.825323557983515895, 0.133862116083687790, 0.071011293597988675, &
                  0.054178477417387376, 0.049451744948102993, 0.050222196224107476, &
                  0.054742913171830353, 0.062746257927001699, 0.074669881043476886, &
                  0.091480845177733472, 0.114705092110997824, 0.146571132581439876, &
                  0.190257137333846284 /)
    real (kind=8), dimension(14) :: &
       D3to4 = (/ 0.919027039242097348, 0.501002159288247514, 0.468831270566456863, &
                  0.517714227776400015, 0.620843391317303107, 0.782364393786869723, &
                  1.019114535076102913, 1.359345202748496052, 1.845717302358827942, &
                  2.541071703153920729, 3.537404655208041337, 4.969296002977425930, &
                  7.033822870030031126, 10.02004322503447140 /)
    real (kind=8), dimension(13) :: &
       B4to5 = (/ 0.839479570270612971, 0.149916440306396336, 0.090831935819428835, &
                  0.080347033483341786, 0.085638440500470454, 0.101954725932990372, &
                  0.130574811533616015, 0.176105076358849928, 0.246835164402955447, &
                  0.356424476867718855, 0.527002562230102743, 0.794389634259304750, &
                  1.216762532429718021 /)
    real (kind=8), dimension(16) :: &
       D4to5 = (/ 0.974404366546369673, 0.613246805394160910, 0.671096669502166996, &
                  0.870727620185086140, 1.229542231202690761, 1.826605967544420569, &
                  2.806934530997762740, 4.418789329084028134, 7.083236057478765325, &
                  11.51508812055758294, 18.93151118599927464, 31.41199693820496388, &
                  52.52072945457582854, 88.38485473506529806, 149.5663744939804784, &
                  254.3179084310411743 /)
    real (kind=8), dimension(14) :: &
       B5to6 = (/ 0.855469615156419991, 0.170896072689739584, 0.121335229026948226, &
                  0.128201883574947410, 0.164687281451527560, 0.237418908749381742, &
                  0.369208104716495452, 0.605658733847927717, 1.033705561557812744, &
                  1.818988489363267885, 3.279377651273850938, 6.029888380717536331, &
                  11.26979685557794172, 21.35457785038283450 /)
    real (kind=8), dimension(17) :: &
       D5to6 = (/ 1.043455295115133534, 0.779625721928504850, 1.029742360932067582, &
                  1.622037223411353130, 2.787989531185347620, 5.048381487372069147, &
                  9.463277611943484295, 18.18148994942766790, 35.58098059117916870, &
                  70.63393546191445013, 141.8285800834330593, 287.4487512501321663, &
                  587.1153846499230762, 1207.065435225480616, 2495.588727248664223, &
                  5184.692429394806441, 10817.21333690413275 /)
    real (kind=8), dimension(16) :: &
       B6to7 = (/ 0.873920061848643136, 0.199814057482376946, 0.172769615878015213, &
                  0.228106913284202167, 0.370468141118071220, 0.679271252884820555, &
                  1.348008496681757302, 2.827670976853820704, 6.179468250123914084, &
                  13.93568601034281150, 32.21892928105972203, 76.00696295922610103, &
                  182.3214490877540696, 443.5150764411264816, 1091.854722902838829, &
                  2715.765866403819588 /)
    real (kind=8), dimension(18) :: &
       D6to7 = (/ 1.133678336575733166, 1.048643173729970391, 1.753465041198464516, &
                  3.523182726803385513, 7.749476413813974582, 17.98645005585073306, &
                  43.25591634623261333, 106.6815344540960170, 268.0984865731174340, &
                  683.6241148502898048, 1763.497085219187407, 4592.374753831163809, &
                  12053.44101904888928, 31846.66302074208170, 84621.22135905680802, &
                  225956.4231829078900, 605941.5172817588600, 1631082.599539268321 /)
    real (kind=8), dimension(19) :: &
       B7to8 = (/ 0.895902820924731621, 0.243140003766786662, 0.273081875594105532, &
                  0.486280007533573324, 1.082747437228230918, 2.743445290986452500, &
                  7.555817828670234627, 22.05194082493752427, 67.15640644740229408, &
                  211.2722537881770962, 681.9037843053270682, 2246.956231592536517, &
                  7531.483865999711792, 25608.51260130241579, 88140.74740089604971, &
                  306564.4242098446591, 1076036.077811072194, 3807218.502573632648, &
                  13566382.24422139551 /)
    real (kind=8), dimension(21) :: &
       D7to8 = (/ 1.260612826574911614, 1.548665638082676581, 3.553669411871607615, &
                  9.900444676104398756, 30.32056661745247199, 98.18025865888308915, &
                  329.7710104345570550, 1136.655989742890393, 3993.834335746229798, &
                  14242.72958655527085, 51394.75729168872096, 187246.7029146231521, &
                  687653.0923753899027, 2542385.535653982270, 9453781.219347490272, &
                  35328363.01797091708, 132593262.3833930149, 499544968.1840548215, &
                  1888409347.294438724, 7160267534.478937192, 27223307946.96339622 /)
    real (kind=8), dimension(15) :: &
       B8to85 = (/ 0.915922052601931494, 0.294714252429483394, 0.435776709264636140, &
                   1.067328246493644239, 3.327844118563268085, 11.90406004445092906, &
                   46.47838820224626394, 192.7556002578809477, 835.3356299261900064, &
                   3743.124548343029103, 17219.07731004063094, 80904.60401669850158, &
                   386808.3292751742460, 1876487.670110449342, 9216559.908641567755 /)
    real (kind=8), dimension(18) :: &
       D8to85 = (/ 1.402200569110579095, 2.322205897861749447, 7.462158366466719683, &
                   29.43506890797307903, 128.1590924337895775, 591.0807036911982326, &
                   2830.546229607726377, 13917.76431889392230, 69786.10525163921228, &
                   355234.1420341879635, 1830019.186413931054, 9519610.812032515607, &
                   49920868.75574849454, 263567700.9826023474, 1399645765.120061119, &
                   7469935792.837635005, 40041555958.35610574, 215463066814.4966654 /)
    real (kind=8), dimension(19) :: &
       B85to9 = (/ 0.931906061029524828, 0.348448029538453861, 0.666809178846938248, &
                   2.210769135708128663, 9.491765048913406881, 47.09304791027740853, &
                   255.9200460211233087, 1480.029532675805408, 8954.040904734313578, &
                   56052.48220982686950, 360395.7241626000917, 2367539.415273216078, &
                   15829949.57277684102, 107415809.3278511100, 738058546.0239595692, &
                   5126022002.555101497, 35935340655.02416589, 253988125761.2812212, &
                   1808180007145.359570 /)
    real (kind=8), dimension(21) :: &
       D85to9 = (/ 1.541690112721819084, 3.379176214579645449, 14.94058385670236672, &
                   81.91773929235074881, 497.4900546551479866, 3205.184010234846235, &
                   21457.32237355321926, 147557.0156564174712, 1035045.290185256525, &
                   7371922.334832212125, 53143443.95142401142, 386882347.5795976313, &
                   2839458401.528033778, 20982661229.43898942, 155961775401.7662418, &
                   1165096220419.884791, 8742012983013.913805, 65847254626723.66919, &
                   497679873706243.4393, 3773018634056605.405, 28682631948378196.60 /)
    real (kind=8) , dimension(14) :: &
       elliptic_BX = (/ 0D0, -0.25D0, -3.125D-2, -1.171875D-2, -6.103515625D-3, -3.7384033203125D-3, &
                      -2.5234222412109375D-3, -1.817464828491210938D-3, -1.371212303638458252D-3, &
                      -1.071259612217545509D-3, -8.599834109190851450D-4, -7.055772985040675849D-4, &
                      -5.893174027278291760D-4, -4.995976058381756957D-4 /), &
       elliptic_B0 = (/ 1D0, -0.25D0, 4.6875D-2, 2.34375D-2, 1.352945963541666667D-2, 8.740743001302083333D-3, &
                       6.0962677001953125D-3, 4.489111900329589844D-3, 3.441497071513107845D-3, &
                       2.721402519715151617D-3, 2.205478451662837336D-3, 1.823336289104075164D-3, &
                       1.532467036966436629D-3, 1.305981575356390531D-3 /)
    real (kind=8) , dimension(13) :: &
       elliptic_DX = (/ 0.5D0, -0.125D0, -2.34375D-2, -9.765625D-3, -5.340576171875D-3, -3.36456298828125D-3, &
                      -2.313137054443359375D-3, -1.687645912170410156D-3, -1.285511534661054611D-3, &
                      -1.011745189316570759D-3, -8.169842403731308877D-4, -6.735056031175190583D-4, &
                      -5.647625109475029603D-4 /), &
       elliptic_D0 = (/ -1D0, 0D0, 3.90625D-2, 2.018229166666666667D-2, 1.202901204427083333D-2, &
                       7.941436767578125000D-3, 5.623292922973632813D-3, 4.187006609780447824D-3, &
                       3.237116100665714060D-3, 2.576826204497751499D-3, 2.099504446134290894D-3, &
                       1.743372975543576159D-3, 1.470660484741195621D-3 /)
    integer, dimension(11) :: &
       elliptic_Bnum = (/ 11, 11, 12, 12, 12, 13, 15, 18, 14, 18, 13 /), &
       elliptic_Dnum = (/ 11, 11, 12, 13, 15, 16, 17, 20, 17, 20, 12 /)

    k_index = floor(k*10D0)
    select case(k_index)
      case(0)
          if( k .eq. 0) then
              ell_B = sml_pi*0.25
              ell_D = sml_pi*0.25
          else if((0D0 .lt. k) .and. (k .le. 0.1)) then
              Bnp1 = elliptic_Bnum(1)+1
              Dnp1 = elliptic_Dnum(1)+1
              ell_B = B0to1(Bnp1)
              ell_D = D0to1(Dnp1)
              kc = k-0.05D0
              do iter=1,elliptic_Bnum(1)
                  ell_B = ell_B*kc+B0to1(Bnp1-iter)
              enddo
              do iter=1,elliptic_Dnum(1)
                  ell_D = ell_D*kc+D0to1(Dnp1-iter)
              enddo
          endif
      case(1)
          Bnp1 = elliptic_Bnum(2)+1
          Dnp1 = elliptic_Dnum(2)+1
          ell_B = B1to2(Bnp1)
          ell_D = D1to2(Dnp1)
          kc = k-0.15D0
          do iter=1,elliptic_Bnum(2)
              ell_B = ell_B*kc+B1to2(Bnp1-iter)
          enddo
          do iter=1,elliptic_Dnum(2)
              ell_D = ell_D*kc+D1to2(Dnp1-iter)
          enddo
      case(2)
          Bnp1 = elliptic_Bnum(3)+1
          Dnp1 = elliptic_Dnum(3)+1
          ell_B = B2to3(Bnp1)
          ell_D = D2to3(Dnp1)
          kc = k-0.25D0
          do iter=1,elliptic_Bnum(3)
              ell_B = ell_B*kc+B2to3(Bnp1-iter)
          enddo
          do iter=1,elliptic_Dnum(3)
              ell_D = ell_D*kc+D2to3(Dnp1-iter)
          enddo
      case(3)
        Bnp1 = elliptic_Bnum(4)+1
        Dnp1 = elliptic_Dnum(4)+1
        ell_B = B3to4(Bnp1)
        ell_D = D3to4(Dnp1)
        kc = k-0.35D0
        do iter=1,elliptic_Bnum(4)
            ell_B = ell_B*kc+B3to4(Bnp1-iter)
        enddo
        do iter=1,elliptic_Dnum(4)
            ell_D = ell_D*kc+D3to4(Dnp1-iter)
        enddo
      case(4)
        Bnp1 = elliptic_Bnum(5)+1
        Dnp1 = elliptic_Dnum(5)+1
        ell_B = B4to5(Bnp1)
        ell_D = D4to5(Dnp1)
        kc = k-0.45D0
        do iter=1, elliptic_Bnum(5)
            ell_B = ell_B*kc+B4to5(Bnp1-iter)
        enddo
        do iter=1, elliptic_Dnum(5)
            ell_D = ell_D*kc+D4to5(Dnp1-iter)
        enddo
      case(5)
        Bnp1 = elliptic_Bnum(6)+1
        Dnp1 = elliptic_Dnum(6)+1
        ell_B = B5to6(Bnp1)
        ell_D = D5to6(Dnp1)
        kc = k-0.55D0
        do iter=1, elliptic_Bnum(6)
            ell_B = ell_B*kc+B5to6(Bnp1-iter)
        enddo
        do iter=1, elliptic_Dnum(6)
            ell_D = ell_D*kc+D5to6(Dnp1-iter)
        enddo
      case(6)
        Bnp1 = elliptic_Bnum(7)+1
        Dnp1 = elliptic_Dnum(7)+1
        ell_B = B6to7(Bnp1)
        ell_D = D6to7(Dnp1)
        kc = k-0.65D0
        do iter=1, elliptic_Bnum(7)
            ell_B = ell_B*kc+B6to7(Bnp1-iter)
        enddo
        do iter=1, elliptic_Dnum(7)
            ell_D = ell_D*kc+D6to7(Dnp1-iter)
        enddo
      case(7)
        Bnp1 = elliptic_Bnum(8)+1
        Dnp1 = elliptic_Dnum(8)+1
        ell_B = B7to8(Bnp1)
        ell_D = D7to8(Dnp1)
        kc = k-0.75D0
        do iter=1, elliptic_Bnum(8)
            ell_B = ell_B*kc+B7to8(Bnp1-iter)
        enddo
        do iter=1, elliptic_Dnum(8)
            ell_D = ell_D*kc+D7to8(Dnp1-iter)
        enddo
      case(8)
          if((0.8 .lt. k) .and. (k .le. 0.85)) then
              Bnp1 = elliptic_Bnum(9)+1
              Dnp1 = elliptic_Dnum(9)+1
              ell_B = B8to85(Bnp1)
              ell_D = D8to85(Dnp1)
              kc = k-0.825D0
              do iter=1, elliptic_Bnum(9)
                  ell_B = ell_B*kc+B8to85(Bnp1-iter)
              enddo
              do iter=1, elliptic_Dnum(9)
                  ell_D = ell_D*kc+D8to85(Dnp1-iter)
              enddo
          else if((0.85 .lt. k) .and. (k .le. 0.9)) then
              Bnp1 = elliptic_Bnum(10)+1
              Dnp1 = elliptic_Dnum(10)+1
              ell_B = B85to9(Bnp1)
              ell_D = D85to9(Dnp1)
              kc = k-0.875D0
              do iter=1, elliptic_Bnum(10)
                  ell_B = ell_B*kc+B85to9(Bnp1-iter)
              enddo
              do iter=1, elliptic_Dnum(10)
                  ell_D = ell_D*kc+D85to9(Dnp1-iter)
              enddo
          endif
      case(9)
        Bnp1 = elliptic_Bnum(11)+1
        Dnp1 = elliptic_Dnum(11)+1
        BX_star = elliptic_BX(Bnp1)
        B0_star = elliptic_B0(Bnp1)
        DX_star = elliptic_DX(Dnp1)
        D0_star = elliptic_D0(Dnp1)
        mc = 1D0-k
    
        do iter=1,elliptic_Bnum(11)
            BX_star = BX_star*mc + elliptic_BX(Bnp1-iter)
            B0_star = B0_star*mc + elliptic_B0(Bnp1-iter)
        enddo
        do iter=1,elliptic_Dnum(11)
            DX_star = DX_star*mc + elliptic_DX(Dnp1-iter)
            D0_star = D0_star*mc + elliptic_D0(Dnp1-iter)
        enddo
        X = -log(mc/16)
        ell_B = (B0_star+BX_star*X)/k
        ell_D = (D0_star+DX_star*X)/k
      case(10)
        ell_B = 0D0
        ell_D = 0D0
        vpic_ierr = 2
      case default
        print *, 'INPUT for elliptic integral is out of range :',k
        vpic_ierr = 2
    end select

    ellip_K = ell_B+ell_D
    ellip_E = ell_B+(1D0-k)*ell_D
end subroutine elliptics 

subroutine get_volume_col3(mode,grid)
  use grid_class
  use sml_module
  use ptl_module
  use eq_module
  use col_module
  use random_xgc

  implicit none
  integer :: mode      ! 0 : reference (or guide) grid, 1 : subgrid
  type(grid_type) :: grid
  integer ::  valid    !! number of valid particle 
  integer (8) ::  total    !! number of total particle generated
  integer :: l
  real (8) :: rdim, zdim  !! R,Z dimension of whole simulation region
  real (8) :: roffset, zoffset !! offset (start position) of simulation region
  real (8) :: mc_den
  real (8) :: r,z,psi,theta,dum1,dum2,rand
  !
  real (8), allocatable :: dum_c_3(:)      ! for NL VPIC
  !
  real (8) , external :: psi_interpol
  integer, external :: col_3_search_indx, col_3_search_indx_subg
  ! 2013-02-21 subgrid correction due to limiter
  integer :: grid_valid_num,i 

  call random_seed()
 
  ! dum_c_3 memory allocation 
  if(mode .eq. 0) then
      allocate(dum_c_3(0:col_3_ntotal_r-1))
  else if (mode .eq. 1) then
      allocate(dum_c_3(0:col_3_total_subg_num-1))
  endif

! simulation boundary is imposed 2001/01/24
  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  valid=0
  total=0
  
  ! initialize volume
  col_3_vol_r=1D-99
  
  do while(valid<sml_monte_num)
     !r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) ) 
     !z=zdim*ranx() + zoffset
     call random_number(rand)
     r=sqrt(roffset**2 + rdim*rand*(2D0*roffset + rdim) ) 
     call random_number(rand)
     z=zdim*rand + zoffset
     psi=psi_interpol(r,z,0,0)
     theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
     if(z < eq_axis_z) then
        theta= sml_2pi - theta    
     endif
     total=total+1
     
     if(mode .eq. 0) then
         l = col_3_search_indx(r,z,psi,theta,grid)
         if(l.lt.-1 .or. l .ge. col_3_ntotal_r) then
             print *, 'col_3_search_indx logic error', r, z, psi, theta, l
             stop
         endif
     else if (mode .eq. 1) then
         l = col_3_search_indx_subg(r,z,psi,theta,grid,0)
     endif
     
     if(l .ne. -1) then
         valid = valid+1 !marker density
         col_3_vol_r(l) = col_3_vol_r(l)+1D0
     endif
  enddo

  ! total sum
  dum1=real(total)
  
  call my_mpi_allreduce(dum1,dum2,1)
  if(mode .eq. 0) then
      call my_mpi_allreduce(col_3_vol_r,dum_c_3,col_3_ntotal_r)
  else if(mode .eq. 1) then
      call my_mpi_allreduce(col_3_vol_r,dum_c_3,col_3_total_subg_num)
  endif

  !subgrid correction
  if(mode .eq. 1) then
      grid_valid_num = 0
      do i=0, col_3_total_subg_num-1
          if(nint(col_3_vol_r(i)) .ne. 0) then
              col_3_grid_valid_list(grid_valid_num) = i
              grid_valid_num = grid_valid_num+1
          endif
      enddo
      col_3_total_subg_num = grid_valid_num
  endif
  

  mc_den=dum2/ (  rdim*zdim*sml_2pi*(roffset+0.5D0*rdim) )
  col_3_vol_r=dum_c_3/(sml_nphi_total*mc_den)

  deallocate(dum_c_3)
end subroutine get_volume_col3

subroutine col_3_setup(grid)
  use grid_class
  use eq_module
  use col_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  real (8) :: min_vol
  integer :: i, j, k
  integer :: sep
  ! 2013-02-13 for collision grid output
  real (kind=8) :: rdim, roffset, zdim, zoffset, spec_dr, spec_dz, spec_r, spec_z, spec_psi, spec_theta
  integer :: resol_r, resol_z, indx
  real (kind=8), external :: psi_interpol
  integer, external :: col_3_search_indx, col_3_search_indx_subg
  integer :: mat_pos_rel(9), mat_pos_rel_indx(9), LU_i_arr(9)
  integer :: elem_n, mat_pos, LU_i, LU_j, incr_LU 
  integer, allocatable, dimension(:) :: LU_colptr_num
  ! 2013-02-21 for subgrid correction


  ! input : col_3_npsi, col_3_ntheta, col_3_npsi_sol, col_3_ntheta_sol, col_period

   
  ! Basic parameter setting
  ! col_3_dt, col_3_sol_solve, col_3_ntotal_r, col_3_ntotal_v
  col_3_dt   =sml_dt*col_period                !forwarded time step (sec) : <pcm_dt> in original version
  if( col_3_npsi_sol * col_3_ntheta_sol .gt. 0 ) then
      col_3_sol_solve = .true.
  else
      col_3_sol_solve = .false.
  endif
  col_3_ntotal_r = col_3_npsi*col_3_ntheta     !Only core simulation
  if(sml_outpsi .gt. eq_x_psi .and. col_3_sol_solve ) &
       col_3_ntotal_r = col_3_ntotal_r + col_3_npsi_sol*col_3_ntheta_sol + 2 !total number of collision cells
                             !Here, +2 means private regions devided into left and right based on x-point line
  col_3_ntotal_v = col_3_nvr  * col_3_nvz



  ! Mesh size setting : col_3_dpsi_core, col_3_dtheta_core, col_3_dpsi_sol, col_3_dtheta_sol, col_3_theta_offset
  if(col_3_npsi .ne. 0) col_3_dpsi_core = (min(eq_x_psi,sml_outpsi)-sml_inpsi)/col_3_npsi
  if(col_3_ntheta .ne. 0) col_3_dtheta_core = sml_2pi/col_3_ntheta
  if(col_3_sol_solve .and. sml_outpsi .gt. eq_x_psi) col_3_dpsi_sol = (sml_outpsi-eq_x_psi)/col_3_npsi_sol
  if(col_3_sol_solve .and. sml_outpsi .gt. eq_x_psi) col_3_dtheta_sol = sml_2pi/col_3_ntheta_sol
  col_3_theta_offset = acos((eq_x_r-eq_axis_r)/sqrt((eq_x_r-eq_axis_r)*(eq_x_r-eq_axis_r)+(eq_x_z-eq_axis_z)*(eq_x_z-eq_axis_z)))



  ! Getting volume of REFERENCE grid : col_3_vol_r
  allocate(col_3_vol_r(0:col_3_ntotal_r-1))
  call get_volume_col3(0,grid)    ! Dependency : col_3_ntotal_r


 
  !====================================SUBGRID SETTING 
  allocate(col_3_subg_num(0:col_3_ntotal_r-1), col_3_refg_indx(0:col_3_ntotal_r-1), &
           col_3_subg_dth(0:col_3_ntotal_r-1) )
  ! col_3_subg_num
  !min_vol = minval(col_3_vol_r)
  min_vol = sum(col_3_vol_r)/(max(1,size(col_3_vol_r)))  !2013-02-20 : I changed this from minimum to mean value
  col_3_subg_num = max(1,nint(col_3_vol_r/min_vol))

  ! col_3_refg_indx (uniform), col_3_total_subg_num 
  col_3_refg_indx(0)=0
  if(col_3_sol_solve) then
      do i=1,col_3_ntotal_r-1-2
          col_3_refg_indx(i) = col_3_refg_indx(i-1) + col_3_subg_num(i-1)
      enddo
      col_3_refg_indx(col_3_ntotal_r-2) = col_3_refg_indx(col_3_ntotal_r-3)+col_3_subg_num(col_3_ntotal_r-3)  !currently, we have fixed 2 private regions for collision
      col_3_refg_indx(col_3_ntotal_r-1) = col_3_refg_indx(col_3_ntotal_r-2)+1
      col_3_total_subg_num = col_3_refg_indx(col_3_ntotal_r-1)+1
  else !core only
      do i=1,col_3_ntotal_r-1
          col_3_refg_indx(i) = col_3_refg_indx(i-1) + col_3_subg_num(i-1)
      enddo
      col_3_total_subg_num = col_3_refg_indx(col_3_ntotal_r-1)+col_3_subg_num(col_3_ntotal_r-1)
  endif
  
  ! col_3_subg_dth 
  col_3_subg_dth(0:col_3_npsi*col_3_ntheta-1) = col_3_dtheta_core / col_3_subg_num(0:col_3_npsi*col_3_ntheta-1)
  col_3_subg_dth(col_3_npsi*col_3_ntheta:col_3_ntotal_r-1) = col_3_dtheta_sol / col_3_subg_num(col_3_npsi*col_3_ntheta:col_3_ntotal_r-1) 

  !=====================================SUBGRID correction due to limiter
  ! Getting volume of SUBgrid : col_3_vol_r
  ! col_3_total_grid_valid_list(new), col_3_total_subg_num(updated) 
  ! dependency : col_3_search_indx_subg <- col_3_subg_dth, col_3_refg_indx, col_3_total_subg_num
  ! NOTE that "col_3_total_subg_num" CHANGES at get_volume_col3
  deallocate(col_3_vol_r)
  allocate(col_3_vol_r(0:col_3_total_subg_num-1))
  allocate(col_3_grid_valid_list(0:col_3_total_subg_num-1))
  call get_volume_col3(1,grid)    ! Dependency : col_3_ntotal_r


 
  ! col_3_nset, col_3_nset_mod, col_3_nset_max, col_3_nset_allocated 
  col_3_nset = col_3_total_subg_num/sml_pe_per_plane   !<colV_num_Vsets> in original version
  col_3_nset_mod = mod(col_3_total_subg_num,sml_pe_per_plane)
  if( col_3_nset_mod .eq. 0) then
     col_3_nset_max = col_3_nset
  else
     col_3_nset_max = col_3_nset + 1
  endif
  if(sml_plane_mype .lt. col_3_nset_mod) then
      col_3_nset_allocated = col_3_nset + 1
  else
      col_3_nset_allocated = col_3_nset
  endif
    

 
  ! col_3_tag
  ! Allocation : collision cells =>> CPU_destination 
  allocate(col_3_tag(0:col_3_nset_allocated-1))
  sep = col_3_nset_max*col_3_nset_mod
  k=0
  do i=0, col_3_total_subg_num-1
     if(i .lt. sep) then
         j = i/col_3_nset_max
     else
         j = col_3_nset_mod + (i-sep)/col_3_nset
     endif

     if(sml_plane_mype .eq. j) then
         col_3_tag(k) = i
         k=k+1
      endif
  enddo


  ! col_3d_scounts, col_3d_sdispls, col_3d_rcounts, col_3d_rdispls
  allocate(col_3d_scounts(0:sml_pe_per_plane-1), col_3d_sdispls(0:sml_pe_per_plane-1), &
           col_3d_rcounts(0:sml_pe_per_plane-1), col_3d_rdispls(0:sml_pe_per_plane-1))
  ! below routine is time independent. i.e. move this to init_col with global memory
  col_3d_scounts = 0
  col_3d_sdispls = 0
  sep = col_3_nset_max*col_3_nset_mod
  do i=0, sml_pe_per_plane-2
     if(i .lt. col_3_nset_mod) then
         col_3d_scounts(i) = col_3_nset_max
         col_3d_sdispls(i+1) = col_3d_sdispls(i)+col_3_nset_max !sdispls(0)=0
     else
         col_3d_scounts(i) = col_3_nset
         col_3d_sdispls(i+1) = col_3d_sdispls(i)+col_3_nset
     endif
  enddo
  col_3d_scounts(sml_pe_per_plane-1) = col_3_nset
  
  col_3d_rcounts = 0
  col_3d_rdispls = 0
  if(sml_plane_mype .lt. col_3_nset_mod) then
      col_3d_rcounts(:) = col_3_nset_max
      do i=0, sml_pe_per_plane-2
          col_3d_rdispls(i+1) = col_3d_rdispls(i)+col_3_nset_max
      enddo
  else 
      col_3d_rcounts(:) = col_3_nset
      do i=0, sml_pe_per_plane-2
          col_3d_rdispls(i+1) = col_3d_rdispls(i)+col_3_nset
      enddo
  endif
        
  ! 2013-02-23 =============================================== SUPER LU
  ! index_map_LU, LU_rowindx, LU_cvalues
  LU_n = col_3_ntotal_v !global
  LU_nnz = 4*4 + 6*2*(col_3_nvr-2)+6*2*(col_3_nvz-2)+9*(col_3_nvr-2)*(col_3_nvz-2) !global 
  LU_nrhs = 1 !global
  LU_ldb = LU_n !global
  allocate(LU_rowindx(LU_nnz), LU_cvalues(LU_n), LU_colptr(LU_n+1), index_map_LU(9,LU_n))       !global
  allocate(LU_colptr_num(LU_n))    !local
  LU_colptr_num = 0   ! number of elements in each column, local
  LU_rowindx = 0         ! global

  !below is time independent. Move to init_col in module
  LU_colptr(1) = 1
  do i=1,LU_n
      !for colptr
      LU_i = (i-1) / col_3_nvr+1
      LU_j = mod(i-1, col_3_nvr)+1
      if( (LU_i .eq. 1) .or. (LU_i .eq. col_3_nvz) ) then
          if( (LU_j .eq. 1) .or. (LU_j .eq. col_3_nvr) ) then
              incr_LU = 4
          else
              incr_LU = 6
          endif
      else
          if( (LU_j .eq. 1) .or. (LU_j .eq. col_3_nvr) ) then
              incr_LU = 6
          else
              incr_LU = 9
          endif
      endif
      LU_colptr(i+1) = LU_colptr(i)+incr_LU
  enddo 

  !===============
  !  3--6--9
  !  |  |  |
  !  2--5--8
  !  |  |  |
  !  1--4--7
  !==============
  mat_pos_rel=(/-col_3_nvr-1,-col_3_nvr, -col_3_nvr+1, -1, 0, 1, col_3_nvr-1, col_3_nvr, col_3_nvr+1/) 
  index_map_LU = 0
  do i=1,col_3_nvz
      do j=1, col_3_nvr
          mat_pos = j+(i-1)*col_3_nvr

          ! INDEXING
          if(i .eq. 1) then
              if(j .eq. 1) then
                  !(J,I)=(1,1)
                  elem_n = 4
                  mat_pos_rel_indx(1:elem_n) =(/5,6,8,9/)
               elseif (j .eq. col_3_nvr) then
                  !(J,I)=(Nr,1)
                  elem_n = 4
                  mat_pos_rel_indx(1:elem_n)=(/4,5,7,8/)
               else
                  !(J,I)=(:,1)
                  elem_n = 6
                  mat_pos_rel_indx(1:elem_n)=(/4,5,6,7,8,9/)
               endif
           elseif(i .eq. col_3_nvz) then
               if(j .eq. 1) then
                  !(J,I)=(1,mesh_Nz)
                  elem_n = 4
                  mat_pos_rel_indx(1:elem_n)=(/2,3,5,6/)
               elseif (j .eq. col_3_nvr) then
                  !(J,I)=(Nr,mesh_Nz)
                  elem_n = 4
                  mat_pos_rel_indx(1:elem_n)=(/1,2,4,5/)
               else
                  !(J,I)=(:,mesh_Nz)
                  elem_n = 6
                  mat_pos_rel_indx(1:elem_n)=(/1,2,3,4,5,6/)
               endif
           else
               if(j .eq. 1) then
                  !(J,I) = (1,:)
                  elem_n = 6
                  mat_pos_rel_indx(1:elem_n)=(/2,3,5,6,8,9/)
               elseif(j .eq. col_3_nvr) then
                  !(J,I) = (mesh_Nr,:)
                  elem_n = 6
                  mat_pos_rel_indx(1:elem_n)=(/1,2,4,5,7,8/)
               else
                  !(J,I) = (:,:)
                  elem_n = 9
                  mat_pos_rel_indx(1:elem_n)=(/1,2,3,4,5,6,7,8,9/)
               endif
           endif
           LU_i_arr(1:elem_n) = mat_pos+mat_pos_rel(mat_pos_rel_indx(1:elem_n))  ! I need to change LU_i to array
           index_map_LU(mat_pos_rel_indx(1:elem_n),mat_pos) = LU_colptr(LU_i_arr(1:elem_n))+LU_colptr_num(LU_i_arr(1:elem_n))
           LU_colptr_num(LU_i_arr(1:elem_n)) = LU_colptr_num(LU_i_arr(1:elem_n))+1
           LU_rowindx(index_map_LU(mat_pos_rel_indx(1:elem_n),mat_pos)) = mat_pos


           LU_cvalues(mat_pos) = index_map_LU(5, mat_pos)  !For implicit time marching
      enddo
  enddo
  deallocate(LU_colptr_num)      
#ifdef COL_3D_DEBUG
        ! collision grid output
        if(sml_mype .eq. 0) then
            rdim=sml_bd_max_r - sml_bd_min_r
            roffset=sml_bd_min_r
            zdim=sml_bd_max_z - sml_bd_min_z
            zoffset=sml_bd_min_z

            resol_r = 1000
            resol_z = 1000

            spec_dr = rdim/resol_r
            spec_dz = zdim/resol_z

            print *, 'printing col grid'

            open(unit=908, status='replace')

            do i=0, resol_r-1 ! arbitrary values of mine for resolution
                spec_r = roffset + i*spec_dr
                do j=0, resol_z-1
                    spec_z = zoffset + j*spec_dz
                    spec_psi = psi_interpol(spec_r, spec_z, 0, 0)
                    spec_theta= acos((spec_r-eq_axis_r)/sqrt((spec_r-eq_axis_r)*(spec_r-eq_axis_r)+(spec_z-eq_axis_z)*(spec_z-eq_axis_z)))
                    if(spec_z .lt. eq_axis_z) then
                        spec_theta= sml_2pi-spec_theta
                    endif

                    !indx = col_3_search_indx(spec_r,spec_z,spec_psi,spec_theta) !for reference grid
                    !indx = col_3_search_indx_subg(spec_r,spec_z,spec_psi,spec_theta)
                    indx = col_3_search_indx_subg(spec_r,spec_z,spec_psi,spec_theta,grid,2)
                    write(908,1002) spec_r, spec_z, indx, spec_psi

                enddo
            enddo
            close(908)
1002 format(2(e19.13,1x),I8,1x,e19.13)  
        endif
#endif
 
end subroutine col_3_setup

! dtheta_sol, ntheta_sol, dpsi_sol, npsi, ntheta, ntheta_sol
! col_3_refg_indx, col_3_subg_dth
integer function col_3_search_indx_subg(r, z, psi, theta, grid, mode)
    use grid_class
    use col_module
    use eq_module, only : eq_x_psi, is_rgn1
    use sml_module, only : sml_2pi, sml_pi, sml_inpsi, sml_outpsi
    implicit none
    integer ::  indx
    real (8) :: r, z, psi, theta, theta_off
    type(grid_type), target :: grid
    integer :: ierror, th_guide, th_subg, ipsi, iref
    !2013-02-20 for domain limiting
    real (8) :: x(2), phi, phi_mid, p(3)
    integer :: itr
    !2013-02-21
    integer :: mode,i
    logical :: domain_in

! mode = 0 : check validity with triangle - output : index based on angle
! mode = 1 : check validity with col_3_grid_valid_list which is usuable after "get_volume_col3"
!            - output : index based on angle
! mode = 2 : check validity with col_3_grid_valid_list which is usuable after "get_volume_col3"
!            - output : index for array (ubound = col_3_total_subg_num, which is effective number)
    
    if(col_3_sol_solve) then
       x(1)=r
       x(2)=z
       call search_tr2(grid,x,itr,p)  !We don't need field_following to check just if particle is on triangle
    else
       itr=1  !just positive is enough 
    endif
    
    theta_off = mod(theta+col_3_theta_offset,sml_2pi)
    
    if(sml_inpsi<psi .and. psi<sml_outpsi .and. itr .ge. 0) then
        if(psi .gt. eq_x_psi) then
            !region 2 - SOL
            th_guide = mod(floor(theta_off/col_3_dtheta_sol),col_3_ntheta_sol)
            ipsi=floor((psi-eq_x_psi)/col_3_dpsi_sol)

            iref = col_3_npsi*col_3_ntheta+ipsi*col_3_ntheta_sol+th_guide 

            th_subg = floor((theta_off-th_guide*col_3_dtheta_sol)/col_3_subg_dth(iref))
            indx=col_3_refg_indx(iref)+th_subg
        elseif( is_rgn1(r,z,psi) ) then
            !region 1 - core
            th_guide = mod(floor(theta_off/col_3_dtheta_core),col_3_ntheta)
            ipsi =floor((psi-sml_inpsi)/col_3_dpsi_core)
            
            iref = ipsi*col_3_ntheta+th_guide
              
            th_subg = floor((theta_off-th_guide*col_3_dtheta_core)/col_3_subg_dth(iref))
            indx=col_3_refg_indx(iref)+th_subg 
        else
            !region 3 - private
            indx = col_3_refg_indx(col_3_npsi*col_3_ntheta+col_3_npsi_sol*col_3_ntheta_sol)
            if(theta_off .lt. sml_pi) indx = indx+1
        endif
        if(mode .eq. 0) call assert(indx .lt. col_3_total_subg_num .and. indx .ge. 0,'col_3 indexing problem in subgrid',ierror)
    else
        indx=-1
    endif

    if(mode .ne. 0) then
        domain_in = .false.
        do i=0, col_3_total_subg_num-1
            if(indx .eq. col_3_grid_valid_list(i)) then
                 domain_in = .true.
                 if(mode .eq. 2) indx = i
                 exit
            endif
        enddo
        if(.not. domain_in) indx = -1
    endif

    col_3_search_indx_subg= indx

   
end function col_3_search_indx_subg

integer function col_3_search_indx_subg_decomp(r, z, psi, theta)
    use col_module
    use eq_module, only : eq_x_psi, is_rgn1
    use sml_module, only : sml_2pi, sml_pi, sml_inpsi, sml_outpsi
    implicit none
    integer ::  indx
    real (8) :: r, z, psi, theta, theta_off
    integer :: ierror, th_guide, th_subg, ipsi, iref

    theta_off = mod(theta+col_3_theta_offset,sml_2pi)
    
    if(sml_inpsi<psi .and. psi<sml_outpsi) then
        if(psi .gt. eq_x_psi) then
            !region 2 - SOL
            th_guide = mod(floor(theta_off/col_3_dtheta_sol),col_3_ntheta_sol)
            ipsi=floor((psi-eq_x_psi)/col_3_dpsi_sol)

            iref = col_3_npsi*col_3_ntheta+ipsi*col_3_ntheta_sol+th_guide 

            th_subg = floor((theta_off-th_guide*col_3_dtheta_sol)/col_3_subg_dth(iref))
            indx=col_3_refg_indx(iref)+th_subg
        elseif( is_rgn1(r,z,psi) ) then
            !region 1 - core
            th_guide = mod(floor(theta_off/col_3_dtheta_core),col_3_ntheta)
            ipsi =floor((psi-sml_inpsi)/col_3_dpsi_core)
            
            iref = ipsi*col_3_ntheta+th_guide
              
            th_subg = floor((theta_off-th_guide*col_3_dtheta_core)/col_3_subg_dth(iref))
            indx=col_3_refg_indx(iref)+th_subg 
        else
            !region 3 - private
            indx = col_3_refg_indx(col_3_npsi*col_3_ntheta+col_3_npsi_sol*col_3_ntheta_sol)
            if(theta_off .lt. sml_pi) indx = indx+1
        endif
        call assert(indx .lt. col_3_total_subg_num .and. indx .ge. 0,'col_3 indexing problem in subgrid',ierror)
    else
        indx=-1
    endif

    col_3_search_indx_subg_decomp = indx
end function col_3_search_indx_subg_decomp

integer function col_3_search_indx(r, z, psi, theta, grid)
    use grid_class
    use sml_module, only : sml_2pi, sml_pi, sml_inpsi, sml_outpsi
    use eq_module, only : eq_x_psi, is_rgn1
    use col_module
    implicit none
    integer ::  indx
    real (8) :: r, z, psi, theta, theta_off
    type(grid_type) :: grid
    integer :: ierror, th_guide, ipsi
    !2013-02-20 for domain limiting
    real (8) :: x(2), phi, phi_mid, p(3)
    integer :: itr

    if(col_3_sol_solve) then
       x(1)=r
       x(2)=z
       call search_tr2(grid,x,itr,p)  !We don't need field_following to check just if particle is on triangle
    else
       itr=1  !just positive is enough 
    endif

    theta_off = mod(theta+col_3_theta_offset,sml_2pi)
    
    if(sml_inpsi<psi .and. psi<sml_outpsi .and. itr .ge. 0) then
        if(psi .gt. eq_x_psi) then
            !region 2 - SOL
            th_guide = mod(floor(theta_off/col_3_dtheta_sol),col_3_ntheta_sol)
            ipsi=floor((psi-eq_x_psi)/col_3_dpsi_sol)

            indx=col_3_npsi*col_3_ntheta+ipsi*col_3_ntheta_sol+th_guide
        elseif( is_rgn1(r,z,psi) ) then
            !region 1 - core
            th_guide = mod(floor(theta_off/col_3_dtheta_core),col_3_ntheta)
            ipsi =floor((psi-sml_inpsi)/col_3_dpsi_core)  

            indx=ipsi*col_3_ntheta+th_guide 
        else
            !region 3 - private
            indx = col_3_npsi*col_3_ntheta+col_3_npsi_sol*col_3_ntheta_sol
            if(theta_off .lt. sml_pi) indx = indx+1
        endif

        call assert(indx .lt. col_3_ntotal_r .and. indx .ge. 0,'col_3 indexing problem in refs grid',ierror)
    else
        indx=-1
    endif

    col_3_search_indx = indx
end function col_3_search_indx


