!*************************************************
! diagnosis main
!************************************************
#include "adios_macro.h"

subroutine diagnosis(istep,irk,grid,psn,spall)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use ptl_module
  use perf_monitor
  implicit none
  integer, intent(in) :: istep,irk
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  !
  integer :: i

  select case( irk )
  case(1) ! IRK = 1

     !for time average variables
     if(diag_tavg_on) then
        diag_1d_tavg_f_pv1=diag_1d_tavg_f_pv1 + diag_1d_f_pv1
        if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
           diag_1d_tavg_df_pv1=diag_1d_tavg_df_pv1 + diag_1d_df_pv1
        endif
     endif

     ! output 1d (+ tavg) variables
     if( mod(istep,diag_1d_period)==0 .and. diag_1d_on) then
        call t_startf("DIAG_1D_OUTPUT")
        call diag_1d_output(grid,psn)
        if(diag_heat_on) call diag_heat_output
        call t_stopf("DIAG_1D_OUTPUT")
        call diag_1d_tavg_port_clear
     endif
     call diag_1d_port_clear  ! requires for tavg.


     if(mod(istep,diag_3d_period)==0 .and. diag_3d_on) then
        call t_startf("DIAG_3D")
        call diag_3d(istep,grid,psn)   ! (some) MPI, I/O
        !if(sml_f0_grid) call diag_f0(istep,grid,psn) !diag_f0 will be called in main
        call t_stopf("DIAG_3D")
     endif

     ! tracer --------------------
     if(diag_tracer_n/=0) then
        if(mod(istep,diag_tracer_period)==0) then
           call t_startf("TRACER")
           ! single particle, (small) computation, I/O
           call tracer(diag_tracer_n,diag_tracer_sp,grid,psn,spall)
           call t_stopf("TRACER")
        endif
     endif

     ! particle dump --------------
     if(diag_particle_mod/=0) then
       if(mod(istep,diag_particle_period)==0) then
          call t_startf("DIAG_PARTICLE")
          call diag_particle(grid,psn,spall)
          call t_stopf("DIAG_PARTICLE")
       endif
      endif

  case(2) ! IRK = 2

  end select


end subroutine diagnosis

subroutine determine_diag_on(istep,ipc,diag_on)
  use diag_module
  implicit none
  integer :: istep, ipc
  logical :: diag_on

  diag_on=.false.
  if(ipc==1) then
     if(diag_tavg_on) then
        diag_on=.true.
     else
        if( mod(istep,diag_1d_period)==0 ) then
           diag_on=.true.
        endif
     endif
  endif

end subroutine determine_diag_on

! port1 called in push
subroutine diag_1d_port1(ptli,derivs,sp_type,vd,ith)
  use ptl_module
  use diag_module
  use sml_module, only : sml_n_vf_diag, sml_deltaf, sml_ev2j, sml_mype, sml_inpsi
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z, eq_axis_r, eq_axis_z,eq_ftn, eq_tempi
  implicit none
  type(ptl_type), intent(in) :: ptli
  real (8), intent(in) :: derivs(ptl_nphase)
  integer, intent(in) :: sp_type
  real (8) :: vd(sml_n_vf_diag) ! 1: R_major, 2: B_toroidal, 3: B_total  4: radial ExB 5: v_para

  integer, intent(in) :: ith
  !
  real (8) :: psi,z, pn, wp, b, r, rho, mu, w, dw, vp, en, den, pe, we
  real (8) ::  diag_1d_de
  integer  :: ip, j, ie
  real (8) :: v(diag_1d_npv1)


  if(ptli%gid > 0) then
     psi = vd(1)
     r=ptli%ph(1)
     z=ptli%ph(2)
     if(.not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) return ! EXIT for private region of lower divertor
     pn=(psi-diag_1d_pin)*diag_1d_dp_inv
     ip=floor(pn)+1
     if(ip <1 .or. diag_1d_npsi <= ip) return     ! EXIT for out of diag_1d_pin/pout range
     wp=1D0 - pn + real(ip-1,8)

     ! local variables for readability
     b=vd(2)

     rho=ptli%ph(4)
     mu=ptli%ct(1)
     w=ptli%ct(2)
     dw=ptli%ph(5)*w
     vp=ptl_c_m(sp_type)*rho*B

     ! obtain variables
     v(1) = 1D0                     ! full weight
     v(2) = derivs(3)*r             ! g.c. toroidal velocity --  dzeta/dt * R_major
     v(3) = vd(7)                   ! g.c. poloidal velocity -- v dot Bp/Bp
     v(4) = vp                      ! v_|| parallel velocity
     v(5) = vp*r*vd(3)/vd(2)        ! v_zeta_|| / R --toroidal angular momentum of v_||   : v_|| * B_zeta /(R*B)
     v(6) = vd(5)                   ! radial drift  - psi_dot
     v(7) = ptl_c2_2m(sp_type)*(rho*B)**2 ! parallel mean energy
     v(8) = mu*B                    ! perp temperature
     v(9) = (v(7)+v(8))*v(6)        ! Energy radial flux
     v(10)= v(4)*vd(4)*vd(3)/vd(2)  ! V_exb * V_|| * B_phi/B
     v(11)= v(4)*vd(5)*vd(3)/vd(2)  ! V_r   * V_|| * B_phi/B
     v(12)= vd(4)                   ! radial drift by exb  - V_exb dot grad_psi
     v(13)= (v(7)+v(8))*vd(4)       ! heat flux by radial exb
     v(14)= vd(8)                   ! poloidal comp. of V_ExB
!     v(15)= vd(6)                   ! grad_psi ^2
     v(15)= sqrt(vd(6))                   ! grad_psi 
     diag_1d_f_pv1(:,ip  ,sp_type,ith)=v(:)*w*wp       +diag_1d_f_pv1(:,ip  ,sp_type,ith)
     diag_1d_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1D0-wp) +diag_1d_f_pv1(:,ip+1,sp_type,ith)

     if(ptl_deltaf_sp(sp_type)) then
        diag_1d_df_pv1(:,ip  ,sp_type,ith)=v(:)*dw*wp       +diag_1d_df_pv1(:,ip  ,sp_type,ith)
        diag_1d_df_pv1(:,ip+1,sp_type,ith)=v(:)*dw*(1D0-wp) +diag_1d_df_pv1(:,ip+1,sp_type,ith)
     endif

     if(diag_tavg_on .and. diag_omid_on) then
        if(  r-eq_axis_r > abs(z-eq_axis_z) ) then
             diag_1d_omid_f_pv1(:,ip  ,sp_type,ith)=v(:)*w*wp       +diag_1d_omid_f_pv1(:,ip  ,sp_type,ith)
             diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)=v(:)*w*(1D0-wp) +diag_1d_omid_f_pv1(:,ip+1,sp_type,ith)
        endif
     endif
     if(diag_eflux_on)then
        diag_1d_emin = 0D0
        diag_1d_emax = 3D0*eq_ftn(sml_inpsi,eq_axis_r,eq_axis_z,eq_tempi)*sml_ev2j
        diag_1d_de=(diag_1d_emax-diag_1d_emin)/diag_1d_ne

        en=v(7)+v(8)
        pe=(en-diag_1d_emin)/diag_1d_de
        ie=floor(pe)+1
     if(ie >=1 .and. diag_1d_ne > ie) then     ! EXIT for out of diag_1d_emin/emax range
        we=1D0 - pe + real(ie-1,8)
        diag_1d_eflux_pv(1,ip  ,ie  ,sp_type,ith)=w*wp*we            +diag_1d_eflux_pv(1,ip  ,ie  ,sp_type,ith)
        diag_1d_eflux_pv(1,ip+1,ie  ,sp_type,ith)=w*(1D0-wp)*we      +diag_1d_eflux_pv(1,ip+1,ie  ,sp_type,ith)
        diag_1d_eflux_pv(1,ip  ,ie+1,sp_type,ith)=w*wp*(1D0-we)      +diag_1d_eflux_pv(1,ip  ,ie+1,sp_type,ith)
        diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith)=w*(1D0-wp)*(1D0-we)+diag_1d_eflux_pv(1,ip+1,ie+1,sp_type,ith)

        diag_1d_eflux_pv(2,ip  ,ie  ,sp_type,ith)=v(9)*w*wp*we            +diag_1d_eflux_pv(2,ip  ,ie  ,sp_type,ith)
        diag_1d_eflux_pv(2,ip+1,ie  ,sp_type,ith)=v(9)*w*(1D0-wp)*we      +diag_1d_eflux_pv(2,ip+1,ie  ,sp_type,ith)
        diag_1d_eflux_pv(2,ip  ,ie+1,sp_type,ith)=v(9)*w*wp*(1D0-we)      +diag_1d_eflux_pv(2,ip  ,ie+1,sp_type,ith)
        diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith)=v(9)*w*(1D0-wp)*(1D0-we)+diag_1d_eflux_pv(2,ip+1,ie+1,sp_type,ith)
        if(ptl_deltaf_sp(sp_type)) then
           diag_2d_dflux_pv(1,ip  ,ie  ,sp_type,ith)=dw*wp*we            +diag_2d_dflux_pv(1,ip  ,ie  ,sp_type,ith)
           diag_2d_dflux_pv(1,ip+1,ie  ,sp_type,ith)=dw*(1D0-wp)*we      +diag_2d_dflux_pv(1,ip+1,ie  ,sp_type,ith)
           diag_2d_dflux_pv(1,ip  ,ie+1,sp_type,ith)=dw*wp*(1D0-we)      +diag_2d_dflux_pv(1,ip  ,ie+1,sp_type,ith)
           diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith)=dw*(1D0-wp)*(1D0-we)+diag_2d_dflux_pv(1,ip+1,ie+1,sp_type,ith)

           diag_2d_dflux_pv(2,ip  ,ie  ,sp_type,ith)=v(9)*dw*wp*we            +diag_2d_dflux_pv(2,ip  ,ie  ,sp_type,ith)
           diag_2d_dflux_pv(2,ip+1,ie  ,sp_type,ith)=v(9)*dw*(1D0-wp)*we      +diag_2d_dflux_pv(2,ip+1,ie  ,sp_type,ith)
           diag_2d_dflux_pv(2,ip  ,ie+1,sp_type,ith)=v(9)*dw*wp*(1D0-we)      +diag_2d_dflux_pv(2,ip  ,ie+1,sp_type,ith)
           diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith)=v(9)*dw*(1D0-wp)*(1D0-we)+diag_2d_dflux_pv(2,ip+1,ie+1,sp_type,ith)
        endif
     endif ! ie
     endif ! diag_eflux_on
  endif

end subroutine diag_1d_port1

subroutine diag_1d_init
  use diag_module
  use sml_module, only : sml_nthreads, sml_deltaf, sml_deltaf_elec, sml_electron_on, sml_f0_grid
  use ptl_module, only : ptl_isp, ptl_nsp
  implicit none
  integer :: np, isp, nsp    ! # of radial grid, species initial, species end
  real(8) :: pmin, pmax  ! minimum q, maximum q
  integer :: i

  diag_1d_isp=ptl_isp
  diag_1d_nsp=ptl_nsp

  np=diag_1d_npsi
  isp=diag_1d_isp
  nsp=diag_1d_nsp

  diag_1d_dp=(diag_1d_pout-diag_1d_pin)/real(np-1,8)
  diag_1d_dp_inv=1D0/diag_1d_dp

  if (allocated(diag_1d_f_pv1)) deallocate(diag_1d_f_pv1)
  if (allocated(diag_1d_vol)) deallocate(diag_1d_vol)
  allocate(diag_1d_f_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads), diag_1d_vol(np))
  diag_1d_f_pv1=0D0

  if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
     if (allocated(diag_1d_df_pv1)) deallocate(diag_1d_df_pv1)
     allocate(diag_1d_df_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads))
     diag_1d_df_pv1=0D0

     if(sml_f0_grid) then
        if(.not. allocated(diag_1d_f0_pv1)) then
           allocate(diag_1d_f0_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads))
           diag_1d_f0_pv1=0D0
        endif
        if(.not. allocated(diag_2d_ef0_pv)) then
           allocate(diag_2d_ef0_pv(2,np,diag_1d_ne,isp:nsp,sml_nthreads))
           diag_2d_ef0_pv=0D0
        endif
     endif
  endif

  if(diag_tavg_on) then
     if (allocated(diag_1d_tavg_f_pv1)) deallocate(diag_1d_tavg_f_pv1)
     allocate(diag_1d_tavg_f_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads))
     diag_1d_tavg_f_pv1=0D0
  endif

  if(diag_tavg_on .and. (sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on))) then
     if (allocated(diag_1d_tavg_df_pv1)) deallocate(diag_1d_tavg_df_pv1)
     allocate(diag_1d_tavg_df_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads))
     diag_1d_tavg_df_pv1=0D0
  endif

  if(diag_omid_on .and. diag_tavg_on) then
     if (allocated(diag_1d_omid_f_pv1)) deallocate(diag_1d_omid_f_pv1)
     allocate(diag_1d_omid_f_pv1(diag_1d_npv1,np,isp:nsp,sml_nthreads))
     diag_1d_omid_f_pv1=0D0
  endif

  if(diag_eflux_on) then
     if (allocated(diag_1d_eflux_pv)) deallocate(diag_1d_eflux_pv)
     allocate(diag_1d_eflux_pv(2,np,diag_1d_ne,isp:nsp,sml_nthreads))
     diag_1d_eflux_pv=0D0
     if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
       if (allocated(diag_2d_dflux_pv)) deallocate(diag_2d_dflux_pv)
       allocate(diag_2d_dflux_pv(2,np,diag_1d_ne,isp:nsp,sml_nthreads))
       diag_2d_dflux_pv=0D0
     endif
  endif

  if(diag_f0_df_on) then
    if (allocated(diag_f0_df_pv1)) deallocate(diag_f0_df_pv1)
    allocate(diag_f0_df_pv1(diag_f0_df_npv1,np,isp:nsp,diag_f0_df_nsource))
    diag_f0_df_pv1=0D0
  endif
end subroutine diag_1d_init

subroutine diag_1d_port_clear
  use diag_module
  use sml_module, only : sml_deltaf, sml_deltaf_elec, sml_electron_on
  implicit none

  diag_1d_f_pv1=0D0
  if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) diag_1d_df_pv1=0D0
  if(diag_eflux_on) then
    diag_1d_eflux_pv=0D0
    if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) diag_2d_dflux_pv=0D0
  endif
end subroutine diag_1d_port_clear

subroutine diag_1d_tavg_port_clear
  use diag_module
  use sml_module
  implicit none

  if(diag_tavg_on) diag_1d_tavg_f_pv1=0D0
  if((sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) .and. diag_tavg_on) diag_1d_tavg_df_pv1=0D0
  if(diag_omid_on)  diag_1d_omid_f_pv1=0D0
  if(diag_eflux_on) then
    diag_1d_eflux_pv=0D0
    if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) diag_2d_dflux_pv=0D0
  endif
end subroutine diag_1d_tavg_port_clear

subroutine diag_f0_df_port_clear
  use diag_module
  use sml_module
  implicit none

  if(diag_f0_df_on) then
    diag_f0_df_pv1=0D0
  endif

end subroutine


subroutine diag_1d_output(grid,psn)
  use grid_class
  use psn_class
  use ptl_module
  use eq_module
  use sml_module
  use diag_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: ip, isp, j
  character (len=256) :: filename
  include 'mpif.h'
  
  real (8) :: psi
  real (8), allocatable :: out1(:,:,:)
  real (8), allocatable :: out1_tavg(:,:,:)
  real (8), allocatable :: out1_omid(:,:,:)
  real (8), allocatable :: out1_eflux(:,:,:)
  real (8), allocatable :: out1_pflux(:,:,:)
  real (8), allocatable :: out1_weight(:,:,:)
  real (8), allocatable :: out1_test(:,:)
  real (8) :: psi00(grid%npsi00)
  real (8) :: dum(diag_1d_npv1,diag_1d_npsi,ptl_isp:ptl_nsp)
  real (8) :: dum1(2,diag_1d_npsi,diag_1d_ne,ptl_isp:ptl_nsp)
  real (8) :: unit1(diag_1d_npv1)
  real (8) ::  diag_1d_de, tmp1, tmp2
  integer :: ith, np, nv1, nv_total, ne, ie
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
!  real (8) :: tmp2(diag_1d_npsi,diag_1d_npv2)
  real (8), allocatable :: tmp(:,:,:), tmp_tavg(:,:,:), tmp_omid(:,:,:)
!  real (8), allocatable :: tmp_eflux(:,:,:), tmp_pflux(:,:,:), tmp_weight(:,:,:)
  real (8) :: pall(diag_1d_npsi), pnorm(diag_1d_npsi)
  character (len=2) :: sp_name(0:6)=(/'e_', 'i_', 'i2', 'i3', 'i4', 'i5', 'i6'/)  ! max 6 ion species
  !character (len=64) :: var_names2(diag_1d_npv2)
  character (len=64), allocatable :: var_names(:)
!  character (len=3), parameter :: postfix='_1d'
#endif


  np=diag_1d_npsi
  nv1=diag_1d_npv1
  ne=diag_1d_ne

  if(.not. (sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on))) then
     nv_total=nv1
  else
     nv_total=nv1*2
  endif

  allocate(out1(nv_total,np,ptl_isp:ptl_nsp))
  if(diag_tavg_on) then
     allocate(out1_tavg(nv_total,np,ptl_isp:ptl_nsp))
  endif
 if(diag_tavg_on .and. diag_omid_on) then
     allocate(out1_omid(nv_total,np,ptl_isp:ptl_nsp))
  endif
  if(diag_eflux_on) then
     diag_1d_emin = 0D0
     diag_1d_emax = 3D0*eq_ftn(sml_inpsi,eq_axis_r,eq_axis_z,eq_tempi)*sml_ev2j
     diag_1d_de=(diag_1d_emax-diag_1d_emin)/diag_1d_ne
     allocate(out1_eflux(np,ne,ptl_isp:ptl_nsp))
     allocate(out1_pflux(np,ne,ptl_isp:ptl_nsp))
     allocate(out1_weight(np,ne,ptl_isp:ptl_nsp))
     allocate(out1_test(np,ptl_isp:ptl_nsp))
  endif

#ifdef ADIOS

  if(sml_f0_grid) call diag_1d_f0(grid,psn)


  allocate(tmp(np,nv_total,ptl_isp:ptl_nsp),var_names(nv_total))
  if(diag_tavg_on) allocate(tmp_tavg(np,nv_total,ptl_isp:ptl_nsp))
  if(diag_tavg_on .and. diag_omid_on) allocate(tmp_omid(np,nv_total,ptl_isp:ptl_nsp))
!  if(diag_eflux_on) allocate(tmp_eflux(np,ne,ptl_isp:ptl_nsp),tmp_pflux(np,ne,ptl_isp:ptl_nsp),&
!     tmp_weight(np,ne,ptl_isp:ptl_nsp))
  var_names(1) ='gc_density'
  var_names(2) ='gc_toroidal_flow'
  var_names(3) ='gc_poloidal_flow'
  var_names(4) ='parallel_flow'
  var_names(5) ='tor_angular_momentum'
  var_names(6) ='radial_flux'
  var_names(7) ='parallel_mean_en'
  var_names(8) ='perp_temperature'
  var_names(9) ='radial_en_flux'
  var_names(10)='radial_mom_flux_ExB'
  var_names(11)='radial_mom_flux'
  var_names(12)='radial_flux_ExB'
  var_names(13)='radial_en_flux_ExB'
  var_names(14)='poloidal_ExB_flow'
  var_names(15)='grad_psi_sqr'

  if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
     do j=1, diag_1d_npv1
        var_names(nv1+j) = trim(var_names(j))//'_df'
     enddo
  endif

 ! var_names2(1)='charge_density00'
#endif


  unit1=1D0;
  unit1(7:8)=sml_j2ev


  ! open mp reduce & mpi reduce
  if (sml_nthreads >= 2) then
  do ith=2, sml_nthreads
     diag_1d_f_pv1(:,:,:,1)=diag_1d_f_pv1(:,:,:,1)+diag_1d_f_pv1(:,:,:,ith)
  enddo
  endif
  call my_mpi_reduce(diag_1d_f_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
  diag_1d_f_pv1(:,:,:,1)=dum

  !delta-f
  if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then

     if (sml_nthreads >= 2) then
     do ith=2, sml_nthreads
        diag_1d_df_pv1(:,:,:,1)=diag_1d_df_pv1(:,:,:,1)+diag_1d_df_pv1(:,:,:,ith)
     enddo
     endif

     call my_mpi_reduce(diag_1d_df_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
     diag_1d_df_pv1(:,:,:,1)=dum
  endif

  !**** time average output ****
  if(diag_tavg_on) then
     if (sml_nthreads >= 2) then
     do ith=2, sml_nthreads
        diag_1d_tavg_f_pv1(:,:,:,1)=diag_1d_tavg_f_pv1(:,:,:,1)+diag_1d_tavg_f_pv1(:,:,:,ith)
     enddo
     endif
     call my_mpi_reduce(diag_1d_tavg_f_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
     diag_1d_tavg_f_pv1(:,:,:,1)=dum

     if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
        if (sml_nthreads >= 2) then
        do ith=2, sml_nthreads
           diag_1d_tavg_df_pv1(:,:,:,1)=diag_1d_tavg_df_pv1(:,:,:,1)+diag_1d_tavg_df_pv1(:,:,:,ith)
        enddo
        endif

        call my_mpi_reduce(diag_1d_tavg_df_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
        diag_1d_tavg_df_pv1(:,:,:,1)=dum
     endif
  endif

  !**** time average outside midplane output ****
  if(diag_tavg_on .and. diag_omid_on) then
     if (sml_nthreads >= 2) then
     do ith=2, sml_nthreads
        diag_1d_omid_f_pv1(:,:,:,1)=diag_1d_omid_f_pv1(:,:,:,1)+diag_1d_omid_f_pv1(:,:,:,ith)
     enddo
     endif
     call my_mpi_reduce(diag_1d_omid_f_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
     diag_1d_omid_f_pv1(:,:,:,1)=dum

!     if(sml_deltaf) then
!        do ith=2, sml_nthreads
!           diag_1d_tavg_df_pv1(:,:,:,1)=diag_1d_tavg_df_pv1(:,:,:,1)+diag_1d_tavg_df_pv1(:,:,:,ith)
!        enddo

!        call my_mpi_reduce(diag_1d_tavg_df_pv1(:,:,:,1),dum,nv1*np*(ptl_nsp-ptl_isp+1))
!        diag_1d_tavg_df_pv1(:,:,:,1)=dum
!     endif

  endif

  !**** Include kinetic ion particle flux ****
  if(diag_eflux_on) then
     if (sml_nthreads >= 2) then
     do ith=2, sml_nthreads
        diag_1d_eflux_pv(:,:,:,:,1)=diag_1d_eflux_pv(:,:,:,:,1)+diag_1d_eflux_pv(:,:,:,:,ith)
     enddo
     endif
     call my_mpi_reduce( diag_1d_eflux_pv(:,:,:,:,1),dum1,2*ne*np*(ptl_nsp-ptl_isp+1))
     diag_1d_eflux_pv(:,:,:,:,1)=dum1
     if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
       if (sml_nthreads >= 2) then
       do ith=2, sml_nthreads
          diag_2d_dflux_pv(:,:,:,:,1)=diag_2d_dflux_pv(:,:,:,:,1)+diag_2d_dflux_pv(:,:,:,:,ith)
       enddo
       endif

       call my_mpi_reduce( diag_2d_dflux_pv(:,:,:,:,1),dum1,2*ne*np*(ptl_nsp-ptl_isp+1))
       diag_2d_dflux_pv(:,:,:,:,1)=dum1
     endif
  endif

  if(sml_mype==0) then
     do isp=ptl_isp, ptl_nsp   ! species loop
        out1(1,:,isp)=diag_1d_f_pv1(1,:,isp,1)/diag_1d_vol ! density
        do ip=1,np
           out1(2:nv1,ip,isp)=diag_1d_f_pv1(2:nv1,ip,isp,1)/diag_1d_f_pv1(1,ip,isp,1)*unit1(2:nv1)
        enddo
        ! delta-f calc. : d(nL) = L*dn+n*dL -> dL = (d(nL)-L*dn)/n where L is a diagnosed quantity
        if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
           if(.not. sml_f0_grid) then
              ! this part is designed to get variation from f0, but not working well
              out1(1+nv1,:,isp)=diag_1d_df_pv1(1,:,isp,1)/diag_1d_vol ! delta density
              do ip=1,np
                 out1(2+nv1:nv1*2,ip,isp)=&
                      (diag_1d_df_pv1(2:nv1,ip,isp,1)) * unit1(2:nv1)
              enddo
           else
              ! add f0 component
              diag_1d_df_pv1(:,:,isp,1)=diag_1d_df_pv1(:,:,isp,1) + diag_1d_f0_pv1(:,:,isp,1)
              !
              out1(1+nv1,:,isp)=(diag_1d_df_pv1(1,:,isp,1))/diag_1d_vol ! density
              do ip=1,np
!                 out1(2+nv1:nv1*2,ip,isp)= diag_1d_df_pv1(2+nv1:nv1*2,ip,isp,1)/diag_1d_df_pv1(1+nv1,ip,isp,1)*unit1(2:nv1)
                 out1(2+nv1:nv1*2,ip,isp)= diag_1d_df_pv1(2:nv1,ip,isp,1)/diag_1d_df_pv1(1,ip,isp,1)*unit1(2:nv1)
              enddo
           endif
        endif
     enddo
!    print *, 'diag_1d_f_pv1=',diag_1d_f_pv1

     !**** time average output ****
     if(diag_tavg_on) then
        do isp=ptl_isp, ptl_nsp   ! species loop
           out1_tavg(1,:,isp)=diag_1d_tavg_f_pv1(1,:,isp,1)/diag_1d_vol/real(diag_1d_period) ! density
           do ip=1,np
              out1_tavg(2:nv1,ip,isp)=diag_1d_tavg_f_pv1(2:nv1,ip,isp,1)/diag_1d_tavg_f_pv1(1,ip,isp,1)*unit1(2:nv1)
           enddo

           ! delta-f calc. : d(nL) = L*dn+n*dL -> dL = (d(nL)-L*dn)/n where L is a diagnosed quantity
           if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
              out1_tavg(1+nv1,:,isp)=diag_1d_tavg_df_pv1(1,:,isp,1)/diag_1d_vol ! delta density
              do ip=1,np
                 out1_tavg(2+nv1:nv1*2,ip,isp)=&
                      (diag_1d_tavg_df_pv1(2:nv1,ip,isp,1)-out1_tavg(2:nv1,ip,isp)*diag_1d_tavg_df_pv1(1,ip,isp,1)) / &
                      diag_1d_tavg_f_pv1(1,ip,isp,1)*unit1(2:nv1)
              enddo
           endif
        enddo
     endif
     !**** time average omid output ****
     if(diag_tavg_on .and. diag_omid_on) then
        do isp=ptl_isp, ptl_nsp   ! species loop
           out1_omid(1,:,isp)=diag_1d_omid_f_pv1(1,:,isp,1)/diag_1d_vol/real(diag_1d_period) ! density
           do ip=1,np
              out1_omid(2:nv1,ip,isp)=diag_1d_omid_f_pv1(2:nv1,ip,isp,1)/diag_1d_omid_f_pv1(1,ip,isp,1)*unit1(2:nv1)
           enddo

           ! delta-f calc. : d(nL) = L*dn+n*dL -> dL = (d(nL)-L*dn)/n where L is a diagnosed quantity
!           if(sml_deltaf) then
!              out1_tavg(1+nv1,:,isp)=diag_1d_tavg_df_pv1(1,:,isp,1)/diag_1d_vol ! delta density
!              do ip=1,np
!                 out1_tavg(2+nv1:nv1*2,ip,isp)=&
!                      (diag_1d_tavg_df_pv1(2:nv1,ip,isp,1)-out1_tavg(2:nv1,ip,isp)*diag_1d_tavg_df_pv1(1,ip,isp,1)) / &
!                      diag_1d_tavg_f_pv1(1,ip,isp,1)*unit1(2:nv1)
!              enddo
!           endif
        enddo
     endif

   !**** kinetic ion eflux
    if(diag_eflux_on) then
       out1_eflux=0D0
       out1_weight=0D0
       out1_pflux=0D0
       out1_test=0D0
       do isp=ptl_isp, ptl_nsp
          do ip=1,np
             do ie=1,ne
                if(abs(diag_1d_eflux_pv(1,ip,ie,isp,1))>1D-20) then
                out1_weight(ip,ie,isp)=diag_1d_eflux_pv(1,ip,ie,isp,1)
                out1_eflux(ip,ie,isp)=diag_1d_eflux_pv(2,ip,ie,isp,1)/out1_weight(ip,ie,isp)
                out1_pflux(ip,ie,isp)=out1_eflux(ip,ie,isp)/(diag_1d_emin+ie*diag_1d_de)
                endif
                out1_test(ip,isp)=out1_test(ip,isp)+out1_eflux(ip,ie,isp)*(diag_1d_emin+ie*diag_1d_de)
             enddo
          enddo
       enddo
       out1_test = out1_test/(diag_1d_emax-diag_1d_emin)
       if(sml_deltaf .or. (sml_deltaf_elec .and. sml_electron_on)) then
         if(.not. sml_f0_grid) then
           print *, "2d particle and heat flux is not ready yet for deltaf-f method"
         else
         ! add f0 component
         diag_2d_dflux_pv(:,:,:,:,1)=diag_2d_dflux_pv(:,:,:,:,1)+diag_2d_ef0_pv(:,:,:,:,1)
         out1_eflux=0D0
         out1_weight=0D0
         out1_pflux=0D0
         out1_test=0D0
         tmp1=0D0
         tmp2=0D0
         do isp=ptl_isp, ptl_nsp
            do ip=1,np
               do ie=1,ne
                  if(abs(diag_2d_dflux_pv(1,ip,ie,isp,1))>1D-10) then
                  out1_weight(ip,ie,isp)=diag_2d_dflux_pv(1,ip,ie,isp,1)
                  out1_eflux(ip,ie,isp)=diag_2d_dflux_pv(2,ip,ie,isp,1)
                  out1_pflux(ip,ie,isp)=out1_eflux(ip,ie,isp)
                  endif
                  tmp1=tmp1+diag_2d_dflux_pv(1,ip,ie,isp,1)
                  tmp2=tmp2+diag_2d_dflux_pv(2,ip,ie,isp,1)
               enddo
               out1_test(ip,isp)=tmp2/tmp1
            enddo
         enddo
         endif ! sml_f0_grid
       endif ! sml_deltaf           
    endif ! diag_eflux_on

  endif !sml_mype

#ifdef ADIOS
  if(sml_mype==0) then
     buf_size=4 + 1000 + 8*np*(ptl_nsp-ptl_isp+1)*(nv_total)*2 + 8*grid%npsi00*(5)
     if(diag_tavg_on) buf_size = buf_size + 8*np*(ptl_nsp-ptl_isp+1)*(nv_total)*2
     if(diag_tavg_on .and. diag_omid_on )  &
                 buf_size = buf_size + 8*np*(ptl_nsp-ptl_isp+1)*(nv_total)*2
     if(diag_eflux_on) buf_size = buf_size + 8*np*ne*(ptl_nsp-ptl_isp+1)*3 +8*np*(ptl_nsp-ptl_isp+1)

     do ip=1,np
        pall(ip)=diag_1d_pin+diag_1d_dp*real(ip-1)
     enddo
     pnorm=pall/eq_x_psi

     !print *, 'adios writing diagnosis.1d group, #bytes = ', buf_size, etag
#ifdef SC17DEMO
     if(sml_gstep/diag_1d_period==1) then
        ADIOS_OPEN(buf_id,'diagnosis.1d','xgc.oneddiag.bp','w',MPI_COMM_SELF,err)
     else
        ADIOS_OPEN(buf_id,'diagnosis.1d','xgc.oneddiag.bp','a',MPI_COMM_SELF,err)
     endif
#else
     if(sml_gstep/diag_1d_period==1) then
        ADIOS_OPEN(buf_id,'diagnosis.1d','xgc.oneddiag.bp','w',sml_comm_null,err)
     else
        ADIOS_OPEN(buf_id,'diagnosis.1d','xgc.oneddiag.bp','a',sml_comm_null,err)
     endif
#endif
     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
     ADIOS_WRITE_LBL(buf_id,'samples',np,err)
     ADIOS_WRITE_LBL(buf_id,'gsamples',grid%npsi00,err)
     ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
     ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
     ADIOS_WRITE_LBL(buf_id,'tindex',sml_gstep/diag_1d_period,err)
     ADIOS_WRITE_LBL(buf_id,'psi_mks',pall,err)
     ADIOS_WRITE_LBL(buf_id,'psi',pnorm,err)
     ADIOS_WRITE_LBL(buf_id,'nenergy',diag_1d_ne,err)

     do isp=ptl_isp, ptl_nsp
        do j=1,nv_total
           tmp(:,j,isp)=out1(j,:,isp)
           ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(j))//'_1d',tmp(:,j,isp),err)
        enddo
! for debug, ASCII output together with BP----------
!       write(700+isp,*) 'isp = ',isp
!       do ip=1,np
!         write(700+isp,*) 'ip = ',ip
!         write(700+isp,1000) pnorm(ip),out1(:,ip,isp)
!       enddo
 1000   format(4(1pe19.6,1x))

     enddo
     !time average
     if(diag_tavg_on) then
        do isp=ptl_isp, ptl_nsp
           do j=1,nv_total
              tmp_tavg(:,j,isp)=out1_tavg(j,:,isp)
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(j))//'_avg',tmp_tavg(:,j,isp),err)
           enddo
        enddo
     endif

     !omid
     if(diag_tavg_on .and. diag_omid_on) then
        do isp=ptl_isp, ptl_nsp
           do j=1,nv_total
              tmp_omid(:,j,isp)=out1_omid(j,:,isp)
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(j))//'_omid',tmp_omid(:,j,isp),err)
           enddo
        enddo
     endif

     !kinetic ion flux
     if(diag_eflux_on) then
        do isp=ptl_isp, ptl_nsp
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim('radial_eflux'),out1_eflux(:,:,isp),err)
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim('radial_pflux'),out1_pflux(:,:,isp),err)
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim('radial_weight'),out1_weight(:,:,isp),err)
              ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim('radial_test'),out1_test(:,isp),err)
        enddo
     endif

     !do j=1,diag_1d_npv2
     !   tmp2(:,j)=diag_1d_pv2(j,:)
     !   ADIOS_WRITE_LBL(buf_id,trim(var_names2(j))//postfix,tmp2(:,j),err)
     !enddo

     do j=1, grid%npsi00
        psi00(j)=grid%psi00min+real(j-1)*grid%dpsi00
     enddo

     ADIOS_WRITE_LBL(buf_id,'psi00', psi00,err)
     ADIOS_WRITE_LBL(buf_id,'iden00_1d', psn%iden00_1d,err)
     ADIOS_WRITE_LBL(buf_id,'eden00_1d', psn%eden00_1d,err)
     ADIOS_WRITE_LBL(buf_id,'cden00_1d', psn%cden00_1d,err)
     ADIOS_WRITE_LBL(buf_id,'pot00_1d', psn%pot00_1d,err)
     !ADIOS_WRITE_LBL(buf_id,'dpotdp'    //postfix, psn%dpotdp    ,err)

     ADIOS_CLOSE(buf_id,err)
  endif
  deallocate(tmp,var_names)
  if(diag_tavg_on) deallocate(tmp_tavg)
  if(diag_tavg_on .and. diag_omid_on) deallocate(tmp_omid)

#endif

  deallocate(out1)
  if(diag_tavg_on) deallocate(out1_tavg)
  if(diag_tavg_on .and. diag_omid_on) deallocate(out1_omid)
  if(diag_eflux_on) then
     deallocate(out1_eflux)
     deallocate(out1_pflux)
     deallocate(out1_weight)
     deallocate(out1_test)
  endif
end subroutine diag_1d_output
!#ifdef DIAG_1D_F0
subroutine diag_1d_f0(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use f0_module
  use diag_module
  use ptl_module
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z

  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: isp, imu, node, ivp
  real (8) :: mu, smu, vp, en, pot, f0, f0_p
  real (8) :: vol_sum(diag_1d_npsi,ptl_isp:ptl_nsp,sml_nthreads)
  real (8) :: dum(diag_1d_npv1,diag_1d_npsi,ptl_isp:ptl_nsp), dum2(diag_1d_npsi,ptl_isp:ptl_nsp)
  real (8) :: dum1(2,diag_1d_npsi,diag_1d_ne,ptl_isp:ptl_nsp)
  real (8) :: v(diag_1d_npv1)
  real (8) :: pvol, pn
  integer :: ip, ie, ierr
  integer :: ith
  real (8) :: mu_vol, vol,vol_rspace, vth, wp, en_th, den, we, pe, diag_1d_de, w, en_p
  ! For the drift velocities on the v-grid
  real (8) :: v_mag(3), v_exb(3), v_pardrift(3), grad_psi_sqr
  real (8) :: psi, r, z

  !rh Initalize with a small value to
  !rh avoid division by zero
  !rh since diag_1d_f0_pv1 will be zero anyway where there
  !rh are no flux surfaces for the uniform psi grid, this
  !rh does not produce large values upon division
  !rh vol_sum=0D0
  vol_sum=1D-15
  diag_1d_f0_pv1=0D0
  diag_2d_ef0_pv=0D0

  ith=1 !open mp thread

  do isp=diag_1d_isp, diag_1d_nsp
     do imu=f0_imu1, f0_imu2

        if(imu==0) then
           mu_vol=0.5D0
        elseif(imu==f0_nmu) then
           mu_vol=0.5D0
        else
           mu_vol=1D0
        endif

        do node=f0_inode1,f0_inode2
           ! volume for data point
           vol=f0_grid_vol(node,isp)*mu_vol
           vol_rspace=grid%node_vol_nearest(node)
           ! normalized v
           en_th=f0_t_ev(node,isp)*sml_ev2j
           vth=sqrt(en_th/ptl_mass(isp))


           !get psi index
           pn=(grid%psi(node)-diag_1d_pin)*diag_1d_dp_inv
           ip=floor(pn)+1
           if(ip <1 .or. diag_1d_npsi <= ip) cycle
           wp=1D0 - pn + real(ip-1,8)


           !vol_sum is summation of real space volume
           !this will be used for converting all values according to diag_1d_vol
           !if(imu==0) then
           if(imu==f0_imu1) then
              vol_sum(ip  ,isp,ith)=vol_sum(ip  ,isp,ith) + vol_rspace*wp
              vol_sum(ip+1,isp,ith)=vol_sum(ip+1,isp,ith) + vol_rspace*(1D0-wp)
           endif


           call get_derivs_at_grid  !local subroutine


           !adding density and temperature of f0 separately
           if(imu==0) then
              den=f0_den(node)

              diag_1d_f0_pv1(1,ip  ,isp,ith)=diag_1d_f0_pv1(1,ip  ,isp,ith) + den*vol_rspace*wp
              diag_1d_f0_pv1(7,ip  ,isp,ith)=diag_1d_f0_pv1(7,ip  ,isp,ith) + den*0.5D0*en_th*vol_rspace*wp
              diag_1d_f0_pv1(8,ip  ,isp,ith)=diag_1d_f0_pv1(8,ip  ,isp,ith) + den*en_th*vol_rspace*wp


              diag_1d_f0_pv1(1,ip+1,isp,ith)=diag_1d_f0_pv1(1,ip+1,isp,ith) + den*vol_rspace*(1D0-wp)
              diag_1d_f0_pv1(7,ip+1,isp,ith)=diag_1d_f0_pv1(7,ip+1,isp,ith) + den*0.5D0*en_th*vol_rspace*(1D0-wp)
              diag_1d_f0_pv1(8,ip+1,isp,ith)=diag_1d_f0_pv1(8,ip+1,isp,ith) + den*en_th*vol_rspace*(1D0-wp)
           endif


           ! r, b, bp, etc - later
           if (imu==0) then
              smu=f0_dsmu/f0_mu0_factor
           else
              smu=imu*f0_dsmu  ! <--- This is v_perp/v_th for v_perp grid!!!!
           endif
           mu=imu*f0_dsmu
           mu=mu*mu  ! <--- This is (v_perp/v_th)^2 for v_perp grid and mu_N for sqrt(mu) grid!!!!

           if(isp==0) pot= 0.5D0 * (psn%dpot_ff(node,0)+psn%dpot_ff(node,1))/f0_t_ev(node,isp)
           do ivp=-f0_nvp, f0_nvp

              !1. get f0 = f0a + f0g
              vp=ivp*f0_dvp

              !2. Get drift velocities
              call get_drift_velocity(grid,node,mu,vp,isp,vth,v_mag,v_exb,v_pardrift,grad_psi_sqr)

              en=0.5d0 * (vp*vp + mu) ! energy normalized by T
              if(isp==0) en = en - pot
              ! No background Maxwellian in the other diagnostics -->
#ifndef F0_TOR_LINEAR
              f0=f0_f0g(ivp,node,imu,isp)
#else
              f0=(f0_f0g(ivp,node,imu,0,isp)+f0_f0g(ivp,node,imu,1,isp))*0.5D0
#endif
              f0_p=f0 + f0_n_Ta(node,isp)*exp(-en)*smu

              !3. get per particle value
              !rh added contribution of f0_f0g to 1D diagnostics, only density (1), temperature (7 and 8) contain f0_ana, too.
              v(1)  = 1D0
              v(2)  = (v_mag(3)+v_exb(3)+v_pardrift(3)) !0D0 ! toroidal velocity -- NOT correct
              v(3)  = (v_mag(2)+v_exb(2)+v_pardrift(2)) !0D0 ! poloidal velocity -- NOT correct
              v(4)  = vp*vth ! parallel velocity MKS
              v(5)  = vp * grid%x(1,node) * grid%bfield(3,node)/grid%bfield(4,node) !0D0 ! tor. ang. momentum ! Not correct
              v(6)  = (v_mag(1)+v_exb(1)+v_pardrift(1)) !0D0 ! radial drift
              v(7)  = 0.5D0*vp*vp*en_th ! par. mean. en.
              v(8)  = 0.5D0*mu*en_th ! perp. en.
              v(9)  = (v(7)+v(8)) * v(6) !0D0 ! radial en. flux
              v(10) = v(4) * v_exb(1) * grid%bfield(3,node)/grid%bfield(4,node) !0D0 ! rad. momentum flux due to ExB
              v(11) = v(4) * v(6) * grid%bfield(3,node)/grid%bfield(4,node) !0D0 ! rad. momentum flux (total)
              v(12) = v_exb(1) !0D0 ! radial ExB drift
              v(13) = (v(7)+v(8)) * v_exb(1) !0D0 ! radial energy flux due to ExB
              v(14) = v_exb(2) !0D0 ! poloidal ExB flow
              v(15) = grad_psi_sqr !0D0 ! |grad(Psi)|^2
              diag_1d_f0_pv1(:,ip  ,isp,ith)=diag_1d_f0_pv1(:,ip  ,isp,ith) + v(:)*f0*vol*wp
              diag_1d_f0_pv1(:,ip+1,isp,ith)=diag_1d_f0_pv1(:,ip+1,isp,ith) + v(:)*f0*vol*(1D0-wp)

              if(diag_eflux_on)then
                diag_1d_emin = 0D0
                diag_1d_emax = 3D0*eq_ftn(sml_inpsi,eq_axis_r,eq_axis_z,eq_tempi)*sml_ev2j
                diag_1d_de=(diag_1d_emax-diag_1d_emin)/diag_1d_ne

                en_p=v(7)+v(8)
                pe=(en_p-diag_1d_emin)/diag_1d_de
                ie=floor(pe)+1

                psi = grid%psi(node)
                r=grid%x(1,node)
                z=grid%x(2,node)
                if(.not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z > eq_x2_z) cycle ! EXIT for private region of lower divertor
                if(ie >=1 .and. diag_1d_ne > ie) then
                  we=1D0 - pe + real(ie-1,8)
                  diag_2d_ef0_pv(1,ip  ,ie  ,isp,ith)=f0_p*wp*we*vol            +diag_2d_ef0_pv(1,ip  ,ie  ,isp,ith)
                  diag_2d_ef0_pv(1,ip+1,ie  ,isp,ith)=f0_p*(1D0-wp)*we*vol      +diag_2d_ef0_pv(1,ip+1,ie  ,isp,ith)
                  diag_2d_ef0_pv(1,ip  ,ie+1,isp,ith)=f0_p*wp*(1D0-we)*vol      +diag_2d_ef0_pv(1,ip  ,ie+1,isp,ith)
                  diag_2d_ef0_pv(1,ip+1,ie+1,isp,ith)=f0_p*(1D0-wp)*(1D0-we)*vol+diag_2d_ef0_pv(1,ip+1,ie+1,isp,ith)

                  diag_2d_ef0_pv(2,ip  ,ie  ,isp,ith)=v(9)*f0*wp*we*vol            +diag_2d_ef0_pv(2,ip  ,ie  ,isp,ith)
                  diag_2d_ef0_pv(2,ip+1,ie  ,isp,ith)=v(9)*f0*(1D0-wp)*we*vol      +diag_2d_ef0_pv(2,ip+1,ie  ,isp,ith)
                  diag_2d_ef0_pv(2,ip  ,ie+1,isp,ith)=v(9)*f0*wp*(1D0-we)*vol      +diag_2d_ef0_pv(2,ip  ,ie+1,isp,ith)
                  diag_2d_ef0_pv(2,ip+1,ie+1,isp,ith)=v(9)*f0*(1D0-wp)*(1D0-we)*vol+diag_2d_ef0_pv(2,ip+1,ie+1,isp,ith)
                endif
              endif
           enddo
        enddo
     enddo
  enddo

  if(sml_nthreads >= 2) then
  do ith=2, sml_nthreads
     diag_1d_f0_pv1(:,:,:,1) = diag_1d_f0_pv1(:,:,:,1)+diag_1d_f0_pv1(:,:,:,ith)
     if(diag_eflux_on) then
       diag_2d_ef0_pv(:,:,:,:,1) = diag_2d_ef0_pv(:,:,:,:,1)+diag_2d_ef0_pv(:,:,:,:,ith)
     endif
     vol_sum(:,:,1)=vol_sum(:,:,1)+vol_sum(:,:,ith)
  enddo
  endif

  !mpi reduce
  call my_mpi_reduce(diag_1d_f0_pv1(:,:,:,1),dum,diag_1d_npv1*diag_1d_npsi*(diag_1d_nsp-diag_1d_isp+1))
  diag_1d_f0_pv1(:,:,:,1)=dum
  call my_mpi_reduce(vol_sum(:,:,1),dum2,diag_1d_npsi*(diag_1d_nsp-diag_1d_isp+1))
  vol_sum(:,:,1)=dum2
  if(diag_eflux_on) then
    dum1=0D0
    call mpi_reduce(diag_2d_ef0_pv(:,:,:,:,1),dum1,2*diag_1d_npsi*diag_1d_ne*(diag_1d_nsp-diag_1d_isp+1), &
         MPI_REAL8, MPI_SUM, 0, sml_COMM,ierr)
    diag_2d_ef0_pv(:,:,:,:,1)=dum1
  endif

  ! correction with diag_1d_vol/vol_sum
  if(sml_mype==0) then
     do isp=diag_1d_isp, diag_1d_nsp
        do ip=1, diag_1d_npsi
           diag_1d_f0_pv1(:,ip,isp,1)=diag_1d_f0_pv1(:,ip,isp,1)*diag_1d_vol(ip)/vol_sum(ip,isp,1)
           if(diag_eflux_on) diag_2d_ef0_pv(:,ip,:,isp,1)=diag_2d_ef0_pv(:,ip,:,isp,1)*diag_1d_vol(ip)/vol_sum(ip,isp,1)
        enddo
     enddo
  endif

contains
  subroutine get_derivs_at_grid
    implicit none





  end subroutine get_derivs_at_grid

  ! Evaluates the drift velocities at given inode, imu, ivp
  ! The basic magnetic field quantities necessary for the calculation of the
  ! magnetic drifts are stored in the grid data structure.
  ! To save memory, this could also be shifted to f0_module.
  subroutine get_drift_velocity(grid,node,mu,vp,isp,v_th,v_mag,v_exb,v_pardrift,grad_psi_sqr)
    use eq_module, only: eq_axis_b
    use ptl_module, only: ptl_mass, ptl_charge
    use f0_module, only: f0_T_ev, f0_dvp, f0_dsmu
    use grid_class
    implicit none
    ! Input/Output
    type(grid_type), intent(in) :: grid
    integer, intent(in) :: node, isp
    real (kind=8), intent(in) :: v_th, mu, vp
    real (kind=8), intent(out), dimension(3) :: v_mag, v_exb, v_pardrift
    real (kind=8), intent(out) :: grad_psi_sqr
    ! Internal variables
    real (kind=8) :: D, rho_mks, mu_mks, E(3), wrho(2), rhoi, rhon, vr, vz
    real (kind=8) :: b, over_B2, cmrho2, cmrho, murho2b, murho2b_c, grad_psi(2)
    real (kind=8) :: mass, charge, over_abs_grad_psi
    integer :: irho
    real (kind=8), external :: psi_interpol

    E=0.D0
    charge=ptl_charge(isp)
    mass=ptl_mass(isp)
    b=grid%bfield(4,node)
    rho_mks=vp*v_th/b*mass/charge
    mu_mks=0.5*mu*f0_T_ev(node,isp)*sml_ev2j/b

    over_B2=1.D0/(b**2)

    cmrho2=(charge/mass)*rho_mks**2
    cmrho =(charge/mass)*rho_mks
    murho2b=(mu_mks+charge*cmrho2 *B)
    murho2b_c=murho2b/charge

    D=1.D0/ ( 1.D0 + rho_mks * grid%nb_curl_nb(node) )

    grad_psi(1)=psi_interpol(grid%x(1,node),grid%x(2,node),1,0)
    grad_psi(2)=psi_interpol(grid%x(1,node),grid%x(2,node),0,1)
    grad_psi_sqr=sum(grad_psi(:)**2)
    over_abs_grad_psi=1.D0/sqrt(grad_psi_sqr)

    ! The magnetic drifts (curvature + grad_b)
    v_mag(:) = D * ( grid%v_gradb(:,node)*murho2b_c*over_B2 &
                    +grid%v_curv(:,node)*cmrho2 )


!    if (isp==1) then
    if (.false.) then
      ! Ions
      rhoi=sqrt(2.D0*mass*mu_mks/grid%bfield(4,node))/charge
      rhon=min(rhoi,grid%rhomax)/grid%drho
      irho=min(floor(rhon),grid%nrho-1)
      wrho(2)=rhon - real(irho)
      wrho(1)=1D0-wrho(2)

      E   = E   + wp*wrho(1)*0.5D0*(psn%E_rho_ff(:,0,irho  ,node)+psn%E_rho_ff(:,1,irho  ,node))
      E   = E   + wp*wrho(2)*0.5D0*(psn%E_rho_ff(:,0,irho+1,node)+psn%E_rho_ff(:,1,irho+1,node))
    else
      ! Electrons
      E = 0.5D0*(psn%E_rho_ff(:,0,0,node)+psn%E_rho_ff(:,1,0,node))
    endif

    ! The ExB drift
    vr       = D * ( E(3)*grid%bfield(2,node) - E(2)*grid%bfield(3,node) ) * over_B2
    vz       = D * ( E(1)*grid%bfield(3,node) - E(3)*grid%bfield(1,node) ) * over_B2
    v_exb(1) = vr * grad_psi(1) + vz * grad_psi(2)
    v_exb(2) = (- vr * grad_psi(2) + vz * grad_psi(1)) * over_abs_grad_psi
    v_exb(3) = D * ( E(2)*grid%bfield(1,node) - E(1)*grid%bfield(2,node) ) * over_B2

    ! The parallel drift
    vr = D * cmrho * grid%bfield(1,node)
    vz = D * cmrho * grid%bfield(2,node)
    v_pardrift(1) = vr * grad_psi(1) + vz * grad_psi(2)
    v_pardrift(2) = (- vr * grad_psi(2) + vz * grad_psi(1)) * over_abs_grad_psi
    v_pardrift(3) = D * cmrho * grid%bfield(3,node)


  end subroutine get_drift_velocity

end subroutine diag_1d_f0
!#endif

!****************** 3D output *************************
subroutine diag_3d(istep,grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use eq_module
  use ptl_module, only: ptl_mass
#ifdef SC17DEMO
  !! jyc: temporary fix for SC17 demo
!  use coupling_core_edge
use new_coupling
#endif
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,istep,iphi
  character (len=30) :: filename, filename2
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  integer, parameter :: nspace=2
  character (len=30) :: pathname
  character (len=40) :: labelname
  real (kind=8) , allocatable :: potential(:), pot3d(:,:), tmp(:,:), iden_f0_approx(:)
#endif
#ifdef XGC1_EM
  real (kind=8), allocatable :: grad_dpot(:,:), grad_apar(:,:), en_den(:), en_dpot(:)
  real (kind=8), allocatable :: en_apar(:), b2(:), den(:), teev(:), eta_loss(:), tearing_drive(:)
  real (kind=8), allocatable :: tearing_drive2(:), gradphi_apar(:), eps0_chi(:), dum(:), tearing_drive3(:)
#endif
  integer   :: nnode,nphi
  integer   :: signal, ierr
  logical, parameter :: handshake = .true.
#ifdef SC17DEMO
  character (len=30) :: filename3, filename4
#endif
  nnode=grid%nnode
  nphi=grid%nphi

1000 format(19(e19.13,' '))
#if defined(ADIOS)
  !if (.false.) then
#ifdef FIELDP
  if (mod(istep,diag_3d_period)==0) then
     ! collect and write out the entire 3d turbulent electrostatic potential
     write(filename2,'("xgc.fieldp",".",i4.4,".bp")') (sml_gstep)/(diag_3d_period)
     if (sml_mype==0) then
        allocate(potential(nnode))
        potential(:) = psn%pot0(:)  ! normalized potential (Volts)
        buf_size = 4 + 4 + 4 + nnode * 8
        !        print *, 'adios writing pot0 group, #bytes = ', buf_size
        ADIOS_OPEN(buf_id,'pot0',filename2,'w',sml_comm_null,err)
        ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
        ADIOS_WRITE_LBL(buf_id,'istep',istep,err)
        ADIOS_WRITE_LBL(buf_id,'nphiP1',sml_nphi_total+1,err)
        ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
        ADIOS_WRITE_LBL(buf_id,'/node_data[0]/values',potential,err)
        ADIOS_CLOSE(buf_id,err)
        deallocate(potential)
     endif
     ! turbulence data only
     if (sml_turb_poisson) then
        if (sml_plane_mype==0) then
           nphi=1
           allocate(pot3d(nnode,1))
           pot3d(:,1) = psn%dpot(:,1)

           buf_size = 4 + 40 + nnode * 8 * 4
	   if (sml_mype==0) then
       !             print *, 'adios writing pot3d group, #bytes = ', buf_size
              ADIOS_OPEN(buf_id,'pot3d',filename2,'a',sml_adios_comm,err)
           else
              ADIOS_OPEN(buf_id,'pot3d',filename2,'w',sml_adios_comm,err)
           endif
           do iphi=1,nphi
              write(pathname,'("/node_data[",i0,"]")') (sml_mype/sml_pe_per_plane)*sml_plane_per_pe+iphi
              write(labelname,'("potential.",i3.3)') (sml_mype/sml_pe_per_plane)*sml_plane_per_pe+iphi-1
              labelname = trim(labelname)//char(0)
              ADIOS_SET_PATH(buf_id,pathname,err)
              ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
              ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
              ADIOS_WRITE_LBL(buf_id,'labelname',labelname,err)
              ADIOS_WRITE_LBL(buf_id,'values',pot3d(:,iphi),err)
           enddo
           ADIOS_CLOSE(buf_id,err)
           deallocate(pot3d)
        endif
     endif
  endif
#endif
  !##################################################################
  ! new file format -- 3d global array


  if (mod(istep,diag_3d_period)==0) then
     ! collect and write out the entire 3d turbulent electrostatic potential
     if (adios_stage_field3d) then
       write(filename2,'("xgc.3d.bp")')
#ifdef SC17DEMO
       !! jyc: temporary fix for SC17 demo
       write(filename2,'("xgc.3d.", A4, ".bp")') cce_my_side
#endif
     else
       write(filename2,'("xgc.3d",".",i5.5,".bp")') (sml_gstep)/(diag_3d_period)
     endif
     ! turbulence data only
     if (sml_turb_poisson) then
#ifdef XGC1_EM
       if (sml_electron_hyb) then
         !rh Simple energy diagnostic
         !rh Only for local delta-f case!

         allocate(grad_dpot(grid%nnode,2),grad_apar(grid%nnode,2))
         allocate(en_den(grid%nnode),en_dpot(grid%nnode),en_apar(grid%nnode))
         allocate(eta_loss(grid%nnode),tearing_drive(grid%nnode))
         allocate(b2(grid%nnode),den(grid%nnode),teev(grid%nnode),eps0_chi(grid%nnode))
         allocate(gradphi_apar(grid%nnode),tearing_drive2(grid%nnode))
         allocate(dum(grid%nnode),tearing_drive3(grid%nnode))

         call grid_deriv(grid,psn%dpot(:,1),grad_dpot(:,1),grad_dpot(:,2))
         call grid_deriv(grid,psn%a_par(:),grad_apar(:,1),grad_apar(:,2))
         call GradPhiX(grid,psn%a_par(:),gradphi_apar(:))

         do i=1,grid%nnode
           b2(i)=grid%bfield(4,i)**2
           den(i)=eq_ftn(grid%psi(i),grid%x(1,i),grid%x(2,i),eq_den)
           teev(i)=eq_ftn(grid%psi(i),grid%x(1,i),grid%x(2,i),eq_tempe)
           eps0_chi(i)=ptl_mass(1)*den(i)/b2(i)
         enddo

         en_den(:)=sml_e_charge*teev(:)/den(:)*psn%eden_hyb(:)**2
         en_dpot(:)=eps0_chi(:)*sum(grad_dpot(:,:)**2,2)
         en_apar(:)=sum(grad_apar(:,:)**2,2)/sml_mu0

         eta_loss(:) = sml_eta*(sml_e_charge*den(:)*psn%u_e(:))**2
         tearing_drive(:) = sml_e_charge*(teev/den*psn%eden_hyb-psn%dpot(:,1)) &
                         * (grid%tearing_drive(:,1)*grad_apar(:,1) + grid%tearing_drive(:,2)*grad_apar(:,2))
         tearing_drive2(:) = sml_e_charge*(teev/den*psn%eden_hyb-psn%dpot(:,1)) * grid%tearing_drive(:,3)*gradphi_apar(:)
         tearing_drive3(:) = sml_e_charge*(teev/den*psn%eden_hyb-psn%dpot(:,1)) &
                            *grid%tearing_drive2*psn%a_par

         deallocate(grad_dpot,grad_apar,gradphi_apar,den,teev,b2,eps0_chi)

         dum(:)=0D0
         call mpi_reduce(en_den(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         en_den(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(en_apar(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         en_apar(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(en_dpot(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         en_dpot(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(eta_loss(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         eta_loss(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(tearing_drive(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         tearing_drive(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(tearing_drive2(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         tearing_drive2(:)=dum(:)
         dum(:)=0D0
         call mpi_reduce(tearing_drive3(:),dum(:),grid%nnode,MPI_REAL8,MPI_SUM,0,sml_intpl_comm,ierr)
         tearing_drive3(:)=dum(:)

         deallocate(dum)
       endif
#endif

        if (sml_plane_mype==0) then
           buf_size = 1000 + nnode * 8 * 2
#ifdef DDEN_DIAG
           buf_size = buf_size + nnode*8
#endif
#ifdef XGC1_EM
           if (sml_electron_hyb) then
             buf_size = buf_size + nnode * 8 * 7
           endif
#endif
           if(sml_electron_on) buf_size = buf_size + nnode*8*2
           if(sml_mype==0) then
              buf_size = buf_size + nnode*8  ! for pot0
              !if(sml_iter_solver .or. sml_poisson_solver_type /=0 ) then
                 buf_size = buf_size + nnode*8  ! for pot0m
              !endif
           endif
           if(diag_3d_more) then
                 buf_size = buf_size + nnode * 8 * diag_3d_nvar
                 if(sml_electron_on)   buf_size = buf_size + nnode * 8 * diag_3d_nvar
                 if(sml_electron_hyb) buf_size = buf_size + nnode*8*3
           endif
#ifdef RHS_DIAG
                 buf_size = buf_size + nnode * 8 * 2
#endif
           if(sml_electron_hyb) then  !why twice????
                 buf_size = buf_size + nnode * 8 * 3
           endif
#ifdef F0_CHARGE_N0
          if(sml_f0_grid) then
              buf_size = buf_size + 8*nnode
              if(sml_electron_on) buf_size = buf_size + 8*nnode
          endif
#endif

           !file open
           ADIOS_OPEN(buf_id,'field3D',filename2,'w',sml_intpl_comm,err)
           ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
           ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
           ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
           ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
           ADIOS_WRITE_LBL(buf_id,'dpot',psn%dpot(:,1),err)
           ADIOS_WRITE_LBL(buf_id,'iden',psn%idensity(:,1),err)
#ifdef DDEN_DIAG
           ADIOS_WRITE_LBL(buf_id,'dden',psn%dden(:,1),err)
#endif
#ifdef RHS_DIAG
           ADIOS_WRITE_LBL(buf_id,'rhs1',psn%rhs1,err)
           ADIOS_WRITE_LBL(buf_id,'rhs2',psn%rhs2,err)
#endif
           if(diag_3d_more) then
              allocate(tmp(grid%nnode,3))
              tmp(:,1)=diag_3d_add(1,:,1)
              ADIOS_WRITE_LBL(buf_id,'iden_real',tmp(:,1),err)
              tmp(:,2)=diag_3d_add(2,:,1)
              ADIOS_WRITE_LBL(buf_id,'itemp_real',tmp(:,2),err)
              tmp(:,3)=diag_3d_add(3,:,1)
              ADIOS_WRITE_LBL(buf_id,'iden_nogyro',tmp(:,3),err)
           endif
           !ADIOS_WRITE_LBL(buf_id,'dpot_real',psn%dpot_real(:,1),err)

           if(sml_electron_on) then
              ADIOS_WRITE_LBL(buf_id,'eden',psn%edensity(:,1),err)
              !ADIOS_WRITE_LBL(buf_id,'ddpotdt',psn%ddpotdt(:,1),err)  ! mpt meeded for usual case
           endif
#ifdef XGC1_EM
           if(sml_electron_hyb) then
              ADIOS_WRITE_LBL(buf_id,'eden_hyb',psn%eden_hyb(:),err)
              ADIOS_WRITE_LBL(buf_id,'a_par',psn%A_par,err)
              ADIOS_WRITE_LBL(buf_id,'u_e',psn%u_e,err)
              if (sml_mype==0) then
                ADIOS_WRITE_LBL(buf_id,'en_den',en_den,err)
                ADIOS_WRITE_LBL(buf_id,'en_dpot',en_dpot,err)
                ADIOS_WRITE_LBL(buf_id,'en_apar',en_apar,err)
                ADIOS_WRITE_LBL(buf_id,'eta_loss',eta_loss,err)
                ADIOS_WRITE_LBL(buf_id,'tearing_drive',tearing_drive,err)
                ADIOS_WRITE_LBL(buf_id,'tearing_drive2',tearing_drive2,err)
                ADIOS_WRITE_LBL(buf_id,'tearing_drive3',tearing_drive3,err)
              endif
           endif
#endif
           ! 0 potential
           if (sml_mype==0) then
              ADIOS_WRITE_LBL(buf_id,'pot0',psn%pot0,err)
              !if(sml_iter_solver .or. sml_poisson_solver_type /=0) then
                 ADIOS_WRITE_LBL(buf_id,'potm0',psn%pot0m,err)
              !endif
           endif


#ifdef F0_CHARGE_N0
           if(sml_f0_grid) then
              ADIOS_WRITE_LBL(buf_id,'iden_f0',psn%idensity_f0(:,1),err)
              if(sml_electron_on) ADIOS_WRITE_LBL(buf_id,'eden_f0',psn%edensity_f0(:,1),err)
           endif
#endif
           ADIOS_CLOSE(buf_id,err)


#ifdef SC17DEMO
          !! jyc: temporary fix for SC17 demo
          if (adios_stage_dpot) then
            write(filename3,'("xgc.dpot.", A4, ".bp")') cce_my_side
          else
            write(filename3,'("xgc.dpot",".",i5.5,".bp")') (sml_gstep)/(diag_3d_period)
          endif
          !file open
          ADIOS_OPEN(buf_id,'dpot',filename3,'w',sml_intpl_comm,err)
          ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
          ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
          ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
          ADIOS_WRITE_LBL(buf_id,'dpot',psn%dpot(:,1),err)
          ADIOS_CLOSE(buf_id,err)

          if (adios_stage_streamer) then
            write(filename4,'("xgc.streamer.", A4, ".bp")') cce_my_side
          else
            write(filename4,'("xgc.streamer",".",i5.5,".bp")') (sml_gstep)/(diag_3d_period)
          endif
          !file open
          ADIOS_OPEN(buf_id,'streamer',filename4,'w',sml_intpl_comm,err)
          ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
          ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
          ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
          ADIOS_WRITE_LBL(buf_id,'dpot',psn%dpot(:,1),err)
          ADIOS_CLOSE(buf_id,err)
#endif


           if(diag_3d_more) deallocate(tmp)
        endif
#ifdef XGC1_EM
        if (sml_electron_hyb) then
          deallocate(en_den,en_dpot,en_apar,eta_loss,tearing_drive,tearing_drive2,tearing_drive3)
        endif
#endif
     endif
  endif
#endif
end subroutine diag_3d

subroutine diag_f0(istep,grid,psn,flag)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use ptl_module, only : ptl_isp, ptl_nsp
  use f0_module
  implicit none
  integer, intent(in) :: istep
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer, intent(in), optional :: flag ! flag is for debug purpose
  !
  integer :: i
  integer :: nnode, inode1m1, ndata
  integer :: nphi, iphi
  integer :: nmup1, imu1m1, mudata
  integer :: vpdata
  character (len=30) :: filename
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: adios_comm, adios_comm_totalpe, adios_comm_mype
  integer :: err
  integer, parameter :: nspace=2
  character (len=30) :: pathname
  character (len=40) :: labelname
  real (8) :: iden_f0_approx(f0_inode1:f0_inode2), tmp(f0_inode1:f0_inode2)
  real (4), allocatable :: tmp_f_real4(:,:,:,:)
#ifndef F0_TOR_LINEAR
  real (4), allocatable :: tmp_g_real4(:,:,:,:)
#else
  real (4), allocatable :: tmp_g_real4(:,:,:,:,:)
#endif
#ifdef DIAG_F0_REAL4
  logical, parameter :: use_real4=.true.
  integer, parameter :: rsize=4
#else
  logical, parameter :: use_real4=.false.
  integer, parameter :: rsize=8
#endif
  include 'mpif.h'
#ifdef DIAG_NOISE
  interface
     subroutine diag_f0_noise(istep,grid,psn,flag)
       use grid_class
       use psn_class
       implicit none
       integer, intent(in) :: istep
       type(grid_type) :: grid
       type(psn_type) :: psn
       integer, intent(in), optional :: flag ! flag is for debug purpose
     end subroutine diag_f0_noise
  end interface
#endif
  if (mod(istep,diag_f0_period)==0) then
      nnode=grid%nnode
      inode1m1=f0_inode1-1
      ndata=f0_inode2-f0_inode1+1

      nphi=sml_nphi_total
      iphi=sml_intpl_mype

      nmup1= f0_nmu+1 ! include 0 index
      imu1m1=f0_imu1 -1 +1 ! additional +1 due to 0 starting index
      mudata=f0_imu2-f0_imu1+1

      vpdata=f0_nvp*2+1

      if (adios_mfilesys) then
         write(filename,'("restart_dir",i1.1,"/xgc.f0.",i5.5,".bp")') &
              & sml_mfilesys_index+1, (sml_gstep)/(diag_f0_period)
         if(present(flag)) then
            if(flag==-1) write(filename,'("restart_dir",i1.1,"/xgc.f0_debug.",i5.5,".bp")') &
                 & sml_mfilesys_index+1, (sml_gstep)/(diag_f0_period)
         endif
         adios_comm = sml_mfilesys_comm
         adios_comm_totalpe = sml_pe_per_filesys
         adios_comm_mype = sml_mfilesys_mype
      else
         if (adios_stage_f0) then
            write(filename,'("restart_dir/xgc.f0.bp")')
         else
            write(filename,'("restart_dir/xgc.f0.",i5.5,".bp")') (sml_gstep)/(diag_f0_period)
         endif
         if(present(flag)) then
            if(flag==-1) write(filename,'("restart_dir/xgc.f0_debug.",i5.5,".bp")') (sml_gstep)/(diag_f0_period)
         endif
         adios_comm = sml_comm
         adios_comm_totalpe = sml_totalpe
         adios_comm_mype = sml_mype
      endif

     if(sml_turb_poisson) then
#ifndef F0_TOR_LINEAR
        tmp=sum(psn%iden_rho_f0(f0_inode1:f0_inode2,:),2)
#else
        tmp=sum(psn%iden_rho_f0(f0_inode1:f0_inode2,1,:),2)
#endif
        iden_f0_approx=tmp  ! mpi reduce



        buf_size = 1000 + ndata*8            !iden_f0_approx
#ifndef F0_TOR_LINEAR
        buf_size = buf_size + ndata*mudata*vpdata*rsize  ! f0_g
#else
        buf_size = buf_size + ndata*mudata*vpdata*rsize*2  ! f0_g
#endif
        buf_size = buf_size + ndata*mudata*vpdata*rsize  ! f0_f
        buf_size = buf_size + ndata*8  !! for psn%dpot -- needs for ion only too?

        if(sml_electron_on) then
           buf_size = buf_size + ndata*8 !eden_f0
#ifndef F0_TOR_LINEAR
           buf_size = buf_size + ndata*mudata*vpdata*rsize  ! f0_g elec
#else
           buf_size = buf_size + ndata*mudata*vpdata*rsize*2  ! f0_g elec
#endif
           buf_size = buf_size + ndata*mudata*vpdata*rsize  ! f0_f elec
        endif

        !file open
        ADIOS_OPEN(buf_id,'f0',filename,'w',adios_comm,err)
        ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)

        ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
        ADIOS_WRITE_LBL(buf_id,'inode1m1',inode1m1, err)
        ADIOS_WRITE_LBL(buf_id,'ndata',ndata, err)

        ADIOS_WRITE_LBL(buf_id,'nphi',nphi,err)
        ADIOS_WRITE_LBL(buf_id,'iphi',iphi,err)


        ADIOS_WRITE_LBL(buf_id,'nmup1',nmup1, err)
        ADIOS_WRITE_LBL(buf_id,'imu1m1',imu1m1, err)
        ADIOS_WRITE_LBL(buf_id,'mudata',mudata, err)

        ADIOS_WRITE_LBL(buf_id,'vpdata',vpdata, err)

        !data
        ADIOS_WRITE_LBL(buf_id,'iden_f0_approx',iden_f0_approx,err)
#ifndef F0_TOR_LINEAR
        if(sml_electron_on) ADIOS_WRITE_LBL(buf_id,'eden_f0',psn%eden_f0(f0_inode1:f0_inode2),err)
#else
        if(sml_electron_on) ADIOS_WRITE_LBL(buf_id,'eden_f0',psn%eden_f0(f0_inode1:f0_inode2,1),err)
#endif

#ifdef F0_CHARGE_N0
        if(sml_mype==0) then
          if(sml_electron_on) then
            ADIOS_WRITE_LBL(buf_id,'f0_edensity_n0_add',f0_density_n0_add(:,0),err)
          endif
          ADIOS_WRITE_LBL(buf_id,'f0_idensity_n0_add',f0_density_n0_add(:,1),err)
        endif
#endif

!        ADIOS_WRITE_LBL(buf_id,'den',f0_den,err)
!        ADIOS_WRITE_LBL(buf_id,'temp_i',f0_temp(:,1),err)
!        if(sml_electron_on) then
!          ADIOS_WRITE_LBL(buf_id,'temp_e',f0_temp(:,0),err)
!        endif

        ! f data
        if(.not. use_real4) then
#ifndef F0_TOR_LINEAR
          ADIOS_WRITE_LBL(buf_id,'i_f0g',f0_f0g(:,:,:,1),err)
#else
          ADIOS_WRITE_LBL(buf_id,'i_f0g' ,f0_f0g(:,:,:,0,1),err)
          ADIOS_WRITE_LBL(buf_id,'i_f0g1',f0_f0g(:,:,:,1,1),err)
#endif
          ADIOS_WRITE_LBL(buf_id,'i_f',f0_f(:,:,:,1),err)
          ADIOS_WRITE_LBL(buf_id,'dpot',psn%dpot(f0_inode1:f0_inode2,1),err)

          if(sml_electron_on) then
#ifndef F0_TOR_LINEAR
            ADIOS_WRITE_LBL(buf_id,'e_f0g',f0_f0g(:,:,:,0),err)
#else
            ADIOS_WRITE_LBL(buf_id,'e_f0g' ,f0_f0g(:,:,:,0,0),err)
            ADIOS_WRITE_LBL(buf_id,'e_f0g1',f0_f0g(:,:,:,1,0),err)
#endif
            ADIOS_WRITE_LBL(buf_id,'e_f',f0_f(:,:,:,0),err)
          endif
          !file close
          ADIOS_CLOSE(buf_id,err)

        else
          allocate(tmp_f_real4(-f0_nvp:f0_nvp, f0_inode1:f0_inode2, f0_imu1:f0_imu2, ptl_isp:ptl_nsp))
#ifndef F0_TOR_LINEAR
          allocate(tmp_g_real4(-f0_nvp:f0_nvp, f0_inode1:f0_inode2, f0_imu1:f0_imu2, ptl_isp:ptl_nsp))
#else
          allocate(tmp_g_real4(-f0_nvp:f0_nvp, f0_inode1:f0_inode2, f0_imu1:f0_imu2, 0:1, ptl_isp:ptl_nsp))
#endif
          tmp_f_real4=f0_f
          tmp_g_real4=f0_f0g

#ifndef F0_TOR_LINEAR
          ADIOS_WRITE_LBL(buf_id,'i_f0g4',tmp_g_real4(:,:,:,1),err)
#else
          ADIOS_WRITE_LBL(buf_id,'i_f0g4' ,tmp_g_real4(:,:,:,0,1),err)
          ADIOS_WRITE_LBL(buf_id,'i_f0g41',tmp_g_real4(:,:,:,1,1),err)
#endif
          ADIOS_WRITE_LBL(buf_id,'i_f4',tmp_f_real4(:,:,:,1),err)
          ADIOS_WRITE_LBL(buf_id,'dpot',psn%dpot(f0_inode1:f0_inode2,1),err)

          if(sml_electron_on) then
#ifndef F0_TOR_LINEAR
            ADIOS_WRITE_LBL(buf_id,'e_f0g4',tmp_g_real4(:,:,:,0),err)
#else
            ADIOS_WRITE_LBL(buf_id,'e_f0g4' ,tmp_g_real4(:,:,:,0,0),err)
            ADIOS_WRITE_LBL(buf_id,'e_f0g41',tmp_g_real4(:,:,:,1,0),err)
#endif
            ADIOS_WRITE_LBL(buf_id,'e_f4',tmp_f_real4(:,:,:,0),err)
          endif
          !file close
          ADIOS_CLOSE(buf_id,err)

          deallocate(tmp_f_real4,tmp_g_real4)
        endif



     endif
  endif
#ifdef DIAG_NOISE
  call diag_f0_noise(istep,grid,psn,flag)
#endif
  
#endif
end subroutine diag_f0

#ifdef DIAG_NOISE
subroutine diag_f0_noise(istep,grid,psn,flag)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use f0_module
  implicit none
  integer, intent(in) :: istep
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer, intent(in), optional :: flag ! flag is for debug purpose
  !
  integer :: i
  integer :: nnode, inode1m1, ndata
  integer :: nphi, iphi
  integer :: nmup1, imu1m1, mudata
  integer :: vpdata
  character (len=30) :: filename
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  integer, parameter :: nspace=2
  character (len=30) :: pathname
  character (len=40) :: labelname
  real (8) :: iden_f0_approx(f0_inode1:f0_inode2), tmp(f0_inode1:f0_inode2)
  include 'mpif.h'



  if (mod(istep,diag_f0_period)==0) then
      nnode=grid%nnode
      inode1m1=f0_inode1-1
      ndata=f0_inode2-f0_inode1+1

      nphi=sml_nphi_total
      iphi=sml_intpl_mype

      nmup1= f0_nmu+1 ! include 0 index
      imu1m1=f0_imu1 -1 +1 ! additional +1 due to 0 starting index
      mudata=f0_imu2-f0_imu1+1

      vpdata=f0_nvp*2+1


     write(filename,'("xgc.f0n",".",i3.3,".",i5.5,".bp")') sml_plane_index, (sml_gstep)/(diag_f0_period)
     if(present(flag)) then
        if(flag==-1) write(filename,'("xgc.f0_debug",".",i3.3,".",i5.5,".bp")') sml_plane_index, (sml_gstep)/(diag_f0_period)
     endif
     if(sml_turb_poisson) then

        buf_size = 1000 
        buf_size = buf_size + ndata*mudata*vpdata*8*5
        
        if(sml_electron_on) then                                 
           buf_size = buf_size + ndata*mudata*vpdata*8*5 
        endif
        
        !file open
        ADIOS_OPEN(buf_id,'f0n',filename,'w',sml_plane_comm,err)
        ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
        
        ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
        ADIOS_WRITE_LBL(buf_id,'inode1m1',inode1m1, err)
        ADIOS_WRITE_LBL(buf_id,'ndata',ndata, err)
        
        ADIOS_WRITE_LBL(buf_id,'nphi',nphi,err)
        ADIOS_WRITE_LBL(buf_id,'iphi',iphi,err)
        
        
        ADIOS_WRITE_LBL(buf_id,'nmup1',nmup1, err)
        ADIOS_WRITE_LBL(buf_id,'imu1m1',imu1m1, err)
        ADIOS_WRITE_LBL(buf_id,'mudata',mudata, err)
        
        ADIOS_WRITE_LBL(buf_id,'vpdata',vpdata, err)
        

        ! f data
        ADIOS_WRITE_LBL(buf_id,'i_f0_n',  f0_diag_n(:,:,:,1) , err)
        ADIOS_WRITE_LBL(buf_id,'i_f0_w0', f0_diag_w0(:,:,:,1), err)
        ADIOS_WRITE_LBL(buf_id,'i_f0_w0s',f0_diag_w0s(:,:,:,1),err)
        ADIOS_WRITE_LBL(buf_id,'i_f0_w',  f0_diag_w(:,:,:,1),  err)
        ADIOS_WRITE_LBL(buf_id,'i_f0_ws', f0_diag_ws(:,:,:,1), err)


        if(sml_electron_on) then
           ADIOS_WRITE_LBL(buf_id,'e_f0_n',  f0_diag_n(:,:,:,0) , err)
           ADIOS_WRITE_LBL(buf_id,'e_f0_w0', f0_diag_w0(:,:,:,0), err)
           ADIOS_WRITE_LBL(buf_id,'e_f0_w0s',f0_diag_w0s(:,:,:,0),err)
           ADIOS_WRITE_LBL(buf_id,'e_f0_w',  f0_diag_w(:,:,:,0),  err)
           ADIOS_WRITE_LBL(buf_id,'e_f0_ws', f0_diag_ws(:,:,:,0), err)
        endif


        !file close
        ADIOS_CLOSE(buf_id,err)
     endif
  endif
#endif
end subroutine diag_f0_noise

#endif

! density, temperature (mean E), parallel flow from f0_F
subroutine diag_3d_f0_f(grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  use f0_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn

  integer :: isp, imu, node, ivp
  !real (8), dimension(f0_inode1:f0_inode2,diag_1d_isp:diag_1d_nsp) :: T1,T2,n,up
  real (8) , allocatable :: T1(:,:),T2(:,:), n(:,:), up(:,:)
  real (8) :: mu_vol, mu, vp, en1, en2, t_norm_ev, en_th, vth, f

  allocate(T1(f0_inode1:f0_inode2,diag_1d_isp:diag_1d_nsp))
  allocate(T2(f0_inode1:f0_inode2,diag_1d_isp:diag_1d_nsp))
  allocate(n(f0_inode1:f0_inode2,diag_1d_isp:diag_1d_nsp))
  allocate(up(f0_inode1:f0_inode2,diag_1d_isp:diag_1d_nsp))

  T1=0D0
  T2=0D0
  n=0D0
  up=0D0



  if(mod(sml_gstep,diag_3d_period)/=0) return  ! do nothing


  do isp=diag_1d_isp, diag_1d_nsp

    ! get sum of f, E*f, vp*f

    do imu=f0_imu1, f0_imu2

      if(imu==0) then
        mu_vol=0.5D0
      elseif(imu==f0_nmu) then
        mu_vol=0.5D0
      else
        mu_vol=1D0
      endif

      mu=imu*f0_dsmu
      mu=mu*mu  ! <--- This is (v_perp/v_th)^2 for v_perp grid and mu_N for sqrt(mu) grid!!!!

      do node=f0_inode1, f0_inode2
        do ivp=-f0_nvp, f0_nvp

          vp=ivp*f0_dvp
          en1= 0.5d0 * mu ! energy normalized by T
          en2= 0.5d0*vp*vp

          f=f0_f(ivp,node,imu,isp)
          T1(node,isp)=T1(node,isp) + f*en1*mu_vol
          T2(node,isp)=T2(node,isp) + f*en2*mu_vol
          n(node,isp) = n(node,isp) + f    *mu_vol
          up(node,isp)=up(node,isp) + f*vp *mu_vol
        enddo   ! ivp
      enddo     ! node
    enddo       ! imu



    ! multiply normalization
    do node=f0_inode1, f0_inode2

      !normalized v
      t_norm_ev=f0_t_ev(node,isp)
      en_th=t_norm_ev*sml_ev2j
      vth=sqrt(en_th/ptl_mass(isp))

      T1(node,isp)=T1(node,isp)/n(node,isp)*t_norm_ev
      T2(node,isp)=T2(node,isp)/n(node,isp)*t_norm_ev
      up(node,isp)=up(node,isp)/n(node,isp)*vth
      n(node,isp)=n(node,isp)*f0_grid_vol_vonly(node,isp)!/grid%node_vol_nearest(node)
    enddo
  enddo         ! isp


  ! write data as file
  call save_as_adios_format


  deallocate(T1,T2,n,up)

contains
  subroutine save_as_adios_format
    implicit none
#ifdef ADIOS
    character (len=30) :: filename
    integer*8 :: buf_id, buf_size, total_size
    integer :: err
    character (len=40) :: labelname
    integer :: nphi, iphim1, nnode, inode1m1, ndata

    nphi=sml_nphi_total
    iphim1=sml_intpl_mype

    nnode=grid%nnode
    inode1m1=f0_inode1-1
    ndata=f0_inode2-f0_inode1+1


    write(filename,'("xgc.f3d",".",i5.5,".bp")') (sml_gstep)/(diag_3d_period)
    buf_size=1000 +  4*ndata*8
    if(diag_1d_isp==0) buf_size=buf_size + 4*ndata*8

    ! Open bp file
    ADIOS_OPEN(buf_id,'f3d',filename,'w',sml_comm,err)
    ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    ADIOS_WRITE_LBL(buf_id,'nphi',nphi,err)
    ADIOS_WRITE_LBL(buf_id,'iphim1',iphim1,err)
    ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
    ADIOS_WRITE_LBL(buf_id,'inode1m1',inode1m1,err)
    ADIOS_WRITE_LBL(buf_id,'ndata',ndata,err)


    ! Write main data
    ADIOS_WRITE_LBL(buf_id,'i_T_perp',T1(:,1),err)
    ADIOS_WRITE_LBL(buf_id,'i_E_para',T2(:,1),err)
    ADIOS_WRITE_LBL(buf_id,'i_u_para',up(:,1),err)
    ADIOS_WRITE_LBL(buf_id,'i_den',n(:,1),err)

    if(diag_1d_isp==0) then
      ADIOS_WRITE_LBL(buf_id,'e_T_perp',T1(:,0),err)
      ADIOS_WRITE_LBL(buf_id,'e_E_para',T2(:,0),err)
      ADIOS_WRITE_LBL(buf_id,'e_u_para',up(:,0),err)
      ADIOS_WRITE_LBL(buf_id,'e_den',n(:,0),err)
    endif


    ! close
    ADIOS_CLOSE(buf_id,err)
#endif
  end subroutine

end subroutine


!subroutine diag_3d_additional(grid,psn,sp)
!   use sml_module
!   use grid_class
!   use psn_class
!   use ptl_module
!   use diag_module
!   use omp_module
!   implicit none
!   include 'mpif.h'
!   type(grid_type) :: grid
!   type(psn_type) :: psn
!   type(species_type) :: sp
!   !
!   real (8), allocatable :: var(:,:,:)
!   integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads), jth
!   integer :: i, j, node
!   real (8) :: b, v(diag_3d_nvar), particle_weight, wp
!   integer :: ierr
!   real (8), external :: b_interpol


!   if(.not. allocated(diag_3d_add)) then
!       allocate(diag_3d_add(diag_3d_nvar,grid%nnode,ptl_isp:ptl_nsp))
!    endif

!    allocate(var(diag_3d_nvar,grid%nnode,sml_nthreads))
!    var=0D0

!    call split_indices(sp%num, sml_nthreads, i_beg, i_end)
   
! !$OMP PARALLEL DO &
! !$OMP PRIVATE( ITH, I, B, V, PARTICLE_WEIGHT, J, NODE, &
! !$OMP          WP)
!     do ith=1, sml_nthreads
!        do i=i_beg(ith),i_end(ith)
!          if(sp%ptl(i)%gid<=0) cycle
!          if((sp%tr_save(i)) <= 0) cycle

!          b=b_interpol(sp%ptl(i)%ph(pir),sp%ptl(i)%ph(piz),sp%ptl(i)%ph(pip))

!          v(1)=1D0
!          v(2)=ptl_c_m(sp%type) * sp%ptl(i)%ph(pirho) * b  ! parallel v
!          v(3)=0.5D0*ptl_mass(sp%type)*v(1)*v(1)  + sp%ptl(i)%ct(pim)*b  ! energy

!          ! particle weight
!          if(ptl_deltaf_sp(sp%type)) then
!            particle_weight=sp%ptl(i)%ph(piw1)*sp%ptl(i)%ct(piw0)
!          else
!            particle_weight=sp%ptl(i)%ct(piw0) ! for full f simulation only
!          endif

!          do j=1, 3
!               node=grid%nd(j,sp%tr_save(i))
!               wp=sp%p_save(j,i)

!               var(:,node,ith)=var(:,node,ith) + v*particle_weight*wp
!          enddo
!        enddo
!     enddo

!     call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

! !$OMP PARALLEL DO &
! !$OMP PRIVATE( JTH, ITH, NODE)
!     do jth=1, sml_nthreads
!      if(sml_nthreads >= 2) then
!        do ith=2, sml_nthreads
!            do node=i_beg(jth), i_end(jth)
!               var(:,node,1)=var(:,node,1)+var(:,node,ith)
!            enddo
!        enddo
!      endif
!     enddo

!    ! mpi_reduce of var(:,:,1) to diag_3d_ad -- plane pe
!     call mpi_reduce(var(:,:,1),diag_3d_add(:,:,sp%type),diag_3d_nvar*grid%nnode,MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)
!     if(sml_plane_mype==0) then
! !$OMP PARALLEL DO &
! !$OMP PRIVATE( JTH, NODE)
!        do jth=1, sml_nthreads
!           do node=i_beg(jth), i_end(jth)
!               diag_3d_add(2,node,sp%type)=diag_3d_add(2,node,sp%type)/diag_3d_add(1,node,sp%type)
!               diag_3d_add(3,node,sp%type)=diag_3d_add(3,node,sp%type)/diag_3d_add(1,node,sp%type)
!           enddo
!        enddo
!     endif
!
!
! end subroutine

subroutine diag_3d_additional(grid,psn,sp)
  use sml_module
  use grid_class
  use psn_class
  use ptl_module
  use diag_module
  use omp_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  !
  real (8), allocatable :: var(:,:,:,:,:)
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads), jth
  integer :: i, j, k, node, dir, irho
  real (8) :: b, v(diag_3d_nvar), particle_weight, wp, &
               phi, wphi(0:1), wrho(0:1)
  integer :: nn, ierr
  real (8), external :: b_interpol
  ! MPI stuff
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus

  if(.not. allocated(diag_3d_add)) then
      allocate(diag_3d_add(diag_3d_nvar,grid%nnode,ptl_isp:ptl_nsp))
   endif

   allocate(var(diag_3d_nvar,grid%nnode,0:1,0:grid%nrho,max(sml_nthreads,2)))
   var=0D0

   call t_startf("DIAG_3D_PART_COLL")

   call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, DIR, I, J, K, B, V, PARTICLE_WEIGHT, NODE, &
!$OMP          WP, PHI, WPHI, WRHO, IRHO)
   do ith=1, sml_nthreads
     do i=i_beg(ith),i_end(ith)
       if(sp%ptl(i)%gid<=0) cycle
       if((sp%tr_save(i)) <= 0) cycle

       phi=modulo(sp%ptl(i)%ph(pip),grid%delta_phi)
       irho=floor(sp%rhoi(i)/grid%drho)
       irho=min(irho,grid%nrho-1)
       wp=min(sp%rhoi(i)/grid%drho-irho,1d0)
       wrho= (/ 1d0-wp, wp /)
       b=b_interpol(sp%ptl(i)%ph(pir),sp%ptl(i)%ph(piz),sp%ptl(i)%ph(pip))
       wp=phi/grid%delta_phi
       wphi= (/ 1d0-wp, wp /)
  ! density
       v(1)=1D0
  ! energy density
       v(2)=ptl_c2_2m(sp%type)*(b*sp%ptl(i)%ph(pirho))**2+sp%ptl(i)%ct(pim)*b

       ! particle weight
       if(ptl_deltaf_sp(sp%type)) then
         particle_weight=sp%ptl(i)%ph(piw1)*sp%ptl(i)%ct(piw0)
       else
         particle_weight=sp%ptl(i)%ct(piw0) ! for full f simulation only
       endif

       do dir=0,1
         do j=1, 3
           node=grid%nd(j,sp%tr_save(i))
           v(3)=psn%pot_rho_ff(dir,0,node)
           do k=0,1
             wp=sp%p_save(j,i)*wphi(dir)*wrho(k)/grid%node_vol_ff(node,dir)
             var(:,node,dir,irho+k,ith)=var(:,node,dir,irho+k,ith)+v*particle_weight*wp
           enddo
         enddo
       enddo
     enddo
   enddo

   call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( JTH, ITH, NODE)
    do jth=1, sml_nthreads
       do ith=2, sml_nthreads
           do node=i_beg(jth), i_end(jth)
              var(:,node,:,:,1)=var(:,node,:,:,1)+var(:,node,:,:,ith)
           enddo
       enddo
    enddo

   call t_stopf("DIAG_3D_PART_COLL")

   call t_startf("DIAG_3D_COMM")

   nn=grid%nnode*diag_3d_nvar
   do dir=0,1
     do irho=0,grid%nrho
       ith=max((irho-1)*2+dir,0)
       call MPI_reduce(var(:,:,dir,irho,1),var(:,:,dir,irho,2),nn, &
            MPI_REAL8,MPI_SUM,ith,SML_PLANE_COMM,ierr)
     enddo
   enddo

   call t_stopf("DIAG_3D_COMM")

   call t_startf("DIAG_3D_GYROFF")

    do j=1,diag_3d_nvar
      do k=0,grid%nrho
        do dir=0,1
          if(sml_plane_mype==max(0,(k-1)*2)+dir) then
            if(k>0) then
              call mat_mult(psn%gyro_avg_mat,var(j,:,dir,k,2),var(j,:,dir,k,1))
            else
              ! here sml_plane_mype==1 doesn't do anything
              var(j,:,dir,k,1)=var(j,:,dir,k,2)
            endif
          endif
        enddo
      enddo
    enddo

   call t_stopf("DIAG_3D_GYROFF")

   call t_startf("DIAG_3D_COMM")

    if(sml_plane_mype/=0) then
       if(sml_plane_mype<2*grid%nrho) then
         isendtag=sml_plane_mype
         irho=sml_plane_mype/2+1
	 dir=mod(sml_plane_mype,2)
         idest=0
         call MPI_send(var(:,:,dir,irho,1),nn,MPI_REAL8,idest,isendtag,SML_PLANE_COMM,ierr)
       endif
    else ! sml_plane_mype==0
       do isource=1,2*grid%nrho-1
         irecvtag=isource
         irho=isource/2+1
	 dir=mod(isource,2)
         call MPI_recv(var(:,:,dir,irho,1),nn,MPI_REAL8,isource,irecvtag,SML_PLANE_COMM,istatus,ierr)
       enddo
    endif

!    I had an OMP loop here, but why bother?

    if(sml_plane_mype==0) then
      var(:,:,:,0,1:2)=sum(var(:,:,:,:,1:2),dim=4)
      do j=1,diag_3d_nvar
        call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,var(j,:,:,0,1),var(j,:,:,1,1))
        call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,var(j,:,:,0,2),var(j,:,:,1,2))
      enddo
    endif
    
    ! mpi_reduce of var(:,:,1) to diag_3d_var
    if(sml_plane_mype==0) then
     nn=grid%nnode*diag_3d_nvar
     idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
     isendtag=sml_intpl_mype

     isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
     irecvtag=isource
     ! this is for gyroaveraged
     call mpi_sendrecv(var(:,:,0,1,1),nn,MPI_REAL8,  idest,isendtag, &
          var(:,:,1,0,1),nn,MPI_REAL8,isource,irecvtag, &
          sml_intpl_comm,istatus,ierr)
     ! this is for gyrocentered
     call mpi_sendrecv(var(:,:,0,1,2),nn,MPI_REAL8,  idest,isendtag, &
          var(:,:,1,0,2),nn,MPI_REAL8,isource,irecvtag, &
          sml_intpl_comm,istatus,ierr)


     diag_3d_add(1:2,:,sp%type)=0.5d0*(var(1:2,:,1,0,1)+var(1:2,:,1,1,1))
     diag_3d_add(3,:,sp%type)=0.5d0*(var(1,:,1,0,2)+var(1,:,1,1,2))

!$OMP PARALLEL DO &
!$OMP PRIVATE( JTH, NODE)
       do jth=1, sml_nthreads
          do node=i_beg(jth), i_end(jth)
              diag_3d_add(2,node,sp%type)=diag_3d_add(2,node,sp%type)/diag_3d_add(1,node,sp%type)
               ! it's density now
!              diag_3d_add(3,node,sp%type)=diag_3d_add(3,node,sp%type)/diag_3d_add(1,node,sp%type)
          enddo
       enddo
    endif

   call t_stopf("DIAG_3D_COMM")

   deallocate(var)
   !
end subroutine

!*********************************************************
! single particle tracer
!*********************************************************
subroutine tracer(n,sp_type,grid,psn,spall)
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  type(fld_type) :: fld
  integer, intent(in) :: n,sp_type
  type(species_type) , pointer :: sp
  integer :: i, m
  real (kind=8) :: r,z,phi,rho,mu,b,psi,psi_c,b_interpol,ev,pitch,br,bz,bphi
  real (kind=8) :: exb1, exb2
  logical :: rz_outside
  real (8), external :: psi_interpol
  integer (kind=8) :: gid






#ifdef ADIOS
  real (kind=8) :: tracerdata(18)
  integer :: ncol_tracer, err
  integer*8:: buf_id, buf_size, total_size
#endif


  !find particle location (in array) to write out
  m=0
  do i=1, spall(sp_type)%num
    gid=spall(sp_type)%ptl(i)%gid
    if(gid==n) then
      m=i
      exit
    endif
  enddo

  if(m==0) return    ! no particle in this mpi processor


  !gid=spall(sp_type)%ptl(m)%gid
  if (gid>0) then
     fld%r=spall(sp_type)%ptl(m)%ph(1)
     fld%z=spall(sp_type)%ptl(m)%ph(2)
     fld%phi=spall(sp_type)%ptl(m)%ph(3)
     fld%psi=psi_interpol(fld%r,fld%z,0,0)

     call field(fld,sml_time,rz_outside)
     call efield(grid,psn,spall(1),m,fld,sml_time)
!     if (sml_deltaf) call f0_info(fld)
     r=fld%r
     z=fld%z
     phi=fld%phi
     rho=spall(sp_type)%ptl(m)%ph(4)
     mu=spall(sp_type)%ptl(m)%ct(pim)
     psi=fld%psi

     b=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     bphi=fld%bphi
     br=fld%br
     bz=fld%bz

     call rho_mu_to_ev_pitch2(rho,mu,b,ev,pitch,sp_type)
     if(sml_bp_sign > 0) then
        psi_c=psi + rho*bphi*r
     else
        psi_c=psi - rho*bphi*r
     endif

     exb1= (fld%bz*fld%Ephi - fld%Bphi * fld%Ez)/B/B
     exb2= (fld%bphi * fld%Er - fld%br * fld%Ephi )/B/B
#ifdef ADIOS
     tracerdata(1)=sml_time/sml_tran
     tracerdata(2)=r
     tracerdata(3)=z
     tracerdata(4)=phi
     tracerdata(5)=rho
     tracerdata(6)=ev
     tracerdata(7)=ev+fld%Epot
     tracerdata(8)=pitch
     tracerdata(9)=psi/eq_x_psi
     tracerdata(10)=psi_c
     tracerdata(11)=mu
     tracerdata(12)=spall(sp_type)%ptl(n)%ph(piw1)
!#     tracerdata(13)=spall(sp_type)%derivs(1,n)
!#     tracerdata(14)=spall(sp_type)%derivs(2,n)
!#     tracerdata(15)=spall(sp_type)%derivs(4,n)
     tracerdata(16)=real(gid)
     tracerdata(17)=(exb1*fld%dpsidr + exb2*fld%dpsidz)/sqrt(fld%dpsidr**2+fld%dpsidz**2)
!#     tracerdata(18)=(spall(sp_type)%derivs(1,n)*fld%Er+spall(sp_type)%derivs(2,n)*fld%Ez+spall(sp_type)%derivs(3,n)*fld%r*fld%Ephi)/fld%f0_temp
     ncol_tracer=18
     buf_size = 4 + ncol_tracer * 8
!    print *, 'adios writing fort.tracer group, #bytes = ', buf_size
     ADIOS_OPEN(buf_id,'fort.tracer','xgc.tracer.bp','a',sml_comm_null,err)
     ADIOS_SET_PATH(buf_id,'/tracer',err)
     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
     ADIOS_WRITE_LBL(buf_id,'ncol_tracer',ncol_tracer,err)
     ADIOS_WRITE_LBL(buf_id,'tracerdata',tracerdata,err)
     ADIOS_CLOSE(buf_id,err)

!#else
!     open(unit=22, file='fort.tracer', position='append')
     open(unit=22, file='tracer.out', position='append')
     write(22,1000) sml_time/sml_tran,r, z, phi, rho,&
          ev,ev+fld%Epot,pitch,psi/eq_x_psi,psi_c,mu,spall(sp_type)%ptl(n)%ph(piw1), &
!#          spall(sp_type)%derivs(1,n),spall(sp_type)%derivs(2,n),spall(sp_type)%derivs(4,n),&
!#!#        real(gid), (exb1*fld%dpsidr + exb2*fld%dpsidz)/sqrt(fld%dpsidr**2+fld%dpsidz**2), &
          real(gid), (exb1*fld%dpsidr + exb2*fld%dpsidz)/sqrt(fld%dpsidr**2+fld%dpsidz**2)
!#          (spall(sp_type)%derivs(1,n)*fld%Er+&
!#          spall(sp_type)%derivs(2,n)*fld%Ez+&
!#          spall(sp_type)%derivs(3,n)*fld%r*fld%Ephi)/fld%f0_temp
     close(22)
#endif
  endif
1000 format(18(e19.13,1x))
end subroutine tracer



#if defined(ADIOS)
subroutine dump_grid(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: m, n, dir, itr, n_n, n_t,  dims(1), n_geo
  integer, parameter :: ndim=1
  real (kind=8), pointer :: coord(:,:)
  integer, pointer :: nodeid(:,:), nextn(:)
  character (len=30) :: meshbp_file, varname
  integer (kind=8) :: meshbp_fileptr
  integer (kind=8) :: buf_id, buf_size, total_size
  integer :: err

  n_n = grid%nnode
  n_t = grid%ntriangle
  allocate(coord(2,n_n),nodeid(3,n_t))
! node coords scaled to physical length unit of meters
  coord(:,:) = grid%x(:,:)
! node id's with zero-base indexing for C code
  nodeid(:,:) = grid%nd(:,:) - 1

  allocate(nextn(n_n))
  if(sml_turb_poisson) then
! use bfollow data arrays to determine node connectivity from one plane to next
    dir = 1
#ifdef OLD_INIT_FF
    if (sml_bt_sign < 0) dir=2
    do n=1,n_n
      itr = psn%bfollow_tr(dir,n)
      m = 1
      if (psn%bfollow_p(2,dir,n) > psn%bfollow_p(m,dir,n)) m=2
      if (psn%bfollow_p(3,dir,n) > psn%bfollow_p(m,dir,n)) m=3
! node id's with zero-base indexing for C code
      nextn(n) = grid%nd(m,itr) - 1
#else
! ff_1dp_p is always phi+ oriented
    do n=1,n_n
      itr = psn%ff_1dp_tr(n,dir)
      m = 1
      if (psn%ff_1dp_p(2,n,dir) > psn%ff_1dp_p(m,n,dir)) m=2
      if (psn%ff_1dp_p(3,n,dir) > psn%ff_1dp_p(m,n,dir)) m=3
! node id's with zero-base indexing for C code
      nextn(n) = grid%nd(m,itr) - 1
#endif
    enddo
  else
! just assume toroidal symmetry in mesh
    do n=1,n_n
! node id's with zero-base indexing for C code
      nextn(n) = n - 1
    enddo
  endif

  buf_size = 4 + 4 + 1000 ! 1KB buffer
  buf_size = buf_size + 2 * n_n * 8 * 2 ! /coordinates/values and rz
  buf_size = buf_size + 3 * n_t * 4 * 2  ! connectivity - nodeid -- zero based.
  buf_size = buf_size + n_n * 4 + 5 * n_n * 8  ! next node, psi, node_vol, ff0, ff1, nearest
  buf_size = buf_size + n_t*8  + n_t*2*3*8 + n_t*4*3*8  ! tr_area , mapping, ff_{h,1}dp_p
  buf_size = buf_size + 4 + psn%nwall*4

  buf_size = buf_size + grid%npsi00 * 5 * 8  !! qsafety, trapped, inv. aspect ratio
  n_geo=grid%npsi00

  buf_size = buf_size + n_n*8 + n_t*4  ! additonal buffer for wrong buf_size calculation -- I could not find it.- sk.

#ifndef OLD_INIT_FF
  buf_size = buf_size + 2*8*n_n
#endif
! if (sml_mype==0) print *, 'adios writing diagnosis.mesh group, #bytes = ', buf_size
  ADIOS_OPEN(buf_id,'diagnosis.mesh','xgc.mesh.bp','w',sml_comm_null,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'n_n',n_n,err)
  ADIOS_WRITE_LBL(buf_id,'n_t',n_t,err)
  ADIOS_WRITE_LBL(buf_id,'n_geo',n_geo,err)
  ADIOS_WRITE_LBL(buf_id,'/coordinates/values',coord,err)  !redundant for compatibility
  ADIOS_WRITE_LBL(buf_id,'rz',coord,err)  
  ADIOS_WRITE_LBL(buf_id,'/cell_set[0]/node_connect_list',nodeid,err) !redunant for compatibility
  ADIOS_WRITE_LBL(buf_id,'nd_connect_list',nodeid,err) 
  ADIOS_WRITE_LBL(buf_id,'nextnode',nextn,err)
#ifdef OLD_INIT_FF
  ADIOS_WRITE_LBL(buf_id,'one_per_dx',psn%bfollow_1_dx,err)
#else
  ADIOS_WRITE_LBL(buf_id,'one_per_dx',psn%ff_1dp_dx,err)
#endif
  ADIOS_WRITE_LBL(buf_id,'psi',grid%psi,err)
  ADIOS_WRITE_LBL(buf_id,'node_vol', grid%node_vol,err)
  ADIOS_WRITE_LBL(buf_id,'node_vol_ff0', grid%node_vol_ff(:,0),err)
  ADIOS_WRITE_LBL(buf_id,'node_vol_ff1', grid%node_vol_ff(:,1),err)
  ADIOS_WRITE_LBL(buf_id,'tr_area',grid%tr_area,err)
  ADIOS_WRITE_LBL(buf_id,'mapping',grid%mapping,err)
  ADIOS_WRITE_LBL(buf_id,'ff_hdp_p',psn%ff_hdp_p,err)
  ADIOS_WRITE_LBL(buf_id,'ff_1dp_p',psn%ff_1dp_p,err)
!  ADIOS_WRITE_LBL(buf_id,'ff_hdp_p1',psn%ff_hdp_p(:,:,1),err)


  if(sml_f0_grid) ADIOS_WRITE_LBL(buf_id,'node_vol_nearest',grid%node_vol_nearest,err)

  !! Minor/major radius, inv. aspect ratio, qsafety, trapped particle fraction
  ADIOS_WRITE_LBL(buf_id,'rmin', grid%epspar(:,1),err)
  ADIOS_WRITE_LBL(buf_id,'rmaj', grid%epspar(:,2),err)
  ADIOS_WRITE_LBL(buf_id,'epsilon', grid%epspar(:,3),err)
  ADIOS_WRITE_LBL(buf_id,'qsafety', grid%qsafety,err)
  ADIOS_WRITE_LBL(buf_id,'trapped', grid%trapped,err)


  if(sml_sheath_mode/=0) then
     ADIOS_WRITE_LBL(buf_id,'psn_nwall', psn%nwall,err)
     ADIOS_WRITE_LBL(buf_id,'psn_wall_nodes', psn%wall_nodes,err)
  endif

  ADIOS_CLOSE(buf_id,err)
  deallocate(coord,nodeid,nextn)
end subroutine dump_grid

subroutine dump_f0_grid(grid)
  use grid_class
  use sml_module
  use f0_module
  use ptl_module, only: ptl_nsp, ptl_isp
  implicit none
  type(grid_type) :: grid
  integer :: n_n, ndata, inode1m1, nsp, i
  integer (kind=8) :: buf_id, buf_size, total_size
  integer :: err
  real (kind=8), allocatable :: gradpsi(:,:)
  real (kind=8), external :: psi_interpol

  n_n=grid%nnode
  inode1m1=f0_inode1-1
  ndata=f0_inode2-f0_inode1+1
  nsp=ptl_nsp-ptl_isp+1

  ! Evaluate grad(Psi)
  allocate(gradpsi(2,f0_inode1:f0_inode2))
  gradpsi=0D0
  do i=f0_inode1,f0_inode2
    gradpsi(1,i)=psi_interpol(grid%x(1,i),grid%x(2,i),1,0)
    gradpsi(2,i)=psi_interpol(grid%x(1,i),grid%x(2,i),0,1)
  enddo

  ! dump f0-grid quantities
  buf_size = 1000 + 6 * 4 + 4 * 8
  buf_size = buf_size + 2 * ndata * nsp * 8 + ndata * 8  !! f0_den, f0_T_ev, f0_grid_vol
  buf_size = buf_size + 2 * ndata * 8 + 2 * 3 * ndata * 8 + ndata * 8 !! v_grad, v_curv, ...
  total_size=0

  ADIOS_OPEN(buf_id,'diagnosis.f0.mesh','xgc.f0.mesh.bp','w',sml_plane_comm,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)

  ADIOS_WRITE_LBL(buf_id,'n_n',n_n,err)
  ADIOS_WRITE_LBL(buf_id,'nsp',nsp,err)
  ADIOS_WRITE_LBL(buf_id,'ndata',ndata,err)
  ADIOS_WRITE_LBL(buf_id,'inode1m1',inode1m1,err)
  ADIOS_WRITE_LBL(buf_id,'f0_nmu',f0_nmu,err)
  ADIOS_WRITE_LBL(buf_id,'f0_nvp',f0_nvp,err)
  ADIOS_WRITE_LBL(buf_id,'f0_smu_max',f0_smu_max,err)
  ADIOS_WRITE_LBL(buf_id,'f0_vp_max',f0_vp_max,err)
  ADIOS_WRITE_LBL(buf_id,'f0_dsmu',f0_dsmu,err)
  ADIOS_WRITE_LBL(buf_id,'f0_dvp',f0_dvp,err)
  ADIOS_WRITE_LBL(buf_id,'f0_den',f0_den(f0_inode1:f0_inode2),err)
  ADIOS_WRITE_LBL(buf_id,'f0_T_ev',f0_T_ev(f0_inode1:f0_inode2,:),err)
  ADIOS_WRITE_LBL(buf_id,'f0_grid_vol_vonly',f0_grid_vol_vonly(f0_inode1:f0_inode2,:),err)

  !! This is needed for post-process evaluation of moments of f0
  ADIOS_WRITE_LBL(buf_id,'v_gradb',grid%v_gradb(:,f0_inode1:f0_inode2),err)
  ADIOS_WRITE_LBL(buf_id,'v_curv',grid%v_curv(:,f0_inode1:f0_inode2),err)
  ADIOS_WRITE_LBL(buf_id,'nb_curl_nb',grid%nb_curl_nb(f0_inode1:f0_inode2),err)
  ADIOS_WRITE_LBL(buf_id,'gradpsi',gradpsi(:,f0_inode1:f0_inode2),err)

  ADIOS_CLOSE(buf_id,err)

  deallocate(gradpsi)

end subroutine dump_f0_grid


subroutine dump_bfield(grid)
  use grid_class
  use eq_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  integer :: i,n_n,dims(2)
  real (kind=8) :: r,z,phi,br,bz,bphi
  real (kind=8), allocatable :: bfld(:,:)
  integer, parameter :: nvar=2, veclen=3, flddim=2, nspace=2, veclen2=1,flddim2=1
  character (len=30) :: filename, fldname, units, varname
  integer (kind=8) :: fileptr

  integer(kind=8) :: buf_id, buf_size, total_size
  integer :: err


  n_n = grid%nnode
  allocate(bfld(3,n_n))
  phi = 3.141592/2.
  do i=1,n_n
     r = grid%x(1,i)
     z = grid%x(2,i)
     call bvec_interpol(r,z,phi,br,bz,bphi)
     bfld(1,i) = br
     bfld(2,i) = bz
     bfld(3,i) = bphi
  enddo


  buf_size = 4 + 3 * n_n * 8 + n_n * 8
! if (sml_mype==0) print *, 'adios writing diagnosis.bfield group, #bytes = ', buf_size
  ADIOS_OPEN(buf_id,'diagnosis.bfield','xgc.bfield.bp','w',sml_comm_null,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'n_n',n_n,err)
  ADIOS_WRITE_LBL(buf_id,'/node_data[0]/values',bfld,err)
  ADIOS_WRITE_LBL(buf_id,'/node_data[1]/values',grid%psi,err)
  ADIOS_CLOSE(buf_id,err)

  deallocate(bfld)
end subroutine dump_bfield
#endif

#if defined(ADIOS)
subroutine restart_write(spall,grid, psn)
  use f0_module
  use sml_module
  use ptl_module
  use psn_class
  use grid_class
  implicit none
  include 'mpif.h'
  type(species_type):: spall(0:ptl_nsp_max)
  type(grid_type) :: grid
  type(psn_type):: psn
  integer :: i,j
  character (len=50) :: filename, dirname
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  real*8 start_time, end_time !SAK added this for timers
  integer :: np
  real (8), allocatable :: phase(:,:), ephase(:,:)
  integer (8), allocatable :: gid(:), egid(:)
  integer, parameter :: ict1=ptl_nphase+1
  integer, parameter :: ict2=ptl_nphase+ptl_nconst
  !
  ! f0 grid
  integer :: nnode, inode1m1
  integer :: nmup1, imu1m1
  integer :: ndata, mudata, vpdata
#ifdef ADIOS_TIME
  real(kind=8) :: t_start, t_end, t_elap, t_sum, buf_sum
  real(kind=8) :: t_elapall(sml_totalpe)
  integer*8 :: buf_sizeall(sml_totalpe)
  logical, save :: isfirst = .true., isfirst_f0 = .true.
#endif

  integer*8 :: inum_total, ioff, enum_total, eoff, inum, enum
  integer*8 :: inum_all(sml_totalpe), enum_all(sml_totalpe)

  integer :: adios_comm, adios_comm_totalpe, adios_comm_mype

  !for EM hybrid
  nnode=grid%nnode

  !SAK DEBUGGING LINE   write(6,'("write",a45,a15,i5)') filename,dirname,sml_mype
  buf_size = 4 + 8 + 4 + 4 + 4 + 8 + 8*spall(1)%num + 8*(ptl_nphase+ptl_nconst)*spall(1)%num + 2*8 + 1000
  if (sml_sheath_mode==1) then
     buf_size = buf_size + psn%nwall * 8 + 4
  endif

  if(sml_electron_on)then
    buf_size = buf_size + 4 + 4 + 8 + 8*spall(0)%num + 8*(ptl_nphase+ptl_nconst)*spall(0)%num + 2*8
  endif

  if (sml_mype==0) print *, 'adios writing restart group, #bytes = ', buf_size
  np=spall(1)%num
  allocate(phase(ict2,np),gid(np))
  ! copy to temp memory

  do i=1, np
     phase(1:ptl_nphase,i)=spall(1)%ptl(i)%ph
     phase(ict1:ict2,i)=spall(1)%ptl(i)%ct
     gid(i)=spall(1)%ptl(i)%gid
  enddo

  if(sml_electron_on)then
     np=spall(0)%num
     allocate(ephase(ict2,np),egid(np))
     ! copy to temp. memory
     do i=1, np
        ephase(1:ptl_nphase,i)=spall(0)%ptl(i)%ph
        ephase(ict1:ict2   ,i)=spall(0)%ptl(i)%ct
        egid(i)=spall(0)%ptl(i)%gid
     enddo
  endif

#ifdef ADIOS_TIME
  call mpi_barrier(sml_comm, err);
  t_start = MPI_Wtime();
#endif
  inum = spall(1)%num
  call mpi_allgather(inum,1,MPI_INTEGER8,inum_all,1,MPI_INTEGER8,sml_comm,err)
  inum_total = sum(inum_all)
  ioff = sum(inum_all(1:(sml_mype+1)))-spall(1)%num

  if (adios_mfilesys) then
     write(filename,'("restart_dir",i1.1,"/xgc.restart.",i5.5,".bp")') sml_mfilesys_index+1, sml_gstep
     ADIOS_OPEN(buf_id,'restart',filename,'w',sml_mfilesys_comm,err)
     adios_comm = sml_mfilesys_comm
     adios_comm_totalpe = sml_pe_per_filesys
     adios_comm_mype = sml_mfilesys_mype
  else
     if (adios_stage_restart) then
        write(filename,'("restart_dir/xgc.restart.bp")')
     else
        write(filename,'("restart_dir/xgc.restart.",i5.5,".bp")') sml_gstep
     endif
     ADIOS_OPEN(buf_id,'restart',filename,'w',sml_comm,err)
     adios_comm = sml_comm
     adios_comm_totalpe = sml_totalpe
     adios_comm_mype = sml_mype
  endif

  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'timestep',sml_gstep,err)
  ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
  if (sml_sheath_mode==1) then
    ADIOS_WRITE_LBL(buf_id,'sheath_nwall',psn%nwall,err)
    ADIOS_WRITE_LBL(buf_id,'sheath_pot',psn%sheath_pot,err)
  endif
  ! ion
  ADIOS_WRITE_LBL(buf_id,'maxnum',spall(1)%maxnum,err)
  ADIOS_WRITE_LBL(buf_id,'inum',spall(1)%num,err)
  ADIOS_WRITE_LBL(buf_id,'inphase',ict2,err)
  ADIOS_WRITE_LBL(buf_id,'imaxgid',spall(1)%maxgid,err)
  
  ADIOS_WRITE_LBL(buf_id,'inum_total',inum_total,err)
  ADIOS_WRITE_LBL(buf_id,'ioff',ioff,err)
  ADIOS_WRITE_LBL(buf_id,'iphase',phase,err)
  ADIOS_WRITE_LBL(buf_id,'igid',gid,err)

  ADIOS_WRITE_LBL(buf_id,'totalpe',adios_comm_totalpe,err)
  ADIOS_WRITE_LBL(buf_id,'mype',adios_comm_mype,err)
  ADIOS_WRITE_LBL(buf_id,'inum_arr',spall(1)%num,err)
  ADIOS_WRITE_LBL(buf_id,'ioff_arr',ioff,err)
  
  if(sml_electron_on)then
     enum = spall(0)%num
     call mpi_allgather(enum,1,MPI_INTEGER8,enum_all,1,MPI_INTEGER8,sml_comm,err)
     enum_total = sum(enum_all)
     eoff = sum(enum_all(1:(adios_comm_mype+1)))-spall(0)%num
     
     ADIOS_WRITE_LBL(buf_id,'enum',spall(0)%num,err)
     ADIOS_WRITE_LBL(buf_id,'enphase',ict2,err)
     ADIOS_WRITE_LBL(buf_id,'emaxgid',spall(0)%maxgid,err)

     ADIOS_WRITE_LBL(buf_id,'enum_total',enum_total,err)
     ADIOS_WRITE_LBL(buf_id,'eoff',eoff,err)
     ADIOS_WRITE_LBL(buf_id,'ephase',ephase,err)
     ADIOS_WRITE_LBL(buf_id,'egid',egid,err)

     ADIOS_WRITE_LBL(buf_id,'enum_arr',spall(0)%num,err)
     ADIOS_WRITE_LBL(buf_id,'eoff_arr',eoff,err)
  endif
  ADIOS_CLOSE(buf_id,err)

  if(sml_electron_hyb) then
    if(sml_plane_mype==0) then
      buf_size=1000 + nnode*8*3
      write(filename,'("restart_dir/xgc.restart_hyb.",i5.5,".bp")') sml_gstep
      write(dirname,'("/node_",i5.5)') sml_mype
      ADIOS_OPEN(buf_id,'restart_hyb',filename,'w',sml_intpl_comm,err)
      ADIOS_SET_PATH(buf_id,dirname,err)
      ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
      ADIOS_WRITE_LBL(buf_id,'timestep',sml_gstep,err)
      ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)

      ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
      ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
      ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
#ifdef XGC1_EM
      ADIOS_WRITE_LBL(buf_id,'eden_hyb',psn%eden_hyb,err)
      ADIOS_WRITE_LBL(buf_id,'a_par',psn%A_par,err)
      ADIOS_WRITE_LBL(buf_id,'u_e',psn%u_e,err)
      ADIOS_CLOSE(buf_id,err)
#endif
    endif
  endif

#ifdef ADIOS_TIME
  call mpi_barrier(sml_comm, err);
  t_end = MPI_Wtime();
  t_elap = t_end - t_start
  
  call mpi_gather(t_elap, 1, MPI_REAL8, t_elapall, 1, MPI_REAL8, 0, sml_comm, err)
  call mpi_gather(buf_size, 1, MPI_INTEGER8, buf_sizeall, 1, MPI_INTEGER8, 0, sml_comm, err)
  
  if (sml_mype==0) then
     if (isfirst) then
        open(unit=99, file='timing_restart.txt',status='replace')
        isfirst=.false.
     else
        open(unit=99, file='timing_restart.txt',position='append')
     endif
     
     do i=1,sml_totalpe
        write (99,'(I8,1X,I8,1X,I15,1X,F15.3)') sml_gstep, i, buf_sizeall(i), t_elapall(i)
     enddo
     
     close(99)
     
     print '(1X, A, 3ES15.3)', "restart: bytes,sec,throughput  ", &
          real(sum(buf_sizeall)), sum(t_elapall)/size(t_elapall), &
          sum(buf_sizeall)/sum(t_elapall)*size(t_elapall)
  endif
#endif

  deallocate(phase,gid)
  if(sml_electron_on)   deallocate(ephase,egid)

  call MPI_BARRIER(sml_comm,err) !SAK


  ! f0_g restart
  if(sml_f0_grid) then
     nnode=grid%nnode
     inode1m1=f0_inode1-1
     ndata=f0_inode2-f0_inode1+1

     nmup1= f0_nmu+1 ! include 0 index
     imu1m1=f0_imu1 -1 +1 ! additional +1 due to 0 starting index
     mudata=f0_imu2-f0_imu1+1

     vpdata=f0_nvp*2+1

#ifndef F0_TOR_LINEAR
     buf_size = 4 + 8 + 4*7 + 8*vpdata*ndata*mudata
     if(sml_electron_on) buf_size = buf_size + 8*vpdata*ndata*mudata
#else
     buf_size = 4 + 8 + 4*7 + 8*vpdata*ndata*mudata*2
     if(sml_electron_on) buf_size = buf_size + 8*vpdata*ndata*mudata*2
#endif
     buf_size = buf_size + 1000
#ifdef F0_CHARGE_N0
     if(sml_mype==0) then
        buf_size = buf_size + 8*nnode
        if(sml_electron_on) buf_size = buf_size + 8*nnode
      endif
#endif
#ifdef ADIOS_TIME
     call mpi_barrier(sml_comm, err);
     t_start = MPI_Wtime();
#endif

     if (adios_mfilesys) then
        write(filename,'("restart_dir",i1.1,"/xgc.restart.f0.",i5.5,".bp")') sml_mfilesys_index+1, sml_gstep

        ADIOS_OPEN(buf_id,'restartf0',filename,'w',sml_mfilesys_comm,err)
     else
        if (adios_stage_restartf0) then
           write(filename,'("restart_dir/xgc.restart.f0.bp")')
        else
           write(filename,'("restart_dir/xgc.restart.f0.",i5.5,".bp")') sml_gstep
        endif

        ADIOS_OPEN(buf_id,'restartf0',filename,'w',sml_comm,err)
     endif

     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)

     ! Common parameters
     ADIOS_WRITE_LBL(buf_id,'timestep',sml_gstep,err)
     ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)

     ADIOS_WRITE_LBL(buf_id,'vpdata',vpdata,err)
     ADIOS_WRITE_LBL(buf_id,'nmup1',nmup1, err)
     ADIOS_WRITE_LBL(buf_id,'imu1m1',imu1m1, err)
     ADIOS_WRITE_LBL(buf_id,'mudata',mudata,err)
     ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
     ADIOS_WRITE_LBL(buf_id,'inode1m1',inode1m1, err)
     ADIOS_WRITE_LBL(buf_id,'ndata',ndata,err)

     ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
     ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)

     ! Ion f0g
#ifndef F0_TOR_LINEAR
     ADIOS_WRITE_LBL(buf_id,'i_f0g',f0_f0g(:,:,:,1),err)
#else
     ADIOS_WRITE_LBL(buf_id,'i_f0g',f0_f0g(:,:,:,0,1),err)
     ADIOS_WRITE_LBL(buf_id,'i_f0g1',f0_f0g(:,:,:,1,1),err)
#endif
     ! Electron f0g
     if(sml_electron_on) then
#ifndef F0_TOR_LINEAR
        ADIOS_WRITE_LBL(buf_id,'e_f0g',f0_f0g(:,:,:,0),err)
#else
        ADIOS_WRITE_LBL(buf_id,'e_f0g',f0_f0g(:,:,:,0,0),err)
        ADIOS_WRITE_LBL(buf_id,'e_f0g1',f0_f0g(:,:,:,1,0),err)
#endif
     endif

     ADIOS_WRITE_LBL(buf_id,'totalpe',adios_comm_totalpe,err)
     ADIOS_WRITE_LBL(buf_id,'mype',adios_comm_mype,err)
     ADIOS_WRITE_LBL(buf_id,'inode1',f0_inode1,err)
     ADIOS_WRITE_LBL(buf_id,'inode2',f0_inode2,err)
     ADIOS_WRITE_LBL(buf_id,'imu1',f0_imu1,err)
     ADIOS_WRITE_LBL(buf_id,'imu2',f0_imu2,err)
     ADIOS_WRITE_LBL(buf_id,'ndata_arr',ndata,err)

#ifdef F0_CHARGE_N0
     if(sml_mype==0) then
     if(sml_electron_on) then
        ADIOS_WRITE_LBL(buf_id,'f0_edensity_n0_add',f0_density_n0_add(:,0),err)
     endif
     ADIOS_WRITE_LBL(buf_id,'f0_idensity_n0_add',f0_density_n0_add(:,1),err)
     endif
#endif

     ADIOS_CLOSE(buf_id,err)

#ifdef ADIOS_TIME
     call mpi_barrier(sml_comm, err);
     t_end = MPI_Wtime();
     t_elap = t_end - t_start

     call mpi_gather(t_elap, 1, MPI_REAL8, t_elapall, 1, MPI_REAL8, 0, sml_comm, err)
     call mpi_gather(buf_size, 1, MPI_INTEGER8, buf_sizeall, 1, MPI_INTEGER8, 0, sml_comm, err)

     if (sml_mype==0) then
        if (isfirst_f0) then
           open(unit=99, file='timing_restartf0.txt',status='replace')
           isfirst_f0=.false.
        else
           open(unit=99, file='timing_restartf0.txt',position='append')
        endif
        
        do i=1,sml_totalpe
           write (99,'(I8,1X,I8,1X,I15,1X,F15.3)') sml_gstep, i, buf_sizeall(i), t_elapall(i)
        enddo
        
        close(99)
        
        print '(1X, A, 3ES15.3)', "restartf0: bytes,sec,throughput", &
             real(sum(buf_sizeall)), sum(t_elapall)/size(t_elapall), &
             sum(buf_sizeall)/sum(t_elapall)*size(t_elapall)
     endif
#endif

     call MPI_BARRIER(sml_comm,err) !SAK
  endif


  if(sml_mype==0) then
     open(unit=1,file='timestep.dat',status='replace')
     write(1,1000) sml_gstep
     write(1,1000) sml_run_count
     close(1)
  endif
1000 format(I8)

  return
end subroutine restart_write

subroutine restart_read(grid,psn,spall)
  use sml_module
  use ptl_module
  use f0_module
  use grid_class
  use psn_class
  use adios_read_mod
  implicit none
  type(grid_type):: grid
  type(psn_type):: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i,ierr,err
  character (len=50) :: filename,dirname
  integer*8 :: buf_id, buf_size, total_size
  integer :: istep_restart
  logical :: fexist
  integer :: np, np2
  real (8), allocatable :: phase(:,:)
  integer (8), allocatable :: gid(:)
  integer, parameter :: ict1=ptl_nphase+1
  integer, parameter :: ict2=ptl_nphase+ptl_nconst
  include "mpif.h"
  integer :: adios_read_method = ADIOS_READ_METHOD_BP
  integer*8 :: sel0=0, sel1=0, sel2=0, sel4=0
  integer*8 :: bb_start1(1), bb_count1(1), bb_start2(2), bb_count2(2)
  integer*8 :: bb_start4(4), bb_count4(4)

  integer :: istep_f0_restart, vpdata, mudata, ndata, idum
  real (8) :: time_f0_restart
  integer :: nnode

  integer :: inum, enum
  integer*8 :: ioff, eoff
  integer :: adios_comm, adios_comm_mype

  nnode=grid%nnode

  ! check existence of timestep.dat
  if(sml_mype==0) then
     inquire(file='timestep.dat',EXIST=fexist)
     if(.not. fexist) then
        print *, ' timestep.dat does not exist. sml_restart is forced to zero'
        sml_restart=.false.
     endif
  endif
  call MPI_BCAST(sml_restart, 1, MPI_LOGICAL, 0, sml_comm, ierr)
  call MPI_BARRIER(sml_comm,err) !SAK
  if(.not. sml_restart) then ! means fexist==.false. in proc 0
     return
  endif

  ! readng timestep.dat
  ! timestep.dat has changed to text format
  !if(sml_mype==0)then
  !  open(1,FILE="timestep.dat",FORM="unformatted")
  !  read(1)sml_istep
  !  close(1)
  !endif
  if(sml_mype==0) then
     open(unit=1,file='timestep.dat',status='old')
     read(1,*) sml_gstep
     read(1,*) sml_run_count
     close(1)
     sml_run_count=sml_run_count+1
  endif

  call MPI_BCAST(sml_gstep, 1, MPI_INTEGER, 0, sml_comm, ierr)
  call MPI_BCAST(sml_run_count, 1, MPI_INTEGER, 0, sml_comm, ierr)

  call MPI_BARRIER(sml_comm,err) !SAK

  if (adios_mfilesys) then
     write(filename,'("restart_dir",i1.1,"/xgc.restart.",i5.5,".bp")') sml_mfilesys_index+1, sml_gstep
     adios_comm = sml_mfilesys_comm
     adios_comm_mype = sml_mfilesys_mype
  else
     write(filename,'("restart_dir/xgc.restart.",i5.5,".bp")') sml_gstep
     adios_comm = sml_comm
     adios_comm_mype = sml_mype
  endif

! SAK DEBUGGING LINE   write(6,'("read ",a45,a15,i5)') filename,dirname,sml_mype

  call adios_read_open_file (buf_id, trim(filename), adios_read_method, adios_comm, err)
   if(err/=0) then
     print *, 'restart_read error: could not open file', filename
     stop
   endif

   call adios_schedule_read (buf_id, sel0, 'timestep', 0, 1, istep_restart, err)
   call adios_schedule_read (buf_id, sel0, 'time', 0, 1, sml_time, err)
   !call adios_schedule_read (buf_id, sel0, 'maxnum', 0, 1, spall(1)%maxnum, err)  ! do not read maxnum -- maxnum is memory size inialized with ptl_maxnum
   !call adios_schedule_read (buf_id, sel0, 'inum', 0, 1, spall(1)%num, err)
   call adios_schedule_read (buf_id, sel0, 'imaxgid', 0, 1, spall(1)%maxgid, err)
   if (sml_sheath_mode==1) then
      ! We know the size from sheath_init because we use the same limiter on restart!!!
      buf_size=psn%nwall * 8
      if (sml_mype==0) print *, 'Reading sheath potential: ',psn%nwall
      !ADIOS_READ_LBL(buf_id,'sheath_pot',psn%sheath_pot,buf_size,err)
      call adios_schedule_read(buf_id, sel0, 'sheath_pot', 0, 1, psn%sheath_pot, err)
   endif
   !call adios_perform_reads (buf_id, err)

   bb_start1(1) = adios_comm_mype
   bb_count1(1) = 1
   call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)
   call adios_schedule_read (buf_id, sel1, 'inum_arr', 0, 1, inum, err)
   call adios_schedule_read (buf_id, sel1, 'ioff_arr', 0, 1, ioff, err)
   call adios_perform_reads (buf_id, err)
   call adios_selection_delete(sel1) 
   spall(1)%num = inum
   
  !
  np=spall(1)%num
  allocate(phase(ict2,np),gid(np))

  bb_start1(1) = ioff
  bb_count1(1) = np
  call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)
  
  bb_start2(1) = 0
  bb_start2(2) = ioff
  bb_count2(1) = ict2
  bb_count2(2) = np
  call adios_selection_boundingbox (sel2, 2, bb_start2, bb_count2)

  call adios_schedule_read (buf_id, sel1, 'igid', 0, 1, gid, err)
  call adios_schedule_read (buf_id, sel2, 'iphase', 0, 1, phase, err)
  call adios_perform_reads (buf_id, err)
  call adios_selection_delete(sel1)
  call adios_selection_delete(sel2)

  if(np>ptl_maxnum) then
    np2=np + np*0.1  ! 10% margin
    print *, 'Not enough ptl_maxnum to load restart data'
    print *, 'array incrased from',ptl_maxnum, 'to', np2, 'mype=',sml_mype

    deallocate(spall(1)%ptl, spall(1)%phase0, spall(1)%tr_save, spall(1)%p_save, spall(1)%rhoi)

    spall(1)%maxnum=np2
    allocate( spall(1)%ptl(np2), spall(1)%phase0(ptl_nphase,np2))
    allocate( spall(1)%tr_save(np2),spall(1)%p_save(3,np2) )
    allocate(spall(1)%rhoi(np2))
   !stop
  endif

  do i=1, np
     spall(1)%ptl(i)%ph=phase(1:ptl_nphase,i)
     spall(1)%ptl(i)%ct=phase(ict1:ict2   ,i)
     spall(1)%ptl(i)%gid=gid(i)
  enddo
  deallocate(phase,gid)

  ! electron
  if(sml_electron_on)then
     !call adios_schedule_read (buf_id, sel0, 'enum', 0, 1, spall(0)%num, err)
     call adios_schedule_read (buf_id, sel0, 'emaxgid', 0, 1, spall(0)%maxgid, err)
     !call adios_perform_reads (buf_id, err)

     bb_start1(1) = adios_comm_mype
     bb_count1(1) = 1
     call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)
     call adios_schedule_read (buf_id, sel1, 'enum_arr', 0, 1, enum, err)
     call adios_schedule_read (buf_id, sel1, 'eoff_arr', 0, 1, eoff, err)
     call adios_perform_reads (buf_id, err)
     call adios_selection_delete(sel1) 
     spall(0)%num = enum

     !
     np=spall(0)%num
     allocate(phase(ict2,np),gid(np))

     bb_start1(1) = eoff
     bb_count1(1) = np
     call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)
     
     bb_start2(1) = 0
     bb_start2(2) = eoff
     bb_count2(1) = ict2
     bb_count2(2) = np
     call adios_selection_boundingbox (sel2, 2, bb_start2, bb_count2)
     
     call adios_schedule_read (buf_id, sel1, 'egid', 0, 1, gid, err)
     call adios_schedule_read (buf_id, sel2, 'ephase', 0, 1, phase, err)
     call adios_perform_reads (buf_id, err)
     call adios_selection_delete(sel1)
     call adios_selection_delete(sel2)

     if(np>ptl_e_maxnum) then
        np2=np + np*0.1  ! 10% margin
        print *, 'Not enough ptl_e_maxnum to load restart data'
        print *, 'array incrased from',ptl_e_maxnum, 'to', np2, 'mype=', sml_mype

        deallocate(spall(0)%ptl, spall(0)%phase0, spall(0)%tr_save, spall(0)%p_save,ptl_ephase_save)

        spall(0)%maxnum=np2
        allocate( spall(0)%ptl(np2), spall(0)%phase0(ptl_nphase,np2))
        allocate( spall(0)%tr_save(np2),spall(0)%p_save(3,np2) )
        allocate(ptl_ephase_save(np2))

     endif

     do i=1, np
        spall(0)%ptl(i)%ph=phase(1:ptl_nphase,i)
        spall(0)%ptl(i)%ct=phase(ict1:ict2   ,i)
        spall(0)%ptl(i)%gid=gid(i)
     enddo
     deallocate(phase,gid)
  endif

  call adios_read_close (buf_id, err)


  if(istep_restart/=sml_gstep) then
     print *, 'Warning: sml_gstep from timestep.dat and xgc.restart.XXXX.bp are different', &
          & sml_gstep, istep_restart,sml_mype, "-----"
  endif

  call validity_check(spall)

#ifdef XGC1_EM
  ! read hyb information
  if(sml_electron_hyb) then
     ! get size
     if(sml_plane_mype==0) then
        write(filename,'("restart_dir/xgc.restart_hyb.",i5.5,".bp")') sml_gstep
        write(dirname,'("/node_",i5.5)') sml_mype

        call adios_read_open_file (buf_id, trim(filename), adios_read_method, sml_intpl_comm, err)
        if(err/=0) then
           print *, 'restart_read error: could not open file', filename
           stop
        endif

        call adios_schedule_read (buf_id, sel1, trim(dirname)//'/eden_hyb', 0, 1, psn%eden_hyb, err)
        call adios_schedule_read (buf_id, sel1, trim(dirname)//'/a_par', 0, 1, psn%A_par, err)
        call adios_schedule_read (buf_id, sel1, trim(dirname)//'/u_e', 0, 1, psn%u_e, err)
        call adios_perform_reads (buf_id, err)
        call adios_read_close (buf_id, err)
     endif
     call MPI_BCAST(psn%eden_hyb, nnode, MPI_REAL8, 0, sml_plane_comm, ierr)
     call MPI_BCAST(psn%a_par, nnode, MPI_REAL8, 0, sml_plane_comm, ierr)
     call MPI_BCAST(psn%u_e, nnode, MPI_REAL8, 0, sml_plane_comm, ierr)
  endif
#endif

  ! read f0 information
  if(sml_f0_grid) then
     call check_point('read restart f0')

     if (adios_mfilesys) then
        write(filename,'("restart_dir",i1.1,"/xgc.restart.f0.",i5.5,".bp")') sml_mfilesys_index+1, sml_gstep
     else
        write(filename,'("restart_dir/xgc.restart.f0.",i5.5,".bp")') sml_gstep
     endif

     call adios_read_open_file (buf_id, trim(filename), adios_read_method, adios_comm, err)
     if(err/=0) then
       print *, 'restart_read error: could not open file', filename
       stop
     endif

     bb_start1(1) = adios_comm_mype
     bb_count1(1) = 1
     call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)

     call adios_schedule_read (buf_id, sel0, 'timestep', 0, 1, istep_f0_restart, err)
     call adios_schedule_read (buf_id, sel0, 'time', 0, 1, time_f0_restart, err)
     call adios_schedule_read (buf_id, sel0, 'vpdata', 0, 1, vpdata, err)
     call adios_schedule_read (buf_id, sel0, 'mudata', 0, 1, mudata, err)
     call adios_schedule_read (buf_id, sel1, 'inode1', 0, 1, f0_inode1, err)
     call adios_schedule_read (buf_id, sel1, 'inode2', 0, 1, f0_inode2, err)
     call adios_schedule_read (buf_id, sel1, 'imu1', 0, 1, f0_imu1, err)
     call adios_schedule_read (buf_id, sel1, 'imu2', 0, 1, f0_imu2, err)
     call adios_schedule_read (buf_id, sel1, 'ndata_arr', 0, 1, ndata, err)
     call adios_perform_reads (buf_id, err)
     call adios_selection_delete(sel1) 

     ! memory allocation


     ! check consistency of time data?
     ! time_f0_restart, istep_f0_restart?

     ! check consistency of data size variable
     if(f0_nvp/=(vpdata-1)/2) then
        print *, 'Wrong vpdata in restart file', f0_nvp, (vpdata-1)/2, vpdata, sml_mype
        stop
     endif
     if(f0_inode1 > f0_inode2) then
        print *, 'Wrong f0_inode1/2 in restart file', f0_inode1, f0_inode2, sml_mype
        stop
     endif
     if(f0_imu1 > f0_imu2 ) then
        print *, 'Wrong f0_inode1/2 in restart file', f0_imu1, f0_imu2, sml_mype
        stop
     endif
     if(f0_imu2-f0_imu1+1 /= mudata) then
        print *, 'Wrong f0_imu1/2 or mudata in restart file', f0_imu1, f0_imu2, mudata, sml_mype
     endif
     if(f0_inode2-f0_inode1+1 /= ndata) then
        print *, 'Wrong f0_inode1/2 or ndata in restart file', f0_inode1, f0_inode2, ndata, sml_mype
     endif
     !if(f0_inode1 <1 .or. f0_inode2 > grid%nnode ) --> need grid data
     !   print *, 'Wrong f0_inode1/2 range', f0_inode1, f0_inode2, sml_mype
     !endif

     if(f0_imu1 <0 .or. f0_imu2 > f0_nmu ) then
        print *, 'Wrong f0_imu1/2 range', f0_imu1, f0_imu2, sml_mype
     endif

     ! it is safer to call f0_initialize than just f0_mem_allocation
     call f0_initialize(grid,f0_inode1,f0_inode2,0,f0_nmu)
     
     ! read f0g data
     bb_start4(1) = 0
     bb_start4(2) = f0_inode1-1
     bb_start4(3) = 0
     bb_start4(4) = sml_intpl_mype
     bb_count4(1) = vpdata
     bb_count4(2) = ndata
     bb_count4(3) = mudata
     bb_count4(4) = 1
     call adios_selection_boundingbox (sel4, 4, bb_start4, bb_count4)
#ifndef F0_TOR_LINEAR
     call adios_schedule_read (buf_id, sel4, 'i_f0g', 0, 1, f0_f0g(:,:,:,1), err)
     call adios_perform_reads (buf_id, err)
#else
     call adios_schedule_read (buf_id, sel4, 'i_f0g', 0, 1, f0_f0g(:,:,:,0,1), err)
     call adios_perform_reads (buf_id, err)

     call adios_schedule_read (buf_id, sel4, 'i_f0g1', 0, 1, f0_f0g(:,:,:,1,1), err)
     call adios_perform_reads (buf_id, err)
#endif
     if(sml_electron_on) then
#ifndef F0_TOR_LINEAR
       call adios_schedule_read (buf_id, sel4, 'e_f0g', 0, 1, f0_f0g(:,:,:,0), err)
       call adios_perform_reads (buf_id, err)
#else
       call adios_schedule_read (buf_id, sel4, 'e_f0g', 0, 1, f0_f0g(:,:,:,0,0), err)
       call adios_perform_reads (buf_id, err)

       call adios_schedule_read (buf_id, sel4, 'e_f0g1', 0, 1, f0_f0g(:,:,:,1,0), err)
       call adios_perform_reads (buf_id, err)
#endif
     endif
     call adios_selection_delete(sel4)


#ifdef F0_CHARGE_N0
     if(sml_mype==0) then
      bb_start1(1)=0
      bb_count1(1)=grid%nnode
      call adios_selection_boundingbox (sel1, 1, bb_start1, bb_count1)
      if(sml_electron_on) then
        call adios_schedule_read(buf_id, sel1, 'f0_edensity_n0_add', 0, 1, &
            f0_density_n0_add(:,0),err)
      endif
      call adios_schedule_read(buf_id, sel1, 'f0_idensity_n0_add', 0, 1, &
            f0_density_n0_add(:,1),err)

     endif
     call adios_perform_reads(buf_id,err)

     call mpi_bcast(f0_density_n0_add, grid%nnode*(ptl_nsp-ptl_isp+1), &
          MPI_DOUBLE_PRECISION, 0, sml_comm,err)

#endif

     call adios_read_close (buf_id, err)
  endif


end subroutine restart_read
#endif

#ifdef ADIOS
subroutine check_adios_err(err,str)
  use sml_module
  implicit none
  integer :: err
  character (len=*) :: str

  return

  if(err/=0) then
     print *, 'ADIOS err:',str, 'pe=',sml_mype
     stop
  endif
end subroutine check_adios_err
#endif

subroutine validity_check(spall)
  use sml_module
  use ptl_module
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i, count, isp
  real (8) :: r, z, phi, mu, phimin, phimax, w

  count = 0
  phimin=sml_2pi_wedge_n/real(sml_nphi_total)* (sml_intpl_mype)
  phimax=sml_2pi_wedge_n/real(sml_nphi_total)* (sml_intpl_mype+1)

  do isp=1,1
     count=0

     if(spall(isp)%maxgid<=0) print *, 'invalid maxgid', spall(isp)%maxgid,sml_mype
     if(spall(isp)%maxnum<=0) print *, 'invalid maxnum', spall(isp)%maxnum,sml_mype
     if(spall(isp)%num<=0 ) print *, 'invalid num', spall(isp)%num, sml_mype
     if(ptl_nphase<=0) print *, 'invalid nphase', ptl_nphase, sml_mype


     do i=1, spall(isp)%num
        if(spall(isp)%ptl(i)%gid<=0) then
           if(count<10) print *, 'invalid gid:',spall(isp)%ptl(i)%gid,i, sml_mype
           count=count+1
        endif
        r  =spall(isp)%ptl(i)%ph(1)
        z  =spall(isp)%ptl(i)%ph(2)
        phi=spall(isp)%ptl(i)%ph(3)
        mu =spall(isp)%ptl(i)%ct(pim)
        w  =spall(isp)%ptl(i)%ct(piw0)
        ! space check
        if(.not. ( sml_bd_min_r <= r .and. r <= sml_bd_max_r)) then
           if(count<10) print *, 'invalid restart r:',r, sml_bd_min_r, sml_bd_max_r, i,sml_mype
           count=count+1
           if(spall(isp)%ptl(i)%gid>0 ) spall(isp)%ptl(i)%gid=-spall(isp)%ptl(i)%gid
        endif
        if(.not. ( sml_bd_min_z <= z .and. z <= sml_bd_max_z)) then
           if(count<10) print *, 'invalid restart z:',z, sml_bd_min_z, sml_bd_max_z, i,sml_mype
           count=count+1
           if(spall(isp)%ptl(i)%gid>0 ) spall(isp)%ptl(i)%gid=-spall(isp)%ptl(i)%gid
        endif
        if(.not. ( phimin <= phi .and. phi <= phimax)) then
           if(count<10) print *, 'invalid restart phi:',phi, phimin, phimax, i,sml_mype
           count=count+1
           if(spall(isp)%ptl(i)%gid>0 ) spall(isp)%ptl(i)%gid=-spall(isp)%ptl(i)%gid
        endif
        if(.not. ( mu >=0 )) then
           if(count<10) print *, 'invalid restart mu:',mu,i,sml_mype
           count=count+1
           if(spall(isp)%ptl(i)%gid>0 ) spall(isp)%ptl(i)%gid=-spall(isp)%ptl(i)%gid
        endif
        if(.not. (w >=0)) then
           if(count<10) print *, 'invalid restart w:',w,i,sml_mype
           count=count+1
           if(spall(isp)%ptl(i)%gid>0 ) spall(isp)%ptl(i)%gid=-spall(isp)%ptl(i)%gid
        endif
     enddo
     if(count>0) print *, 'Warning : invalid data in restart file', sml_mype, count
  enddo

  call memory_cleaning_simple(spall)
end subroutine validity_check

subroutine charging_test(istep,grid,psn)
  use sml_module
  use diag_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,istep,iphi
  real (kind=8) :: sum1(sml_nphi_total),tmp(sml_nphi_total)

  if(sml_mype==0) then
     do i=1, grid%nnode
        if(psn%idensity0(i)/=0D0) then
           write(1222,1000) sml_time/sml_tran,grid%x(1,i),grid%x(2,i), psn%idensity0(i)
        endif
     enddo
  endif
  write(1222,*) ' '

  do i=1, grid%nphi
     sum1(i+grid%iphi_offset)= sum(psn%idensity(:,i))
  enddo

  call my_mpi_reduce(sum1,tmp,grid%nphi)

  do i=1, sml_nphi_total
     write(1333,2000) sml_time/sml_tran,real(i),tmp(i)
  enddo

  write(1333,*) ' '

1000 format(4(e19.13,' '))
2000 format(3(e19.13,' '))

end subroutine charging_test




#ifdef ADIOS
subroutine background_edensity0_output(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  integer (8) :: buf_size, buf_id, total_size
  integer :: err

  ! dump 2d field diagnostics in adios bp format
  buf_size = 4 +  grid%nnode * 8 + 4 + grid%npsi00 * 8
  !          print *, 'adios writing diag_2d group, #bytes = ', buf_size
  ADIOS_OPEN(buf_id,'edensity0','xgc.edensity0.bp','w',sml_comm_null,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'nnode',grid%nnode,err)
  ADIOS_WRITE_LBL(buf_id,'edensity0_node',psn%edensity0,err)
  ADIOS_WRITE_LBL(buf_id,'npsi',grid%npsi00,err)
  ADIOS_WRITE_LBL(buf_id,'eden00_psi',psn%eden00_1d,err)
  ADIOS_CLOSE(buf_id,err)

end subroutine background_edensity0_output

subroutine background_edensity0_read(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use adios_read_mod
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer (8) :: buf_size, buf_id, total_size
  integer :: err, nnode, npsi
  integer :: adios_read_method = ADIOS_READ_METHOD_BP
  integer*8 :: sel1=0

  if(sml_mype==0) then
    call adios_read_open_file (buf_id, 'xgc.edensity0.bp', adios_read_method, MPI_COMM_SELF, err)
    if(err/=0) then
      print *, 'restart_read error: could not open file', 'xgc.edensity0.bp' 
      stop
    endif
    
    call adios_schedule_read (buf_id, sel1, 'nnode', 0, 1, nnode, err)
    call adios_schedule_read (buf_id, sel1, 'edensity0_node', 0, 1, psn%edensity0, err)
    call adios_schedule_read (buf_id, sel1, 'npsi', 0, 1, npsi, err)
    call adios_schedule_read (buf_id, sel1, 'eden00_psi', 0, 1, psn%eden00_1d, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)
    if(nnode/=grid%nnode) then
       print *, 'invalid nnode value in xgc.edensity0.bp', nnode
       stop
     endif
    if(npsi/=grid%npsi00) then
       print *, 'invalid npsi value in xgc.edensity0.bp', npsi
       stop
     endif
  endif

  ! broad case

  call mpi_bcast(psn%edensity0, grid%nnode, MPI_DOUBLE_PRECISION, 0, sml_comm,err)
  call mpi_bcast(psn%eden00_1d, grid%npsi00, MPI_DOUBLE_PRECISION, 0, sml_comm,err)

end subroutine

#endif

! initialize diag_neu routines
subroutine diag_neu_init
  use diag_module
  implicit none

  !diag_neu_npsi=neu_grid_mpsi
  !diag_neu_pin =
  !diag_neu_pout=

  diag_neu_dp_inv=(diag_neu_pout-diag_neu_pin)/real(diag_neu_npsi-1)
  diag_neu_dp_inv=1D0/diag_neu_dp_inv

  allocate(diag_neu_vol(diag_neu_npsi))
  allocate(diag_neu_port(diag_neu_nv,diag_neu_npsi,3))
  diag_neu_port=0D0

end subroutine diag_neu_init

subroutine diag_neu_ionization(ptli,psi,b,sp_type)
  use ptl_module
  use diag_module
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z
  implicit none
  type(ptl_type) , intent(in) :: ptli
  real (8), intent(in) :: psi, b
  integer, intent(in) :: sp_type
  !
  real (8) :: r, z, rho, mu, w, pn, wp
  integer :: ip
  real (8) :: v(diag_neu_nv)


  r=ptli%ph(pir)
  z=ptli%ph(piz)
  rho=ptli%ph(pirho)
  mu=ptli%ct(pim)
  w=ptli%ct(piw0)
  if( .not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z>eq_x2_z ) return ! Exit for private region

  pn=(psi-diag_neu_pin)*diag_neu_dp_inv
  ip=floor(pn)+1
  if(ip < 1 .or. diag_neu_npsi <= ip) return ! Exit for out of psi range
  wp=1D0 - pn + real(ip-1,8)

  v(1)= 1D0                                   ! full weight
  v(2)= ptl_c2_2m(sp_type)*(rho*B)**2 + mu*B  ! energy
  v(3)= ptl_c_m(sp_type)*rho*B                ! parallel momentum - without mass (v_parallel)
  v(4)= v(3)*r                                ! angular momentum - without B_phi / B
                                       !--> large q gives small poloidal flow--> ignorable?

  diag_neu_port(:,ip,1)=v(:)*w*wp      + diag_neu_port(:,ip,1)
  diag_neu_port(:,ip,1)=v(:)*w*(1D0-wp)+ diag_neu_port(:,ip,1)

end subroutine diag_neu_ionization


! diag_neu_elastic
subroutine diag_neu_elastic(ptli,rho2,mu2,psi,b,sp_type)
  use ptl_module
  use diag_module
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z
  implicit none
  type(ptl_type) , intent(in) :: ptli
  real (8), intent(in) :: rho2,mu2,psi, b
  integer, intent(in) :: sp_type
  !
  real (8) :: r, z, rho, mu, w, pn, wp
  integer :: ip
  real (8) :: v(diag_neu_nv)


  r=ptli%ph(pir)
  z=ptli%ph(piz)
  rho=ptli%ph(pirho)
  mu=ptli%ct(pim)
  w=ptli%ct(piw0)
  if( .not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z>eq_x2_z ) return ! Exit for private region

  pn=(psi-diag_neu_pin)*diag_neu_dp_inv
  ip=floor(pn)+1
  if(ip < 1 .or. diag_neu_npsi <= ip) return ! Exit for out of psi range
  wp=1D0 - pn + real(ip-1,8)

  v(1)= 1D0                                   ! full weight
  v(2)= ptl_c2_2m(sp_type)*B*B*(rho2*rho2-rho*rho) + (mu2-mu)*B  ! energy
  v(3)= ptl_c_m(sp_type)*(rho2-rho)*B                ! parallel momentum - without mass (v_parallel)
  v(4)= v(3)*r                                ! angular momentum - without B_phi / B
                                       !--> large q gives small poloidal flow--> ignorable?

  diag_neu_port(:,ip,2)=v(:)*w*wp      + diag_neu_port(:,ip,2)
  diag_neu_port(:,ip,2)=v(:)*w*(1D0-wp)+ diag_neu_port(:,ip,2)

end subroutine diag_neu_elastic


! diag_neu_cx
subroutine diag_neu_cx(ptli,rho2,mu2,psi,b,sp_type)
  use ptl_module
  use diag_module
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z
  implicit none
  type(ptl_type) , intent(in) :: ptli
  real (8), intent(in) :: rho2,mu2,psi, b
  integer, intent(in) :: sp_type
  !
  real (8) :: r, z, rho, mu, w, pn, wp
  integer :: ip
  real (8) :: v(diag_neu_nv)


  r=ptli%ph(pir)
  z=ptli%ph(piz)
  rho=ptli%ph(pirho)
  mu=ptli%ct(pim)
  w=ptli%ct(piw0)
  if( .not. is_rgn12(r,z,psi) .or. z < eq_x_z .or. z>eq_x2_z ) return ! Exit for private region

  pn=(psi-diag_neu_pin)*diag_neu_dp_inv
  ip=floor(pn)+1
  if(ip < 1 .or. diag_neu_npsi <= ip) return ! Exit for out of psi range
  wp=1D0 - pn + real(ip-1,8)

  v(1)= 1D0                                   ! full weight
  v(2)= ptl_c2_2m(sp_type)*B*B*(rho2*rho2-rho*rho) + (mu2-mu)*B  ! energy
  v(3)= ptl_c_m(sp_type)*(rho2-rho)*B                ! parallel momentum - without mass (v_parallel)
  v(4)= v(3)*r                                ! angular momentum - without B_phi / B
                                       !--> large q gives small poloidal flow--> ignorable?

  diag_neu_port(:,ip,3)=v(:)*w*wp      + diag_neu_port(:,ip,3)
  diag_neu_port(:,ip,3)=v(:)*w*(1D0-wp)+ diag_neu_port(:,ip,3)

end subroutine diag_neu_cx

#ifdef ADIOS
subroutine diag_neu_output(new_n0)
  use sml_module
  use diag_module
  use neu_module
  use eq_module
  implicit none
  real (8) :: new_n0
  real (8) :: tmp(diag_neu_nv,diag_neu_npsi,3)
  integer :: icase, ivar, ip
  real (8) :: pall(diag_neu_npsi), pnorm(diag_neu_npsi)
  real (8) :: tmp2(diag_neu_npsi)

  !ADIOS
  integer (8) :: buf_id, buf_size, total_size
  integer :: err
  character (len=64) :: var_names(diag_neu_nv)
  character (len=64) :: case_names(3)

  var_names(1) = 'dn'
  var_names(2) = 'denergy'
  var_names(3) = 'dv_para'
  var_names(4) = 'dv_para_r'

  case_names(1)='_ionize'
  case_names(2)='_elastic'
  case_names(3)='_cx'


  ! prepare output -- mpisum, normalization
  call my_mpi_reduce(diag_neu_port,tmp,diag_neu_nv*diag_neu_npsi*3)

  if(sml_mype==0) then

     do icase=1,3
        do ivar=1,diag_neu_nv
           tmp(ivar,:,icase)=tmp(ivar,:,icase)/diag_neu_vol
        enddo
     enddo

     do ip=1, diag_neu_npsi
        pall(ip)=diag_neu_pin + real(ip-1)/diag_neu_dp_inv
     enddo
     pnorm=pall/eq_x_psi


     ! write
     buf_size = 1000 + 8*diag_neu_nv*diag_neu_npsi*3  + 8*diag_neu_npsi
     if(sml_gstep/neu_col_period==1) then
        ADIOS_OPEN(buf_id,'diagnosis.neu','xgc.neudiag.bp','w',sml_comm_null,err)
     else
        ADIOS_OPEN(buf_id,'diagnosis.neu','xgc.neudiag.bp','a',sml_comm_null,err)
     endif
     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
     ADIOS_WRITE_LBL(buf_id,'samples',diag_neu_npsi,err)
     ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
     ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
     ADIOS_WRITE_LBL(buf_id,'tindex',sml_gstep/neu_col_period,err)
     ADIOS_WRITE_LBL(buf_id,'psi_mks',pall,err)
     ADIOS_WRITE_LBL(buf_id,'psi',pnorm,err)
     ADIOS_WRITE_LBL(buf_id,'volume_neu',diag_neu_vol,err)
     ADIOS_WRITE_LBL(buf_id,'neu_den0',new_n0,err)

     do icase=1,3
        do ivar=1,diag_neu_nv
           tmp2=tmp(ivar,:,icase)
           ADIOS_WRITE_LBL(buf_id,trim(var_names(ivar))//trim(case_names(icase)),tmp2,err)
        enddo
     enddo

     ! neutral 2d density, temperature, flow ??

     ADIOS_CLOSE(buf_id,err)
  endif

  ! initialize arrays
  diag_neu_port=0D0

end subroutine diag_neu_output
#endif

subroutine diag_heat_init
    use diag_module
    use sml_module
    use ptl_module
    implicit none

    ! setup parameters
    diag_heat_rmax= (/diag_heat_rmax1, diag_heat_rmax2, diag_heat_rmax3 /)
    diag_heat_rmin= (/diag_heat_rmin1, diag_heat_rmin2, diag_heat_rmin3 /)
    diag_heat_zmax= (/diag_heat_zmax1, diag_heat_zmax2, diag_heat_zmax3 /)
    diag_heat_zmin= (/diag_heat_zmin1, diag_heat_zmin2, diag_heat_zmin3 /)

    diag_heat_dr  = (diag_heat_rmax-diag_heat_rmin)/real(diag_heat_nr-1)
    diag_heat_dz  = (diag_heat_zmax-diag_heat_zmin)/real(diag_heat_nz-1)

    ! memory allocation
    allocate(diag_heat_pv(diag_heat_nvar,diag_heat_nr,diag_heat_nz,diag_heat_nsection,&
        ptl_isp:ptl_nsp,sml_nthreads))
    allocate(diag_heat_pv_psi(diag_heat_nvar,diag_heat_npsi,diag_heat_nsection,ptl_isp:ptl_nsp,sml_nthreads))
    diag_heat_pv=0D0
    diag_heat_pv_psi=0D0


    !finding pmax and pmin
    call find_pmax_pmin

    diag_heat_dp = (diag_heat_pmax-diag_heat_pmin)/real(diag_heat_npsi-1)

  contains
    !finding maximum/minimum psi of diagnostic box
    subroutine find_pmax_pmin
      implicit none
      integer :: i, j, k
      real (8) :: psi, pmin, pmax, r, z
      real (8), external :: psi_interpol
      
      ! for each section
      do i=1, diag_heat_nsection

         !initialize extream value
         pmin= 1D90
         pmax=-1D90

         do j=1, diag_heat_nr
            do k=1, diag_heat_nz
               r=diag_heat_dr(i)*real(j-1) + diag_heat_rmin(i)
               z=diag_heat_dz(i)*real(k-1) + diag_heat_zmin(i)
               
               psi=psi_interpol(r,z,0,0)
               
               pmin=min(pmin,psi)
               pmax=max(pmax,psi)
               
            enddo
         enddo
         diag_heat_pmin(i)=pmin
         diag_heat_pmax(i)=pmax
      enddo      

    end subroutine find_pmax_pmin

end subroutine diag_heat_init

subroutine diag_heat_port(w, pot, epara, eperp, ct, old_ph, new_ph,dphi, stype, ith)
  use sml_module
  use ptl_module
  use diag_module
  implicit none
  real (8), intent(in) :: w, pot, epara, eperp
  real (8), intent(in) :: ct(ptl_nconst), old_ph(ptl_nphase), new_ph(ptl_nphase), dphi
  integer, intent(in) :: stype, ith
  real (8) :: wp, r, z, psi, rn, zn, pn, wr,wz, ws, v(diag_heat_nvar)
  integer :: i, ir,iz, ip
  real (8), external :: psi_interpol
  real (8) :: x(2), phi, xff(2), phi_mid
  
  ! chracteristic r and z - use mean value for simplicity
!  r=(old_ph(pir)+new_ph(pir))*0.5D0
!  z=(old_ph(piz)+new_ph(piz))*0.5D0
   
  !get field following
  x=new_ph(1:2)
  phi=new_ph(3)
  phi_mid=(floor(phi/dphi) + 0.5D0) * dphi
  call field_following_pos2(x,phi,phi_mid,xff)
  r=xff(1)
  z=xff(2)
  


  wp=w*ct(piw0)
  
  v(1)=1D0
  v(2)=epara
  v(3)=eperp
  v(4)=pot
  v(5)=sqrt((old_ph(pir)-new_ph(pir))**2+(old_ph(piz)-old_ph(piz))**2)
  
  ! for all sections
  do i=1, diag_heat_nsection
     ! check range
     if(diag_heat_rmin(i) < r .and. r < diag_heat_rmax(i) &
          .and. diag_heat_zmin(i) < z .and. z < diag_heat_zmax(i)) then
        
        !r index
        rn=(r-diag_heat_rmin(i))/diag_heat_dr(i)
        ir=floor(rn)+1
        if(ir<1 .or. diag_heat_nr<= ir) cycle
        wr=1D0 - rn + real(ir-1,8)
        
        !z index
        zn=(z-diag_heat_zmin(i))/diag_heat_dz(i)
        iz=floor(zn)+1
        if(iz<1 .or. diag_heat_nz<= iz) cycle
        wz=1D0 - zn + real(iz-1,8)
                
        
        diag_heat_pv(:,ir  ,iz  ,i,stype,ith)=diag_heat_pv(:,ir  ,iz  ,i,stype,ith) + v*wp*wr      *wz
        diag_heat_pv(:,ir+1,iz  ,i,stype,ith)=diag_heat_pv(:,ir+1,iz  ,i,stype,ith) + v*wp*(1D0-wr)*wz
        diag_heat_pv(:,ir  ,iz+1,i,stype,ith)=diag_heat_pv(:,ir  ,iz+1,i,stype,ith) + v*wp*wr      *(1D0-wz)
        diag_heat_pv(:,ir+1,iz+1,i,stype,ith)=diag_heat_pv(:,ir+1,iz+1,i,stype,ith) + v*wp*(1D0-wr)*(1D0-wz)

        !psi 
        psi=psi_interpol(r,z,0,0)
        pn=(psi-diag_heat_pmin(i))/diag_heat_dp(i)
        ip=floor(pn)+1
        if(ip<1 .or. diag_heat_npsi<= ip) cycle
        ws=1D0 - pn + real(ip-1,8)
        
        diag_heat_pv_psi(:,ip  ,i,stype,ith)=diag_heat_pv_psi(:,ip  ,i,stype,ith) + v*wp*ws
        diag_heat_pv_psi(:,ip+1,i,stype,ith)=diag_heat_pv_psi(:,ip+1,i,stype,ith) + v*wp*(1D0-ws)

     endif
  enddo
  !
    
end subroutine diag_heat_port


subroutine diag_heat_output
  use sml_module
  use ptl_module
  use diag_module
  implicit none
  integer :: ith, n, n2
  real (8) :: pv_sum(diag_heat_nvar,diag_heat_nr,diag_heat_nz,diag_heat_nsection,ptl_isp:ptl_nsp)
  real (8) :: pv_psi_sum(diag_heat_nvar,diag_heat_npsi,diag_heat_nsection,ptl_isp:ptl_nsp)
  real (8) :: out1(diag_heat_nr, diag_heat_nz,diag_heat_nsection, diag_heat_nvar, ptl_isp:ptl_nsp)
  real (8) :: out2(diag_heat_npsi,diag_heat_nsection, diag_heat_nvar, ptl_isp:ptl_nsp)
  integer :: i, ir, iz,ip, isp, j, k
  real (8) :: rall(diag_heat_nr, diag_heat_nsection), zall(diag_heat_nz,diag_heat_nsection)
  real (8) :: psiall(diag_heat_npsi,diag_heat_nsection)
  ! ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  character (len=2) :: sp_name(0:6)=(/'e_', 'i_', 'i2', 'i3', 'i4', 'i5', 'i6'/)  ! max 6 ion species
  character (len=64):: var_names(diag_heat_nvar)
  
  var_names(1) ='number'
  var_names(2) ='para_energy'
  var_names(3) ='perp_energy'
  var_names(4) ='potential'
  var_names(5) ='distance'
  
  ! OpenMP summation
  if(sml_nthreads >= 2) then
  do ith=2, sml_nthreads
     diag_heat_pv(:,:,:,:,:,1)=diag_heat_pv(:,:,:,:,:,1)+diag_heat_pv(:,:,:,:,:,ith)
     diag_heat_pv_psi(:,:,:,:,1)=diag_heat_pv_psi(:,:,:,:,1)+diag_heat_pv_psi(:,:,:,:,ith)
  enddo
  endif
  ! MPI summation
  n=diag_heat_nvar*diag_heat_nr*diag_heat_nz*diag_heat_nsection*(ptl_nsp-ptl_isp+1)
  call my_mpi_reduce(diag_heat_pv(:,:,:,:,:,1), pv_sum, n)
  n2=diag_heat_nvar*diag_heat_npsi*diag_heat_nsection*(ptl_nsp-ptl_isp+1)
  call my_mpi_reduce(diag_heat_pv_psi(:,:,:,:,1), pv_psi_sum, n2)

  
  ! obtain mean values

#ifdef ADIOS
    ! ADIOS
  if(sml_mype==0) then
     buf_size=4 + 1000 +(diag_heat_nr+diag_heat_nz+diag_heat_npsi)*diag_heat_nsection*8+ n*8 + n2*8
     do i=1, diag_heat_nsection
        do ir=1,diag_heat_nr
           rall(ir,i)=diag_heat_rmin(i)+diag_heat_dr(i)*real(ir-1)
        enddo
        do iz=1,diag_heat_nz
           zall(iz,i)=diag_heat_zmin(i)+diag_heat_dz(i)*real(iz-1)
        enddo
        do ip=1,diag_heat_npsi
           psiall(ip,i)=diag_heat_pmin(i)+diag_heat_dp(i)*real(ip-1)
        enddo
     enddo
     
     ! open ADIOS
     if(sml_gstep/diag_1d_period==1) then
        ADIOS_OPEN(buf_id,'diagnosis.heat','xgc.heatdiag.bp','w',sml_comm_null,err)
     else
        ADIOS_OPEN(buf_id,'diagnosis.heat','xgc.heatdiag.bp','a',sml_comm_null,err)
     endif
     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
     
     ! write data
     ADIOS_WRITE_LBL(buf_id,'rsamples',diag_heat_nr,err)
     ADIOS_WRITE_LBL(buf_id,'zsamples',diag_heat_nz,err)
     ADIOS_WRITE_LBL(buf_id,'psamples',diag_heat_npsi,err)
     ADIOS_WRITE_LBL(buf_id,'section',diag_heat_nsection,err)
     ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
     ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
     ADIOS_WRITE_LBL(buf_id,'tindex',sml_gstep/diag_1d_period,err)
     ADIOS_WRITE_LBL(buf_id,'r',rall,err)
     ADIOS_WRITE_LBL(buf_id,'z',zall,err)
     ADIOS_WRITE_LBL(buf_id,'psi',psiall,err)
     
     do isp=ptl_isp, ptl_nsp
        do j=1,diag_heat_nvar
           do i=1, diag_heat_nsection

              !rz data
              do k=1, diag_heat_nz
                 out1(:,k,i,j,isp)=pv_sum(j,:,k,i,isp)
              enddo

              !psi data
              out2(:,i,j,isp)=pv_psi_sum(j,:,i,isp)                               

           enddo
           ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(j)),out1(:,:,:,j,isp),err)
           ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(j))//'_psi',out2(:,:,j,isp),err)

        enddo
     enddo
     
     !close ADIOS
     ADIOS_CLOSE(buf_id,err)
     
  endif
#endif

  
  ! clear diag_heat_pv
  diag_heat_pv=0D0
  diag_heat_pv_psi=0D0
end subroutine diag_heat_output



subroutine diag_sheath(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use diag_module, only : diag_1d_period
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  !
  integer*8 :: buf_id, buf_size, total_size
  integer :: err


  if(sml_mype==0 ) then
     buf_size= 4 + psn%nwall * 8 * 3 + 100 ! 100 is for safety
#ifdef ADIOS
     ! open ADIOS
     if(sml_gstep/diag_1d_period==1) then
        ADIOS_OPEN(buf_id,'diagnosis.sheath','xgc.sheathdiag.bp','w',sml_comm_null,err)
     else
        ADIOS_OPEN(buf_id,'diagnosis.sheath','xgc.sheathdiag.bp','a',sml_comm_null,err)
     endif
     ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
     ADIOS_WRITE_LBL(buf_id,'nwall',psn%nwall,err)
     ADIOS_WRITE_LBL(buf_id,'sheath_pot',psn%sheath_pot,err)
     ADIOS_WRITE_LBL(buf_id,'sheath_lost',psn%sheath_lost(:,1),err)
     ADIOS_WRITE_LBL(buf_id,'sheath_ilost',psn%sheath_ilost(:,1),err)
     
     ADIOS_CLOSE(buf_id,err)
#endif
  endif
end subroutine diag_sheath

!save change of densiy, energy and momentum
subroutine diag_f0_df_port1(idx,grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use f0_module
  use ptl_module
  use diag_module
  use eq_module, only : is_rgn12, eq_x_z, eq_x2_z
  implicit none
  integer, intent(in) :: idx
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: isp, imu, node, ivp, ip
  real (8) :: mu_vol, smu, mu, vol, en_th, vth, pn, wp, vp, en, mass
  real (8) :: var(diag_f0_df_npv1)

  if(.not. diag_f0_df_on) return   ! No diagnostic output



  ! start of calculation

  do isp=diag_1d_isp, diag_1d_nsp
    mass=ptl_mass(isp)
    do imu=f0_imu1, f0_imu2

      if(imu==0) then
        mu_vol=0.5D0
      elseif(imu==f0_nmu) then
        mu_vol=0.5D0
      else
        mu_vol=1D0
      endif

      if (imu==0) then
        smu=f0_dsmu/f0_mu0_factor
      else
        smu=imu*f0_dsmu  ! <--- This is v_perp/v_th for v_perp grid!!!!
      endif
      mu=imu*f0_dsmu ! <--- This is (v_perp/v_th)^2 for v_perp grid and mu_N for sqrt(mu) grid!!!!
      mu=mu*mu


      do node=f0_inode1, f0_inode2
        ! volume for data point
        vol=f0_grid_vol(node,isp)*mu_vol
        ! normalized v
        en_th=f0_t_ev(node,isp)*sml_ev2j
        vth=sqrt(en_th/mass)


        !get psi index
        pn=(grid%psi(node)-diag_1d_pin)*diag_1d_dp_inv
        ip=floor(pn)+1
        if(ip <1 .or. diag_1d_npsi <= ip) cycle
        wp=1D0 - pn + real(ip-1,8)


        do ivp=-f0_nvp, f0_nvp
          vp=ivp*f0_dvp

          en=0.5d0 * (vp*vp + mu) ! energy normalized by T

          if(imu==1 .and. ivp==0) then   ! zero imu has half volume
            var(1)=1D0/f0_grid_vol_vonly(node,isp)   ! volume calculation
          else
            var(1)=0D0
          endif
          var(2)=f0_f(ivp,node,imu,isp)              ! for density
          var(3)=f0_df0g(ivp,node,imu,isp)           ! for density change
          var(4)=f0_df0g(ivp,node,imu,isp)*en*en_th  ! for energy change
          var(5)=f0_df0g(ivp,node,imu,isp)*vp*vth*mass    ! for momentum change - mass will be multiplied later

          diag_f0_df_pv1(:,ip  ,isp,idx)=diag_f0_df_pv1(:,ip  ,isp,idx) + var(:)*vol*wp
          diag_f0_df_pv1(:,ip+1,isp,idx)=diag_f0_df_pv1(:,ip+1,isp,idx) + var(:)*vol*(1D0-wp)
        enddo ! ivp
      enddo ! node
    enddo ! imu
  enddo ! isp

end subroutine


subroutine diag_f0_df(istep,grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use ptl_module
  use diag_module
  use f0_module
  implicit none
  include 'mpif.h'
  integer, intent(in) :: istep
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (8) :: sum1(diag_f0_df_npv1,diag_1d_npsi,ptl_isp:ptl_nsp,diag_f0_df_nsource)
  integer :: np, nvar_total
  !for Output
  character (len=256) :: filename
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  real (8) :: tmp(diag_1d_npsi)
  real (8) :: pall(diag_1d_npsi), pnorm(diag_1d_npsi)
  character (len=2) :: sp_name(0:6)=(/'e_', 'i_', 'i2', 'i3', 'i4', 'i5', 'i6'/)  ! max 6 ion species
  character (len=64) ::  var_names(diag_f0_df_npv1)
  character (len=64) ::  source_names(diag_f0_df_nsource)
#endif
  integer :: isource, ip, isp, ivar

  if(.not. diag_f0_df_on) return   ! No diagnostic output

  if(mod(istep,sml_f_source_period)/=0) then
    if(sml_mype==0) print *, 'mismatch istep in diag_f0_df'
    return  ! only when f_source is called.
  endif

  np=diag_1d_npsi
  nvar_total=diag_f0_df_npv1
  !mpi reduce
  call my_mpi_reduce(diag_f0_df_pv1,sum1,diag_f0_df_npv1*np*(ptl_nsp-ptl_isp+1)*diag_f0_df_nsource)


#ifdef ADIOS
  if(sml_mype==0) then
    buf_size=4 + 1000 + 4*np*2 + 8*np*(ptl_nsp-ptl_isp+1)*(diag_f0_df_npv1)*diag_f0_df_nsource

    do ip=1,np
      pall(ip)=diag_1d_pin+diag_1d_dp*real(ip-1)
    enddo
    pnorm=pall/eq_x_psi

    !print *, 'adios writing diagnosis.fsource group, #bytes = ', buf_size, etag
    if(sml_gstep/sml_f_source_period==1) then
      ADIOS_OPEN(buf_id,'diagnosis.fsource','xgc.fsourcediag.bp','w',sml_comm_null,err)
    else
      ADIOS_OPEN(buf_id,'diagnosis.fsource','xgc.fsourcediag.bp','a',sml_comm_null,err)
    endif
    ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    ADIOS_WRITE_LBL(buf_id,'samples',np,err)
    ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
    ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
    ADIOS_WRITE_LBL(buf_id,'tindex',sml_gstep/sml_f_source_period,err)
    ADIOS_WRITE_LBL(buf_id,'psi_mks',pall,err)
    ADIOS_WRITE_LBL(buf_id,'psi',pnorm,err)


    var_names(1)='volume'
    var_names(2)='density'
    var_names(3)='density_change'
    var_names(4)='energy_change'
    var_names(5)='momentum_change'

    source_names(1)='_collision'
    source_names(2)='_heat_torque'
    source_names(3)='_neutral'
    source_names(4)='_radiation'

    do isource=1,diag_f0_df_nsource

      if(sum(sum1(1,:,1,isource))==0D0) cycle  ! this source type is not used

      do isp=ptl_isp,ptl_nsp
        do ivar=1,nvar_total

          !prepare data
          tmp=sum1(ivar,:,isp,isource)

          if(ivar==2) then ! make volume
            tmp=tmp/sum1(1,:,isp,isource)  ! weight sum / volume
          elseif(ivar>=3) then ! make relative change
            tmp=tmp/sum1(2,:,isp,isource)  ! [weichg change / weight sum] or [energy change / weight sum] or [momentum change / weight sum]
          endif

          !write adios file
          ADIOS_WRITE_LBL(buf_id,sp_name(isp)//trim(var_names(ivar))//trim(source_names(isource)),tmp,err)
        enddo
      enddo
    enddo
    !close adios file

    ADIOS_CLOSE(buf_id,err)
  endif
#endif

  call diag_f0_df_port_clear

end subroutine

subroutine diag_particle(grid, psn, spall)
  use sml_module
  use grid_class
  use psn_class
  use ptl_module
  use diag_module
#ifdef SC17DEMO
  !! jyc: temporary fix for SC17 demo
!  use coupling_core_edge
use new_coupling
#endif
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  !
  integer :: inum ! ion particle number to be saved per proc
  integer*8 :: ioffset, itotal
  integer :: enum ! electron particle number to be saved per proc
  integer*8 :: eoffset, etotal
  !
  real (4), allocatable :: iphase(:,:), ephase(:,:)  ! 32-bit real
  integer (8), allocatable :: igid(:), egid(:)
  integer, parameter :: ict1=ptl_nphase+1
  integer, parameter :: ict2=ptl_nphase+ptl_nconst

  integer :: isp

#ifdef ADIOS
  character (len=40) :: filename
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
#endif

  ! count numbers of particles to be saved for each processor and allocate memory------

  !ion
  isp=1
  call count_num(isp,inum)
  allocate(iphase(ict2,inum),igid(inum))


  !electron
  if(sml_electron_on) then
    isp=0
    call count_num(isp,enum)
    allocate(ephase(ict2,enum),egid(enum))
  endif


  ! get off set numbers in global array --------------------------------
  call get_offset(inum,ioffset,itotal)
  if(sml_electron_on) call get_offset(enum, eoffset,etotal)


  ! prepare ion data to be saved and write------------------------------

  !ion
  isp=1
  call prep_data(isp,inum,iphase,igid)

  !electron
  if(sml_electron_on) then
    isp=0
    call prep_data(isp,enum,ephase,egid)
  endif


  ! Write to file -----------------------------------------------
#ifdef ADIOS
  buf_size = 20 + 4 + 8   ! 20 for safety
  buf_size = buf_size + 8 + 8 + 8 + 4 + ict2*inum*4 + inum*8
  if(sml_electron_on) then
    buf_size = buf_size + 8 + 8 + 8 + 4 + ict2*enum*4 + enum*8
  endif

  if (adios_stage_particle) then
    write(filename,'("xgc.particle.bp")')
#ifdef SC17DEMO
    !! jyc: temporary fix for SC17 demo
    write(filename,'("xgc.particle.", A4, ".bp")') cce_my_side
#endif
  else
    write(filename,'("xgc.particle.",i5.5,".bp")') sml_gstep
  endif
  ADIOS_OPEN(buf_id,'particle',filename,'w',sml_comm,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'timestep',sml_gstep,err)
  ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)

  ADIOS_WRITE_LBL(buf_id,'inum_total',itotal,err)
  ADIOS_WRITE_LBL(buf_id,'inum',inum,err)
  ADIOS_WRITE_LBL(buf_id,'ioff',ioffset,err)
  ADIOS_WRITE_LBL(buf_id,'inphase',ict2,err)
  ADIOS_WRITE_LBL(buf_id,'igid',igid,err)
  ADIOS_WRITE_LBL(buf_id,'iphase',iphase,err)

  if(sml_electron_on) then
    ADIOS_WRITE_LBL(buf_id,'enum_total',etotal,err)
    ADIOS_WRITE_LBL(buf_id,'enum',enum,err)
    ADIOS_WRITE_LBL(buf_id,'eoff',eoffset,err)
    ADIOS_WRITE_LBL(buf_id,'enphase',ict2,err)
    ADIOS_WRITE_LBL(buf_id,'egid',egid,err)
    ADIOS_WRITE_LBL(buf_id,'ephase',ephase,err)
  endif

  ADIOS_CLOSE(buf_id,err)
#endif


!deallocate -----------------------------------------------------
deallocate(iphase,igid)
if(sml_electron_on) deallocate(ephase,egid)


contains
  subroutine count_num(isp,num)
    integer, intent(in) :: isp
    integer, intent(out) :: num
    integer :: i
    integer (8) :: gid1

    num=0 ! initialize

    do i=1, spall(isp)%num
      gid1=spall(isp)%ptl(i)%gid
      if(gid1>0) then
        if(mod(gid1-1,diag_particle_mod)==0) then
          num=num+1
        endif
      endif
    enddo
  end subroutine

  subroutine get_offset(num,offset,total)
    integer, intent(in) :: num
    integer (8), intent(out) :: offset, total
    integer :: err
    include 'mpif.h'
    !
    integer*8 :: num8, inum_all(sml_intpl_totalpe) ! inum array per plane
    integer*8 :: inumsum, inumsum_all(sml_plane_totalpe) ! sum of inums in each plane

    ! two step mpi_allgather to avoid mpi_allgather in com_world
    num8 = num
    call mpi_allgather(num8,1,MPI_INTEGER8,inum_all,1,MPI_INTEGER8,sml_intpl_comm,err)
    inumsum = sum(inum_all)
    call mpi_allgather(inumsum,1,MPI_INTEGER8,inumsum_all,1,MPI_INTEGER8,sml_plane_comm,err)
    offset = sum(inumsum_all(1:sml_plane_mype)) + sum(inum_all(1:sml_intpl_mype))  !mype has zero base
    total = sum(inumsum_all)
  end subroutine


  subroutine prep_data(isp,num,phase,gid)
    integer, intent(in) :: isp, num
    real (4), intent(out) :: phase(ict2,num)
    integer (8), intent(out) :: gid(num)
    integer :: i,n
    integer (8) :: gid1

    !
    n=0  !initialize

    do i=1, spall(isp)%num
      gid1=spall(isp)%ptl(i)%gid
      if(gid1>0) then
        if(mod(gid1-1,diag_particle_mod)==0) then
          n=n+1
          phase(1:ptl_nphase,n)=spall(isp)%ptl(i)%ph
          phase(ict1:ict2,n)   =spall(isp)%ptl(i)%ct
          gid(n)=gid1
        endif
      endif
    enddo

  end subroutine

end subroutine

