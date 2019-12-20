!! Diagnosis module
module diag_module_gpu
   use dimensions_mod_gpu, only : nthreads_dim
   use ptl_module_gpu,only : ptl_nsp_max

   use util_mod_gpu 
   use precision_mod_gpu
   use diag_module, only :  &
      diag_max_sp_num_host => diag_max_sp_num, &
      diag_1d_npv1_host => diag_1d_npv1, &
      diag_heat_nvar_host => diag_heat_nvar
   implicit none

   integer, parameter :: diag_max_sp_num=diag_max_sp_num_host  ! two species

  ! tracer 
!  integer :: diag_tracer_period !! Period for tracer
!  integer :: diag_tracer_n     !! Particle index for tracer routine
!  integer :: diag_tracer_sp     ! Species index for tracer 

  ! for 1D
   integer, parameter :: diag_1d_npv1=diag_1d_npv1_host

!  Integer parameters for heat load diagnosis
   integer, parameter :: diag_heat_nvar=diag_heat_nvar_host

   integer :: diag_heat_nr, diag_heat_nz, diag_heat_npsi, diag_heat_nsection 

!  logical  :: diag_1d_on
   logical  :: diag_eflux_on
!  integer  :: diag_1d_period
   integer  :: diag_1d_npsi  !
   integer  :: diag_1d_ne
   real (kind=work_p) :: diag_1d_pin 
!  real (8) :: diag_1d_pout
!  real (8) :: diag_1d_dp 
   real (kind=work_p) :: diag_1d_dp_inv
   real (kind=work_p) :: diag_1d_emin
   real (kind=work_p) :: diag_1d_emax
!  real (8),allocatable :: diag_1d_vol(:)
!  integer :: diag_1d_isp,diag_1d_nsp  ! same as ptl_isp/nsp
!  logical  :: diag_tavg_on
!  
!   integer, parameter :: ptl_isp = 0
   integer, parameter :: diag_1d_isp = 0

   ! ------------------------------------------
   ! value for diag_1d_npsi is set in setup.F90
   ! ------------------------------------------
   integer, parameter :: diag_1d_npsi_dim = 80 
   integer, parameter :: diag_1d_nsp = ptl_nsp_max

   integer, parameter :: np = diag_1d_npsi_dim
   integer, parameter :: nsp = diag_1d_nsp
   integer, parameter :: isp = diag_1d_isp
   logical :: diag_heat_on 

   private :: diag_1d_isp, diag_1d_npsi_dim, diag_1d_nsp
   private :: np,nsp,isp

#ifdef ORIGINAL
   integer, parameter :: diag_1d_dim4 = nthreads_dim
#else
   integer, parameter :: diag_1d_dim4 = 16
#endif
   real(kind=work_p), allocatable :: diag_1d_f_pv1(:,:,:,:)
   real(kind=work_p), allocatable :: diag_1d_df_pv1(:,:,:,:)
   real(kind=work_p), allocatable :: diag_1d_eflux_pv(:,:,:,:,:)
   real(kind=work_p), allocatable :: diag_2d_dflux_pv(:,:,:,:,:)

! for heat load diagnosis
   real(kind=work_p), allocatable :: diag_heat_pv(:,:,:,:,:,:), diag_heat_pv_psi(:,:,:,:,:)
   real(kind=work_p), allocatable :: diag_heat_rmax(:), diag_heat_rmin(:), diag_heat_zmax(:), diag_heat_zmin(:), &
                                     diag_heat_dr(:), diag_heat_dz(:), diag_heat_pmax(:), diag_heat_pmin(:), &
                                     diag_heat_dp(:)
!  real (8), allocatable :: diag_1d_tavg_f_pv1(:,:,:,:)
!  real (8), allocatable :: diag_1d_tavg_df_pv1(:,:,:,:)
!
!  logical :: diag_3d_on
!  integer :: diag_3d_period  

   attributes(constant) :: diag_1d_npsi
   attributes(constant) :: diag_1d_pin, diag_1d_dp_inv
   attributes(device) :: diag_1d_f_pv1, diag_1d_df_pv1
   attributes(device) :: diag_heat_nr, diag_heat_nz, diag_heat_npsi, diag_heat_nsection , diag_heat_on
!   attributes(device) :: diag_heat_rmax1, diag_heat_rmax2, diag_heat_rmax3, diag_heat_rmin1, diag_heat_rmin2, &
!                         diag_heat_rmin2, diag_heat_rmin3, diag_heat_zmax1, diag_heat_zmax2, diag_heat_zmax3, &
!                         diag_heat_zmin1, diag_heat_zmin2, diag_heat_zmin3
   attributes(device) :: diag_heat_rmax, diag_heat_rmin, diag_heat_zmax, diag_heat_zmin, diag_heat_dr, &
                         diag_heat_dz, diag_heat_pmax, diag_heat_pmin, diag_heat_dp
   attributes(device) :: diag_heat_pv, diag_heat_pv_psi
   attributes(device) :: diag_1d_eflux_pv, diag_2d_dflux_pv
   ! These get assigned (non-automatically) in a kernel. Why are they module variables?
   attributes(device) :: diag_1d_emin, diag_1d_emax
   attributes(constant) :: diag_1d_ne, diag_eflux_on
   contains



   attributes(host) &
   subroutine update_device_diag()
   use assert_mod
   use sml_module, only : sml_mype

   use diag_module, only : &
      diag_1d_npsi_host => diag_1d_npsi, &
      diag_1d_dp_inv_host => diag_1d_dp_inv, &
      diag_1d_pin_host => diag_1d_pin, &
      diag_1d_df_pv1_host => diag_1d_df_pv1, &
      diag_1d_f_pv1_host => diag_1d_f_pv1, &
      diag_heat_nr_host => diag_heat_nr, &
      diag_heat_nz_host => diag_heat_nz, &
      diag_heat_npsi_host => diag_heat_npsi, &
      diag_heat_nsection_host => diag_heat_nsection, &
      diag_heat_on_host => diag_heat_on, &
!      diag_heat_rmax1_host => diag_heat_rmax1, &
!      diag_heat_rmax2_host => diag_heat_rmax2, &
!      diag_heat_rmax3_host => diag_heat_rmax3, &
!      diag_heat_rmin1_host => diag_heat_rmin1, &
!      diag_heat_rmin2_host => diag_heat_rmin2, &
!      diag_heat_rmin3_host => diag_heat_rmin3, &
!      diag_heat_zmax1_host => diag_heat_zmax1, &
!      diag_heat_zmax2_host => diag_heat_zmax2, &
!      diag_heat_zmax3_host => diag_heat_zmax3, &
!      diag_heat_zmin1_host => diag_heat_zmin1, &
!      diag_heat_zmin2_host => diag_heat_zmin2, &
!      diag_heat_zmin3_host => diag_heat_zmin3, &
      diag_heat_rmax_host => diag_heat_rmax, &
      diag_heat_rmin_host => diag_heat_rmin, &
      diag_heat_zmax_host => diag_heat_zmax, &
      diag_heat_zmin_host => diag_heat_zmin, &
      diag_heat_dr_host => diag_heat_dr, &
      diag_heat_dz_host => diag_heat_dz, &
      diag_heat_pmax_host => diag_heat_pmax, &
      diag_heat_pmin_host => diag_heat_pmin, &
      diag_heat_dp_host => diag_heat_dp, &
      diag_heat_pv_host => diag_heat_pv, &
      diag_heat_pv_psi_host => diag_heat_pv_psi, &
      diag_eflux_on_host => diag_eflux_on, &
      diag_1d_ne_host => diag_1d_ne, &
      diag_1d_emin_host => diag_1d_emin, &
      diag_1d_emax_host => diag_1d_emax, &
      diag_1d_eflux_pv_host => diag_1d_eflux_pv, &
      diag_2d_dflux_pv_host => diag_2d_dflux_pv

   integer, parameter :: idebug = 2

   integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb,ub
   integer :: ierr
   logical :: need_diag_1d_f_pv1
   logical :: need_diag_1d_df_pv1
   logical :: need_diag_heat_pv
   logical :: need_diag_heat_pv_psi
   logical :: need_diag_1d_eflux_pv
   logical :: need_diag_2d_dflux_pv

   diag_1d_npsi = diag_1d_npsi_host
   diag_1d_dp_inv = diag_1d_dp_inv_host
   diag_1d_pin = diag_1d_pin_host

   diag_heat_npsi = diag_heat_npsi_host
   diag_heat_nr = diag_heat_nr_host
   diag_heat_nz = diag_heat_nz_host
   diag_heat_nsection = diag_heat_nsection_host
   diag_heat_on = diag_heat_on_host

!   diag_heat_rmax1 = diag_heat_rmax1_host
!   diag_heat_rmax2 = diag_heat_rmax2_host
!   diag_heat_rmax3 = diag_heat_rmax3_host
!   diag_heat_rmin1 = diag_heat_rmin1_host
!   diag_heat_rmin2 = diag_heat_rmin2_host
!   diag_heat_rmin3 = diag_heat_rmin3_host
!   diag_heat_zmax1 = diag_heat_zmax1_host
!   diag_heat_zmax2 = diag_heat_zmax2_host
!   diag_heat_zmax3 = diag_heat_zmax3_host
!   diag_heat_zmin1 = diag_heat_zmin1_host
!   diag_heat_zmin2 = diag_heat_zmin2_host
!   diag_heat_zmin3 = diag_heat_zmin3_host
   diag_eflux_on =diag_eflux_on_host
   diag_1d_ne =diag_1d_ne_host
   diag_1d_emin =diag_1d_emin_host
   diag_1d_emax =diag_1d_emax_host

   need_diag_1d_f_pv1 = allocated(diag_1d_f_pv1_host)
   if (need_diag_1d_f_pv1) then

!       ------------------
!       clear previous storage
!       ------------------
     if (allocated(diag_1d_f_pv1)) then
      deallocate( diag_1d_f_pv1, stat=ierr)
      call assert(ierr.eq.0,'dealloc(diag_1d_f_pv1),ierr',ierr)
     endif

   lb1 = lbound(diag_1d_f_pv1_host,1)
   lb2 = lbound(diag_1d_f_pv1_host,2)
   lb3 = lbound(diag_1d_f_pv1_host,3)

   ub1 = ubound(diag_1d_f_pv1_host,1)
   ub2 = ubound(diag_1d_f_pv1_host,2)
   ub3 = ubound(diag_1d_f_pv1_host,3)

    
     allocate( diag_1d_f_pv1(lb1:ub1,lb2:ub2,lb3:ub3,diag_1d_dim4),     &
             stat=ierr)
     call assert(ierr.eq.0,'alloc(diag_1d_f_pv1),diag_1d_dim4',         &
             diag_1d_dim4)

   endif

   need_diag_1d_df_pv1 = allocated(diag_1d_df_pv1_host)
   if (need_diag_1d_df_pv1) then
!       ------------------
!       clear previous storage
!       ------------------
     if (allocated(diag_1d_df_pv1)) then
       deallocate( diag_1d_df_pv1,stat=ierr)
       call assert(ierr.eq.0,'dealloc(diag_1d_df_pv1,ierr=',ierr)
     endif

   lb1 = lbound(diag_1d_df_pv1_host,1)
   lb2 = lbound(diag_1d_df_pv1_host,2)
   lb3 = lbound(diag_1d_df_pv1_host,3)

   ub1 = ubound(diag_1d_df_pv1_host,1)
   ub2 = ubound(diag_1d_df_pv1_host,2)
   ub3 = ubound(diag_1d_df_pv1_host,3)


     allocate( diag_1d_df_pv1(lb1:ub1,lb2:ub2,lb3:ub3,diag_1d_dim4), &
             stat=ierr)
     call assert(ierr.eq.0,'alloc(diag_1d_df_pv1),diag_1d_dim4',  &
             diag_1d_dim4)
   endif


   if (allocated(diag_1d_f_pv1)) then
      diag_1d_f_pv1 = 0
   endif

   if (allocated(diag_1d_df_pv1)) then
      diag_1d_df_pv1 = 0
   endif
   
   if (diag_heat_on_host) then
     if (allocated(diag_heat_rmax)) then
        deallocate(diag_heat_rmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_rmax),ierr=',ierr) 
     endif
     lb = lbound(diag_heat_rmax_host,1)
     ub = ubound(diag_heat_rmax_host,1)
     allocate( diag_heat_rmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_rmax),ierr=',ierr) 
     diag_heat_rmax=diag_heat_rmax_host


     if (allocated(diag_heat_rmin)) then      
        deallocate(diag_heat_rmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_rmin),ierr=',ierr)
     endif        
     lb = lbound(diag_heat_rmin_host,1)
     ub = ubound(diag_heat_rmin_host,1)
     allocate( diag_heat_rmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_rmin),ierr=',ierr)
     diag_heat_rmin=diag_heat_rmin_host


     if (allocated(diag_heat_zmax)) then
        deallocate(diag_heat_zmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_zmax),ierr=',ierr)
     endif
     lb = lbound(diag_heat_zmax_host,1)
     ub = ubound(diag_heat_zmax_host,1)
     allocate( diag_heat_zmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_zmax),ierr=',ierr)
     diag_heat_zmax=diag_heat_zmax_host


     if (allocated(diag_heat_zmin)) then
        deallocate(diag_heat_zmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_zmin),ierr=',ierr)
     endif
     lb = lbound(diag_heat_zmin_host,1)
     ub = ubound(diag_heat_zmin_host,1)
     allocate( diag_heat_zmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_zmin),ierr=',ierr)
     diag_heat_zmin=diag_heat_zmin_host


     if (allocated(diag_heat_dr)) then
        deallocate(diag_heat_dr, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dr),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dr_host,1)
     ub = ubound(diag_heat_dr_host,1)
     allocate( diag_heat_dr(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dr),ierr=',ierr)
     diag_heat_dr=diag_heat_dr_host


     if (allocated(diag_heat_dz)) then
        deallocate(diag_heat_dz, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dz),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dz_host,1)
     ub = ubound(diag_heat_dz_host,1)
     allocate( diag_heat_dz(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dz),ierr=',ierr)
    diag_heat_dz=diag_heat_dz_host


     if (allocated(diag_heat_pmax)) then
        deallocate(diag_heat_pmax, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_pmax),ierr=',ierr)
     endif
     lb = lbound(diag_heat_pmax_host,1)
     ub = ubound(diag_heat_pmax_host,1)
     allocate( diag_heat_pmax(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_pmax),ierr=',ierr)
     diag_heat_pmax=diag_heat_pmax_host


     if (allocated(diag_heat_pmin)) then
        deallocate(diag_heat_pmin, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_pmin),ierr=',ierr)
     endif
     lb = lbound(diag_heat_pmin_host,1)
     ub = ubound(diag_heat_pmin_host,1)
     allocate( diag_heat_pmin(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_pmin),ierr=',ierr)
     diag_heat_pmin=diag_heat_pmin_host


     if (allocated(diag_heat_dp)) then
        deallocate(diag_heat_dp, stat=ierr)
        call assert(ierr.eq.0,'dealloc(diag_heat_dp),ierr=',ierr)
     endif
     lb = lbound(diag_heat_dp_host,1)
     ub = ubound(diag_heat_dp_host,1)
     allocate( diag_heat_dp(lb:ub), stat=ierr )
     call assert(ierr.eq.0, 'alloc(diag_heat_dp),ierr=',ierr)
     diag_heat_dp=diag_heat_dp_host


     need_diag_heat_pv = allocated(diag_heat_pv_host)
     if (need_diag_heat_pv) then
        if (allocated(diag_heat_pv)) then
           deallocate( diag_heat_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_heat_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_heat_pv_host,1)
        lb2 = lbound(diag_heat_pv_host,2)
        lb3 = lbound(diag_heat_pv_host,3)
        lb4 = lbound(diag_heat_pv_host,4)
        lb = lbound(diag_heat_pv_host,5)

        ub1 = ubound(diag_heat_pv_host,1)
        ub2 = ubound(diag_heat_pv_host,2)
        ub3 = ubound(diag_heat_pv_host,3)
        ub4 = ubound(diag_heat_pv_host,4)
        ub = ubound(diag_heat_pv_host,5)
        allocate( diag_heat_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb:ub,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_heat_pv),diag_1d_dim4',  &
             diag_1d_dim4,size(diag_heat_pv_host))
     endif


     need_diag_heat_pv_psi = allocated(diag_heat_pv_psi_host)
     if (need_diag_heat_pv_psi) then
        if (allocated(diag_heat_pv_psi)) then
           deallocate( diag_heat_pv_psi,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_heat_pv_psi,ierr=',ierr)
        endif
        lb1 = lbound(diag_heat_pv_psi_host,1)
        lb2 = lbound(diag_heat_pv_psi_host,2)
        lb3 = lbound(diag_heat_pv_psi_host,3)
        lb4 = lbound(diag_heat_pv_psi_host,4)

        ub1 = ubound(diag_heat_pv_psi_host,1)
        ub2 = ubound(diag_heat_pv_psi_host,2)
        ub3 = ubound(diag_heat_pv_psi_host,3)
        ub4 = ubound(diag_heat_pv_psi_host,4)
        allocate( diag_heat_pv_psi(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_heat_pv_psi),diag_1d_dim4',  &
             diag_1d_dim4,size(diag_heat_pv_psi_host))
     endif

     if (allocated(diag_heat_pv)) then
        diag_heat_pv = 0_work_p
     endif

     if (allocated(diag_heat_pv_psi)) then
        diag_heat_pv_psi = 0_work_p
     endif

   endif

   if(diag_eflux_on_host) then

     need_diag_1d_eflux_pv = allocated(diag_1d_eflux_pv_host)
  
     if (need_diag_1d_eflux_pv) then
        if (allocated(diag_1d_eflux_pv)) then
           deallocate( diag_1d_eflux_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_1d_eflux_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_1d_eflux_pv_host,1)
        lb2 = lbound(diag_1d_eflux_pv_host,2)
        lb3 = lbound(diag_1d_eflux_pv_host,3)
        lb4 = lbound(diag_1d_eflux_pv_host,4)

        ub1 = ubound(diag_1d_eflux_pv_host,1)
        ub2 = ubound(diag_1d_eflux_pv_host,2)
        ub3 = ubound(diag_1d_eflux_pv_host,3)
        ub4 = ubound(diag_1d_eflux_pv_host,4)
        allocate( diag_1d_eflux_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_1d_eflux_pv),diag_1d_dim4',  &
             diag_1d_dim4,size(diag_1d_eflux_pv_host))
     endif

     need_diag_2d_dflux_pv = allocated(diag_2d_dflux_pv_host)
     if (need_diag_2d_dflux_pv) then
        if (allocated(diag_2d_dflux_pv)) then
           deallocate( diag_2d_dflux_pv,stat=ierr)
           call assert(ierr.eq.0,'dealloc(diag_2d_dflux_pv,ierr=',ierr)
        endif
        lb1 = lbound(diag_2d_dflux_pv_host,1)
        lb2 = lbound(diag_2d_dflux_pv_host,2)
        lb3 = lbound(diag_2d_dflux_pv_host,3)
        lb4 = lbound(diag_2d_dflux_pv_host,4)

        ub1 = ubound(diag_2d_dflux_pv_host,1)
        ub2 = ubound(diag_2d_dflux_pv_host,2)
        ub3 = ubound(diag_2d_dflux_pv_host,3)
        ub4 = ubound(diag_2d_dflux_pv_host,4)
        allocate( diag_2d_dflux_pv(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,diag_1d_dim4), &
                  stat=ierr)
        call assert(ierr.eq.0,'alloc(diag_2d_dflux_pv),diag_1d_dim4',  &
             diag_1d_dim4,size(diag_2d_dflux_pv))
     endif

     if (allocated(diag_1d_eflux_pv)) then
        diag_1d_eflux_pv = 0_work_p
     endif

     if (allocated(diag_2d_dflux_pv)) then
        diag_2d_dflux_pv = 0_work_p
     endif

   endif ! diag_eflux_on_host

   return
   end subroutine update_device_diag



   attributes(host) &
   subroutine update_host_diag()
   use sml_module, only : sml_mype
   use precision_mod_gpu
   use diag_module, only : &
      diag_1d_npsi_host => diag_1d_npsi, &
      diag_1d_dp_inv_host => diag_1d_dp_inv, &
      diag_1d_pin_host => diag_1d_pin, &
      diag_1d_df_pv1_host => diag_1d_df_pv1, &
      diag_1d_f_pv1_host => diag_1d_f_pv1, &
      diag_heat_pv_host => diag_heat_pv, &
      diag_heat_pv_psi_host => diag_heat_pv_psi, &
      diag_1d_eflux_pv_host => diag_1d_eflux_pv, &
      diag_2d_dflux_pv_host => diag_2d_dflux_pv, &
      diag_eflux_on_host => diag_eflux_on, &
      diag_heat_on_host => diag_heat_on

   integer, parameter :: idebug = 2

   real (kind=work_p), allocatable, dimension(:,:,:,:) :: p_diag_1d_f_pv1
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: p_diag_1d_df_pv1

   real (kind=work_p), allocatable, dimension(:,:,:,:,:,:) :: p_diag_heat_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_heat_pv_psi

   real (kind=work_p), allocatable, dimension(:,:,:) :: dsum_f, dsum_df
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: dsum_heat_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_heat_psi
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_eflux_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:) :: dsum_dflux_pv

   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_1d_eflux_pv
   real (kind=work_p), allocatable, dimension(:,:,:,:,:) :: p_diag_2d_dflux_pv

   integer :: i1,i2,i3,i4,i5,i6
   integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4
   integer :: llb1,uub1,llb2,uub2,llb3,uub3,llb4,uub4

   integer :: lb1_df,ub1_df,lb2_df,ub2_df,lb3_df,ub3_df,lb4_df,ub4_df
   integer :: llb1_df,uub1_df,llb2_df,uub2_df,llb3_df,uub3_df,llb4_df,uub4_df

   integer :: lb1_pv,ub1_pv,lb2_pv,ub2_pv,lb3_pv,ub3_pv,lb4_pv,ub4_pv,lb5_pv,ub5_pv,lb6_pv,ub6_pv
   integer :: lb1_psi,ub1_psi,lb2_psi,ub2_psi,lb3_psi,ub3_psi,lb4_psi,ub4_psi,lb5_psi,ub5_psi

   integer :: lb1_ef,ub1_ef,lb2_ef,ub2_ef,lb3_ef,ub3_ef,lb4_ef,ub4_ef,lb5_ef,ub5_ef
   integer :: lb1_2df,ub1_2df,lb2_2df,ub2_2df,lb3_2df,ub3_2df,lb4_2df,ub4_2df,lb5_2df,ub5_2df

   integer :: ierr

   logical,parameter :: use_dsum = .true.
   logical,parameter :: use_p_diag = .true.
   logical :: need_diag_1d_df_pv1
   logical :: need_diag_heat_pv
   logical :: need_diag_heat_pv_psi
   logical :: need_diag_1d_eflux_pv
   logical :: need_diag_2d_dflux_pv
   integer :: n1,n2,n3,n4,n5,n6
   integer :: l1,u1,l2,u2,l3,u3,l4,u4,l5,u5,l6,u6

   diag_1d_npsi_host = diag_1d_npsi
   diag_1d_dp_inv_host = diag_1d_dp_inv
   diag_1d_pin_host = diag_1d_pin

   need_diag_1d_df_pv1 = allocated(diag_1d_df_pv1_host)
   need_diag_heat_pv = allocated(diag_heat_pv_host)   
   need_diag_heat_pv_psi = allocated(diag_heat_pv_psi_host)
   need_diag_1d_eflux_pv = allocated(diag_1d_eflux_pv_host)
   need_diag_2d_dflux_pv = allocated(diag_2d_dflux_pv_host)

   lb1 = lbound(diag_1d_f_pv1_host,1)
   lb2 = lbound(diag_1d_f_pv1_host,2)
   lb3 = lbound(diag_1d_f_pv1_host,3)
   lb4 = lbound(diag_1d_f_pv1_host,4)

   ub1 = ubound(diag_1d_f_pv1_host,1)
   ub2 = ubound(diag_1d_f_pv1_host,2)
   ub3 = ubound(diag_1d_f_pv1_host,3)
   ub4 = ubound(diag_1d_f_pv1_host,4)


   lb1_df = lbound(diag_1d_df_pv1_host,1)
   lb2_df = lbound(diag_1d_df_pv1_host,2)
   lb3_df = lbound(diag_1d_df_pv1_host,3)
   lb4_df = lbound(diag_1d_df_pv1_host,4)

   ub1_df = ubound(diag_1d_df_pv1_host,1)
   ub2_df = ubound(diag_1d_df_pv1_host,2)
   ub3_df = ubound(diag_1d_df_pv1_host,3)
   ub4_df = ubound(diag_1d_df_pv1_host,4)

   llb1 = lbound(diag_1d_f_pv1,1)
   llb2 = lbound(diag_1d_f_pv1,2)
   llb3 = lbound(diag_1d_f_pv1,3)
   llb4 = lbound(diag_1d_f_pv1,4)

   uub1 = ubound(diag_1d_f_pv1,1)
   uub2 = ubound(diag_1d_f_pv1,2)
   uub3 = ubound(diag_1d_f_pv1,3)
   uub4 = ubound(diag_1d_f_pv1,4)

   if (idebug >= 2) then
   call assert( llb1 <= lb1,'diag_module_gpu: invalid llb1 ',lb1,llb1)
   call assert( llb2 <= lb2,'diag_module_gpu: invalid llb2 ',lb2,llb2)
   call assert( llb3 <= lb3,'diag_module_gpu: invalid llb3 ',lb3,llb3)

   call assert( ub1 <= uub1,'diag_module_gpu: invalid uub1 ',ub1,uub1)
   call assert( ub2 <= uub2,'diag_module_gpu: invalid uub2 ',ub2,uub2)
   call assert( ub3 <= uub3,'diag_module_gpu: invalid uub3 ',ub3,uub3)
   endif

   if (use_p_diag) then

   allocate(p_diag_1d_f_pv1(llb1:uub1,llb2:uub2,llb3:uub3,llb4:uub4), &
            stat=ierr)
   if (idebug >= 1) then
     call assert(ierr.eq.0,'alloc(p_diag_1d_f_pv1) ',ierr)
   endif

   p_diag_1d_f_pv1 = diag_1d_f_pv1

   endif

   if (need_diag_1d_df_pv1) then
     llb1_df = lbound(diag_1d_df_pv1,1)
     llb2_df = lbound(diag_1d_df_pv1,2)
     llb3_df = lbound(diag_1d_df_pv1,3)
     llb4_df = lbound(diag_1d_df_pv1,4)

     uub1_df = ubound(diag_1d_df_pv1,1)
     uub2_df = ubound(diag_1d_df_pv1,2)
     uub3_df = ubound(diag_1d_df_pv1,3)
     uub4_df = ubound(diag_1d_df_pv1,4)

     if (use_p_diag) then

     allocate(p_diag_1d_df_pv1(llb1_df:uub1_df, &
                llb2_df:uub2_df,llb3_df:uub3_df,llb4_df:uub4_df),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_1d_df_pv1) ',ierr)
     endif
     p_diag_1d_df_pv1 = 0

     p_diag_1d_df_pv1 = diag_1d_df_pv1

     endif
   endif

  if(diag_heat_on_host) then 
   if (need_diag_heat_pv) then
     lb1_pv = lbound(diag_heat_pv,1)
     lb2_pv = lbound(diag_heat_pv,2)
     lb3_pv = lbound(diag_heat_pv,3)
     lb4_pv = lbound(diag_heat_pv,4)
     lb5_pv = lbound(diag_heat_pv,5)
     lb6_pv = lbound(diag_heat_pv,6)

     ub1_pv = ubound(diag_heat_pv,1)
     ub2_pv = ubound(diag_heat_pv,2)
     ub3_pv = ubound(diag_heat_pv,3)
     ub4_pv = ubound(diag_heat_pv,4)
     ub5_pv = ubound(diag_heat_pv,5)
     ub6_pv = ubound(diag_heat_pv,6)

     if (use_p_diag) then

     allocate(p_diag_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv, &
            lb4_pv:ub4_pv,lb5_pv:ub5_pv,lb6_pv:ub6_pv),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_heat_pv) ',ierr)
     endif
     p_diag_heat_pv = 0
     p_diag_heat_pv = diag_heat_pv
     endif
   endif

   if (need_diag_heat_pv_psi) then
     lb1_psi = lbound(diag_heat_pv_psi,1)
     lb2_psi = lbound(diag_heat_pv_psi,2)
     lb3_psi = lbound(diag_heat_pv_psi,3)
     lb4_psi = lbound(diag_heat_pv_psi,4)
     lb5_psi = lbound(diag_heat_pv_psi,5)

     ub1_psi = ubound(diag_heat_pv_psi,1)
     ub2_psi = ubound(diag_heat_pv_psi,2)
     ub3_psi = ubound(diag_heat_pv_psi,3)
     ub4_psi = ubound(diag_heat_pv_psi,4)
     ub5_psi = ubound(diag_heat_pv_psi,5)

     if (use_p_diag) then

     allocate(p_diag_heat_pv_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi, &
            lb3_psi:ub3_psi,lb4_psi:ub4_psi,lb5_psi:ub5_psi),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_heat_pv_psi) ',ierr)
     endif
     p_diag_heat_pv_psi = 0
     p_diag_heat_pv_psi = diag_heat_pv_psi

     endif
   endif
 endif ! diag_heat_on

  if(diag_eflux_on_host) then
   if (need_diag_1d_eflux_pv) then
     lb1_ef = lbound(diag_1d_eflux_pv,1)
     lb2_ef = lbound(diag_1d_eflux_pv,2)
     lb3_ef = lbound(diag_1d_eflux_pv,3)
     lb4_ef = lbound(diag_1d_eflux_pv,4)
     lb5_ef = lbound(diag_1d_eflux_pv,5)

     ub1_ef = ubound(diag_1d_eflux_pv,1)
     ub2_ef = ubound(diag_1d_eflux_pv,2)
     ub3_ef = ubound(diag_1d_eflux_pv,3)
     ub4_ef = ubound(diag_1d_eflux_pv,4)
     ub5_ef = ubound(diag_1d_eflux_pv,5)

     if (use_p_diag) then

     allocate(p_diag_1d_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef, &
            lb4_ef:ub4_ef,lb5_ef:ub5_ef),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_1d_eflux_pv) ',ierr)
     endif
     p_diag_1d_eflux_pv = 0
     p_diag_1d_eflux_pv = diag_1d_eflux_pv

     endif
   endif

   if (need_diag_2d_dflux_pv) then
     lb1_2df = lbound(diag_2d_dflux_pv,1)
     lb2_2df = lbound(diag_2d_dflux_pv,2)
     lb3_2df = lbound(diag_2d_dflux_pv,3)
     lb4_2df = lbound(diag_2d_dflux_pv,4)
     lb5_2df = lbound(diag_2d_dflux_pv,5)

     ub1_2df = ubound(diag_2d_dflux_pv,1)
     ub2_2df = ubound(diag_2d_dflux_pv,2)
     ub3_2df = ubound(diag_2d_dflux_pv,3)
     ub4_2df = ubound(diag_2d_dflux_pv,4)
     ub5_2df = ubound(diag_2d_dflux_pv,5)

     if (use_p_diag) then

     allocate(p_diag_2d_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df, &
            lb3_2df:ub3_2df,lb4_2df:ub4_2df,lb5_2df:ub5_2df),&
            stat=ierr)
     if (idebug >= 1) then
       call assert(ierr.eq.0,'alloc(p_diag_2d_dflux_pv) ',ierr)
     endif
     p_diag_2d_dflux_pv = 0
     p_diag_2d_dflux_pv = diag_2d_dflux_pv

     endif
   endif
 endif ! diag_eflux_on

!  ------------------
!  perform sum reduce
!  into diag_1d_f_pv1(:,:,:,1) so
!  that it is independent of sml_nthreads
!  ------------------
   if (use_dsum) then
     allocate( dsum_f(llb1:uub1,llb2:uub2,llb3:uub3), stat=ierr)
     call assert(ierr.eq.0,'alloc(dsum_f) ',ierr)

     dsum_f = 0

     if (use_p_diag) then
     do i4=llb4,uub4
!$omp parallel do private(i1,i2,i3)
     do i3=lb3,ub3
     do i2=lb2,ub2
     do i1=lb1,ub1
       dsum_f(i1,i2,i3) = dsum_f(i1,i2,i3) + &
          p_diag_1d_f_pv1(i1,i2,i3,i4)
     enddo
     enddo
     enddo
     enddo

     else
       n1  = size(diag_1d_f_pv1,1)
       n2  = size(diag_1d_f_pv1,2)
       n3  = size(diag_1d_f_pv1,3)
       n4  = size(diag_1d_f_pv1,4)

       l1 = lbound(diag_1d_f_pv1,1)
       l2 = lbound(diag_1d_f_pv1,2)
       l3 = lbound(diag_1d_f_pv1,3)
       l4 = lbound(diag_1d_f_pv1,4)

       u1 = ubound(diag_1d_f_pv1,1)
       u2 = ubound(diag_1d_f_pv1,2)
       u3 = ubound(diag_1d_f_pv1,3)
       u4 = ubound(diag_1d_f_pv1,4)
       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_f_pv1, dsum_f )
       call sum2_gpu( n1*n2*n3, n4, diag_1d_f_pv1, dsum_f )
     endif

     i4 = 1
!$omp parallel do private(i1,i2,i3)
     do i3=lb3,ub3
     do i2=lb2,ub2
     do i1=lb1,ub1
       diag_1d_f_pv1_host(i1,i2,i3,1) = &
         diag_1d_f_pv1_host(i1,i2,i3,1) + dsum_f(i1,i2,i3)
     enddo
     enddo
     enddo

     if (need_diag_1d_df_pv1) then
       allocate( dsum_df(llb1:uub1,llb2:uub2,llb3:uub3), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_df) ',ierr)
       endif

       dsum_df = 0


       if (use_p_diag) then

       do i4=llb4_df,uub4_df
!$omp parallel do private(i1,i2,i3)
       do i3=lb3_df,ub3_df
       do i2=lb2_df,ub2_df
       do i1=lb1_df,ub1_df
        dsum_df(i1,i2,i3) = dsum_df(i1,i2,i3) + &
          p_diag_1d_df_pv1(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo

       else


       n1  = size(diag_1d_df_pv1,1)
       n2  = size(diag_1d_df_pv1,2)
       n3  = size(diag_1d_df_pv1,3)
       n4  = size(diag_1d_df_pv1,4)

       l1 = lbound(diag_1d_df_pv1,1)
       l2 = lbound(diag_1d_df_pv1,2)
       l3 = lbound(diag_1d_df_pv1,3)
       l4 = lbound(diag_1d_df_pv1,4)

       u1 = ubound(diag_1d_df_pv1,1)
       u2 = ubound(diag_1d_df_pv1,2)
       u3 = ubound(diag_1d_df_pv1,3)
       u4 = ubound(diag_1d_df_pv1,4)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3, n4, diag_1d_df_pv1, dsum_df )

       endif


       i4 = 1
!$omp parallel do private(i1,i2,i3)
       do i3=lb3_df,ub3_df
       do i2=lb2_df,ub2_df
       do i1=lb1_df,ub1_df
         diag_1d_df_pv1_host(i1,i2,i3,1) = &
           diag_1d_df_pv1_host(i1,i2,i3,1) + dsum_df(i1,i2,i3)
       enddo
       enddo
       enddo
     endif

    if(diag_heat_on_host) then
     if (need_diag_heat_pv) then
       allocate( dsum_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv, &
                 lb4_pv:ub4_pv,lb5_pv:ub5_pv), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_heat_pv) ',ierr)
       endif

       dsum_heat_pv = 0


       if (use_p_diag) then

       do i6=lb6_pv,ub6_pv
!$omp parallel do private(i1,i2,i3,i4,i5)
       do i5=lb5_pv,ub5_pv
       do i4=lb4_pv,ub4_pv
       do i3=lb3_pv,ub3_pv
       do i2=lb2_pv,ub2_pv
       do i1=lb1_pv,ub1_pv
        dsum_heat_pv(i1,i2,i3,i4,i5) = dsum_heat_pv(i1,i2,i3,i4,i5) + &
          p_diag_heat_pv(i1,i2,i3,i4,i5,i6)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1  = size(diag_heat_pv,1)
       n2  = size(diag_heat_pv,2)
       n3  = size(diag_heat_pv,3)
       n4  = size(diag_heat_pv,4)
       n5  = size(diag_heat_pv,5)
       n6  = size(diag_heat_pv,6)

       l1 = lbound(diag_heat_pv,1)
       l2 = lbound(diag_heat_pv,2)
       l3 = lbound(diag_heat_pv,3)
       l4 = lbound(diag_heat_pv,4)
       l5 = lbound(diag_heat_pv,5)
       l6 = lbound(diag_heat_pv,6)

       u1 = ubound(diag_heat_pv,1)
       u2 = ubound(diag_heat_pv,2)
       u3 = ubound(diag_heat_pv,3)
       u4 = ubound(diag_heat_pv,4)
       u5 = ubound(diag_heat_pv,5)
       u6 = ubound(diag_heat_pv,6)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4*n5, n6, diag_heat_pv, dsum_heat_pv )

       endif

       i6 = 1
!$omp parallel do private(i1,i2,i3,i4,i5)
       do i5=lb5_pv,ub5_pv
       do i4=lb4_pv,ub4_pv
       do i3=lb3_pv,ub3_pv
       do i2=lb2_pv,ub2_pv
       do i1=lb1_pv,ub1_pv
         diag_heat_pv_host(i1,i2,i3,i4,i5,1) = &
           diag_heat_pv_host(i1,i2,i3,i4,i5,1) + dsum_heat_pv(i1,i2,i3,i4,i5)
       !if(diag_heat_pv_host(i1,i2,i3,i4,0,1)/=0) print *, "nonzero diag_heat_pv_host(i1,i2,i3,i4,0,1)"
       enddo
       enddo
       enddo
       enddo
       enddo
     endif


     if (need_diag_heat_pv_psi) then
       allocate( dsum_heat_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi, &
                 lb4_psi:ub4_psi), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_heat_psi) ',ierr)
       endif

       dsum_heat_psi = 0

       if (use_p_diag) then
       do i5=lb5_psi,ub5_psi
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_psi,ub4_psi
       do i3=lb3_psi,ub3_psi
       do i2=lb2_psi,ub2_psi
       do i1=lb1_psi,ub1_psi
        dsum_heat_psi(i1,i2,i3,i4) = dsum_heat_psi(i1,i2,i3,i4) + &
          p_diag_heat_pv_psi(i1,i2,i3,i4,i5)
          !if(p_diag_heat_pv_psi(i1,i2,i3,0,i5)/=0) print *, "nonzero p_diag_heat_pv_psi(i1,i2,i3,i4,i5)"
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1  = size(diag_heat_pv_psi,1)
       n2  = size(diag_heat_pv_psi,2)
       n3  = size(diag_heat_pv_psi,3)
       n4  = size(diag_heat_pv_psi,4)
       n5  = size(diag_heat_pv_psi,5)

       l1 = lbound(diag_heat_pv_psi,1)
       l2 = lbound(diag_heat_pv_psi,2)
       l3 = lbound(diag_heat_pv_psi,3)
       l4 = lbound(diag_heat_pv_psi,4)
       l5 = lbound(diag_heat_pv_psi,5)

       u1 = ubound(diag_heat_pv_psi,1)
       u2 = ubound(diag_heat_pv_psi,2)
       u3 = ubound(diag_heat_pv_psi,3)
       u4 = ubound(diag_heat_pv_psi,4)
       u5 = ubound(diag_heat_pv_psi,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_heat_pv_psi, dsum_heat_psi )

       endif

       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_psi,ub4_psi
       do i3=lb3_psi,ub3_psi
       do i2=lb2_psi,ub2_psi
       do i1=lb1_psi,ub1_psi
         diag_heat_pv_psi_host(i1,i2,i3,i4,1) = &
           diag_heat_pv_psi_host(i1,i2,i3,i4,1) + dsum_heat_psi(i1,i2,i3,i4)
           !if(abs(diag_heat_pv_psi_host(i1,i2,i3,0,1))>0) print *,"nonzero diag_heat_psi_cpu"
       enddo
       enddo
       enddo
       enddo
     endif
    endif ! diag_heat_on

    if(diag_eflux_on_host) then
     if (need_diag_1d_eflux_pv) then
       allocate( dsum_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef, &
                 lb4_ef:ub4_ef), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_eflux_pv) ',ierr)
       endif

       dsum_eflux_pv = 0

       if (use_p_diag) then

       do i5=lb5_ef,ub5_ef
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_ef,ub4_ef
       do i3=lb3_ef,ub3_ef
       do i2=lb2_ef,ub2_ef
       do i1=lb1_ef,ub1_ef
        dsum_eflux_pv(i1,i2,i3,i4) = dsum_eflux_pv(i1,i2,i3,i4) + &
          p_diag_1d_eflux_pv(i1,i2,i3,i4,i5)
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1  = size(diag_1d_eflux_pv,1)
       n2  = size(diag_1d_eflux_pv,2)
       n3  = size(diag_1d_eflux_pv,3)
       n4  = size(diag_1d_eflux_pv,4)
       n5  = size(diag_1d_eflux_pv,5)

       l1 = lbound(diag_1d_eflux_pv,1)
       l2 = lbound(diag_1d_eflux_pv,2)
       l3 = lbound(diag_1d_eflux_pv,3)
       l4 = lbound(diag_1d_eflux_pv,4)
       l5 = lbound(diag_1d_eflux_pv,5)

       u1 = ubound(diag_1d_eflux_pv,1)
       u2 = ubound(diag_1d_eflux_pv,2)
       u3 = ubound(diag_1d_eflux_pv,3)
       u4 = ubound(diag_1d_eflux_pv,4)
       u5 = ubound(diag_1d_eflux_pv,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_1d_eflux_pv, dsum_eflux_pv )

       endif
       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_ef,ub4_ef
       do i3=lb3_ef,ub3_ef
       do i2=lb2_ef,ub2_ef
       do i1=lb1_ef,ub1_ef
         diag_1d_eflux_pv_host(i1,i2,i3,i4,1) = &
           diag_1d_eflux_pv_host(i1,i2,i3,i4,1) + dsum_eflux_pv(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo
     endif

     if (need_diag_2d_dflux_pv) then
       allocate( dsum_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df, &
                 lb4_2df:ub4_2df), stat=ierr)
       if (idebug >= 1) then
         call assert(ierr.eq.0,'alloc(dsum_dflux_pv) ',ierr)
       endif

       dsum_dflux_pv = 0

       if (use_p_diag) then

       do i5=lb5_2df,ub5_2df
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_2df,ub4_2df
       do i3=lb3_2df,ub3_2df
       do i2=lb2_2df,ub2_2df
       do i1=lb1_2df,ub1_2df
        dsum_dflux_pv(i1,i2,i3,i4) = dsum_dflux_pv(i1,i2,i3,i4) + &
          p_diag_2d_dflux_pv(i1,i2,i3,i4,i5)
       enddo
       enddo
       enddo
       enddo
       enddo

       else

       n1  = size(diag_2d_dflux_pv,1)
       n2  = size(diag_2d_dflux_pv,2)
       n3  = size(diag_2d_dflux_pv,3)
       n4  = size(diag_2d_dflux_pv,4)
       n5  = size(diag_2d_dflux_pv,5)

       l1 = lbound(diag_2d_dflux_pv,1)
       l2 = lbound(diag_2d_dflux_pv,2)
       l3 = lbound(diag_2d_dflux_pv,3)
       l4 = lbound(diag_2d_dflux_pv,4)
       l5 = lbound(diag_2d_dflux_pv,5)

       u1 = ubound(diag_2d_dflux_pv,1)
       u2 = ubound(diag_2d_dflux_pv,2)
       u3 = ubound(diag_2d_dflux_pv,3)
       u4 = ubound(diag_2d_dflux_pv,4)
       u5 = ubound(diag_2d_dflux_pv,5)

       !call sum4_gpu( l1,u1,l2,u2,l3,u3,l4,u4, diag_1d_df_pv1, dsum_df )
       call sum2_gpu( n1*n2*n3*n4, n5, diag_2d_dflux_pv, dsum_dflux_pv )

       endif
       i5 = 1
!$omp parallel do private(i1,i2,i3,i4)
       do i4=lb4_2df,ub4_2df
       do i3=lb3_2df,ub3_2df
       do i2=lb2_2df,ub2_2df
       do i1=lb1_2df,ub1_2df
         diag_2d_dflux_pv_host(i1,i2,i3,i4,1) = &
           diag_2d_dflux_pv_host(i1,i2,i3,i4,1) + dsum_dflux_pv(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo
     endif
    endif ! diag_eflux_on

     if (allocated(dsum_f)) then
       deallocate( dsum_f, stat=ierr )
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_f) ',ierr)
       endif
     endif

     if (allocated(dsum_df)) then
       deallocate( dsum_df, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_df) ',ierr)
       endif
     endif

   if(diag_heat_on_host) then
     if (allocated(dsum_heat_pv)) then
       deallocate( dsum_heat_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_heat_pv) ',ierr)
       endif
     endif

     if (allocated(dsum_heat_psi)) then
       deallocate( dsum_heat_psi, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_heat_psi) ',ierr)
       endif
     endif
   endif ! diag_heat_on

   if(diag_eflux_on_host) then
     if (allocated(dsum_eflux_pv)) then
       deallocate( dsum_eflux_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_eflux_pv) ',ierr)
       endif
     endif

     if (allocated(dsum_dflux_pv)) then
       deallocate( dsum_dflux_pv, stat=ierr)
       if (idebug >= 2) then
         call assert( ierr.eq.0,'dealloc(dsum_dflux_psi) ',ierr)
       endif
     endif
   endif ! diag_eflux_on

   else
   do i4=llb4,uub4
     diag_1d_f_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) =  &
       diag_1d_f_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) + &
         p_diag_1d_f_pv1(lb1:ub1,lb2:ub2,lb3:ub3, i4 )
   enddo

    if (need_diag_1d_df_pv1) then
     do i4=llb4_df,uub4_df
      diag_1d_df_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) =  &
       diag_1d_df_pv1_host(lb1:ub1,lb2:ub2,lb3:ub3,1) + &
         p_diag_1d_df_pv1(lb1:ub1,lb2:ub2,lb3:ub3, i4 )
     enddo
    endif

   if(diag_heat_on_host) then
    if (need_diag_heat_pv) then
     do i6=lb6_pv,ub6_pv
      diag_heat_pv_host(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,1) =  &
       diag_heat_pv_host(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,1) + &
         p_diag_heat_pv(lb1_pv:ub1_pv,lb2_pv:ub2_pv,lb3_pv:ub3_pv,lb4_pv:ub4_pv,lb5_pv:ub5_pv,i6 )
     enddo
    endif

    if (need_diag_heat_pv_psi) then
     do i5=lb5_psi,ub5_psi
      diag_heat_pv_psi_host(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,1) =  &
       diag_heat_pv_psi_host(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,1) + &
       p_diag_heat_pv_psi(lb1_psi:ub1_psi,lb2_psi:ub2_psi,lb3_psi:ub3_psi,lb4_psi:ub4_psi,i5 )
     enddo
    endif
   endif ! diag_heat_on

   if(diag_eflux_on_host) then
    if (need_diag_1d_eflux_pv) then
     do i5=lb5_ef,ub5_ef
      diag_1d_eflux_pv_host(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,1) =  &
       diag_1d_eflux_pv_host(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,1) + &
         p_diag_1d_eflux_pv(lb1_ef:ub1_ef,lb2_ef:ub2_ef,lb3_ef:ub3_ef,lb4_ef:ub4_ef,i5 )
     enddo
    endif

    if (need_diag_2d_dflux_pv) then
     do i5=lb5_2df,ub5_2df
      diag_2d_dflux_pv_host(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,1) =  &
       diag_2d_dflux_pv_host(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,1) + &
       p_diag_2d_dflux_pv(lb1_2df:ub1_2df,lb2_2df:ub2_2df,lb3_2df:ub3_2df,lb4_2df:ub4_2df,i5 )
     enddo
    endif

   endif ! diag_eflux_on

   endif

! -------
! clean up
! -------

   if (allocated(p_diag_1d_f_pv1)) then
     deallocate( p_diag_1d_f_pv1, stat=ierr )
     call assert( ierr.eq.0,'dealloc(p_diag_1d_f_pv1 ',ierr)
   endif


   if (allocated(p_diag_1d_df_pv1)) then
     deallocate( p_diag_1d_df_pv1, stat=ierr )
     call assert( ierr.eq.0,'dealloc(p_diag_1d_df_pv1 ',ierr)
   endif

   if (allocated(diag_1d_f_pv1)) then
     deallocate(diag_1d_f_pv1,stat=ierr)
     call assert(ierr.eq.0,'final dealloc(diag_1d_f_pv1),ierr=',ierr)
   endif

   if (allocated(diag_1d_df_pv1)) then
     deallocate( diag_1d_df_pv1,stat=ierr)
     call assert(ierr.eq.0,'final dealloc(diag_1d_df_pv1),ierr=',ierr)
   endif

   if(diag_eflux_on_host)then
     if (allocated(diag_1d_eflux_pv)) then
       deallocate(diag_1d_eflux_pv,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(diag_1d_eflux_pv),ierr=',ierr)
     endif
     if (allocated(diag_2d_dflux_pv)) then
       deallocate(diag_2d_dflux_pv,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(diag_2d_dflux_pv),ierr=',ierr)
     endif
   endif ! diag_eflux_on

   if(diag_heat_on_host)then
     if (allocated(diag_heat_pv)) then
        deallocate(diag_heat_pv,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(diag_heat_pv),ierr=',ierr)
     endif
     if (allocated(diag_heat_pv_psi)) then
        deallocate(diag_heat_pv_psi,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(diag_heat_pv_psi),ierr=',ierr)
     endif
     if (allocated(p_diag_heat_pv)) then
        deallocate(p_diag_heat_pv,stat=ierr)
        call assert(ierr.eq.0,'final dealloc(p_diag_heat_pv),ierr=',ierr)
     endif
     if (allocated(p_diag_heat_pv_psi)) then
        deallocate(p_diag_heat_pv_psi,stat=ierr)
       call assert(ierr.eq.0,'final dealloc(p_diag_heat_pv_psi),ierr=',ierr)
     endif

   endif ! diag_eflux_on


   return
   end subroutine update_host_diag
  
end module diag_module_gpu
