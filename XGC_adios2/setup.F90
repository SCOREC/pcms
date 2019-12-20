subroutine setup(grid,psn,spall)
  use sml_module
  use ptl_module
  use psn_class
  use itp_module
  use eq_module
  use col_module
  use diag_module
  use bnc_module
  use grid_class
  use smooth_module
  use pol_decomp_module
  use src_module
  use rad_module
  use perf_monitor
  use random_xgc
  use neu_module
  use xgc_interfaces
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h>
#else
#include <petsc/finclude/petsc.h>
#endif
  type(grid_type), target :: grid
  type(psn_type), target :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  !
  PetscErrorCode::ierr
#ifdef _OPENMP
  integer omp_get_max_threads
  external omp_get_max_threads
#endif
  
!pw  call t_startf("SETUP_INIT")
#ifdef _OPENMP
!$OMP PARALLEL
!$OMP MASTER
  sml_nthreads = omp_get_max_threads()
!$OMP END MASTER
!$OMP END PARALLEL
#else
  sml_nthreads = 1
#endif
  if(sml_mype==0) print *, '# of OMP threads =',sml_nthreads
  if(sml_mype==0) print *, 'Total simulation processors : ',sml_totalpe


  if(sml_mype==0) print *, 'setting up.........',sml_mype
!pw initialize random number seeds moved to beginning of post_process_input
!pw call init_ranx()

  ! setup constant
  sml_pi=4D0*datan(1D0)
  sml_2pi=8D0*datan(1D0)
  sml_sqrtpi=sqrt(sml_pi)
  sml_sqrt2pi=sqrt(sml_2pi)
  !sml_2pi_inv=1D0/sml_2pi
!pw  call t_stopf("SETUP_INIT")

  ! Setup parameters and poisson array
  if(sml_mype==0) print *, 'Read input file'
  call t_startf("READ_INPUT")
  call read_input
  call t_stopf("READ_INPUT")
  call t_startf("POST_PROCESS_INPUT")
  call post_process_input(spall)
  call t_stopf("POST_PROCESS_INPUT")
  if(sml_mype==0) print *, 'setting up grid system.....',sml_mype
  call t_startf("SETUP_GRID")
  call setup_grid(grid,psn)
  call t_stopf("SETUP_GRID")

  call check_point('diagnosis_1d init')

!pw  !This routine should be done before volume calculation..
!pw  call t_startf("DIAG_1D_INIT")
!pw  (called twice) call diag_1d_init
!pw  call t_stopf("DIAG_1D_INIT")

  ! Get space volume using Monte-Carlo method
  if(sml_mype==0) print *, 'Monte-Carlo volume calculation...'
  call t_startf("SETUP_GET_VOLUME")
  call t_adj_detailf(+1)
  call get_volume(grid,psn,spall)
  call t_adj_detailf(-1)
  call t_stopf("SETUP_GET_VOLUME")

  call check_point('smoothing initializing')
!pw  call t_startf("SMOOTH_R_INITS")
  call smooth_r_init1(smooth_r1,smooth_r1_d0_in,smooth_r1_n_in,smooth_r1_type_in,grid%nnode)
  call smooth_r_init2(smooth_r1,grid)
!pw  call t_stopf("SMOOTH_R_INITS")
!  call smooth_nearx_init(smooth_nearx1,grid,smooth_nearx_d0,smooth_nearx_dr,smooth_nearx_dz,smooth_nearx_nsmth)
!  call check_point('init_f0')
!  call init_f0(grid)
!pw  call t_stopf("INIT_F0")

  ! particle loading
  if(sml_mype==0) print *, 'load',sml_mype
  call t_startf("LOAD")
  ! Sheath pot is in restart file, need to allocate memory before call to load
  if (sml_sheath_mode /= 0) call init_sheath(grid,psn)

  call extend_boundary2(grid,psn)  ! moved after init_sheath
  

  call load(grid,psn,spall)
  call t_stopf("LOAD")
  if(sml_mype==0)  print *, 'loading end', sml_mype
  if(sml_sheath_mode/=0 .and. .not. sml_restart) then
     if(sml_mype==0)  print *, 'sheath_pot_init',sml_mype
     call t_startf("SHEATH_POT_INIT")
     call sheath_pot_init(grid,psn)
     call t_stopf("SHEATH_POT_INIT")
  endif

#if defined(ADIOS)
  ! Dump grid description data in binary file for post-processing and viz
  ! sheath diag information is included. it need to be called after sheath_pot_init when sml_sheat_mode/=0
  if (sml_mype==0) then
    call t_startf("DUMP_GRID")
    call dump_grid(grid,psn)
    call t_stopf("DUMP_GRID")
  endif
#endif

  ! setup convert 2d grid to 1d00
  call t_startf("SETUP_CONV_GRID_INIT")
  call convert_grid_init(grid,10)   ! nsub should be input parameter
  call t_stopf("SETUP_CONV_GRID_INIT")

  ! init poisson solver moved from setup_grid -- need to be called after convert_grid_init
  call check_point('Init poisson')
  ! PETSc bug on OSX
  if (PETSC_COMM_SELF==0) call MPI_Comm_Dup(MPI_COMM_SELF,PETSC_COMM_SELF,PETSC_NULL_INTEGER)
#ifdef XGC1_EM
  if (sml_use_ts_solver) then
     call ts_create(psn%ts) ! just initilizes to zero, etc.
     call ts_init(grid,psn,psn%pbd0_2,psn%ts,ierr);CHKERRQ(ierr)
     !rh We need the grid information in ts and bfollow information
     !rh this needs to go somewhere else, has to be done only once
     psn%ts%grid => grid
     psn%ts%bd => psn%pbd0_2
#ifdef OLD_INIT_FF
     psn%ts%tr => psn%bfollow_tr
     psn%ts%p => psn%bfollow_p
     psn%ts%dx => psn%bfollow_1_dx
#else
     psn%ts%tr => psn%ff_1dp_tr
     psn%ts%p => psn%ff_1dp_p
     psn%ts%dx => psn%ff_1dp_dx
#endif
  else
     ! explicit solver init
     call init_poisson(grid,psn,psn%pbd0_2)
     ! create mass KSP solver for J
     call KSPCreate(psn%solver00%comm,psn%kspmass,ierr);CHKERRQ(ierr)
     call KSPSetOperators(psn%kspmass,psn%solver00%rhs_mat,psn%solver00%rhs_mat,ierr);CHKERRQ(ierr)
     call KSPSetOptionsPrefix(psn%kspmass,'mass_',ierr)
     call KSPSetFromOptions(psn%kspmass,ierr);CHKERRQ(ierr)
  end if
#else
  call t_startf("SETUP_POISSON")
  call t_adj_detailf(+1)
  call init_poisson(grid,psn,3)
  !initialize local variable for poisson subroutine
  call solve_poisson(grid,psn,-1)  ! init local vars in solve routine
  call t_adj_detailf(-1)
  call t_stopf("SETUP_POISSON")
#endif

  if(sml_pol_decomp) then
     call check_point('setup poloidal decomposition')
     call t_startf("SETUP_POL_DECOMP")
     call setup_pol_decomp(grid,spall)
     call t_stopf("SETUP_POL_DECOMP")
  endif

  if(sml_mype==0) then
     print *, 'initial diagnosis'
     call t_startf("OUTPUT_BFIELD")
     call output_bfield
     call output_volumes(grid)
     call t_stopf("OUTPUT_BFIELD")
#if defined(ADIOS)
     ! Dump bfield component data in binary file for post-processing and viz
     print *, 'dump_bfield'
     call t_startf("DUMP_BFIELD")
     call dump_bfield(grid)
     call t_stopf("DUMP_BFIELD")
#endif
     !     call output_psi_der
  endif

  call check_point('background elec charge')
  if(.not. sml_electron_on) then
    call t_startf("CHARGEE_BACKGROUND")
    call t_adj_detailf(+1)
    call chargee_background(grid,psn,spall)
    call t_adj_detailf(-1)
    call t_stopf("CHARGEE_BACKGROUND")
  endif
  
  !initialize col_3 variables
  if(col_mode==3) then
#ifdef VPIC_COL
     call col_3_setup(grid)
#endif
  endif !for col_3

  ! ES : Below col_mode==4 doesn't have to be here.
  ! Plz move to proper place if you want...
  if(col_mode==4) then
     call check_point('col_f_setup')
     call col_f_setup
  endif

  if(sml_mype==0) then
!pw    call t_startf("WRITE_RUNINFO")
    call write_runinfo
!pw    call t_stopf("WRITE_RUNINFO")
  endif
  call check_point('End of setup')
end subroutine setup

subroutine post_process_input(spall)
  use sml_module
  use ptl_module
  use itp_module
  use eq_module
  use col_module
  use diag_module
  use bnc_module
  use grid_class
  use smooth_module
  use pol_decomp_module
  use src_module
  use rad_module
  use neu_module
  use random_xgc
  use perf_monitor
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i

  sml_plane_per_pe=1
  sml_pe_per_plane=sml_totalpe/sml_nphi_total

  if(sml_totalpe*sml_plane_per_pe/=sml_nphi_total*sml_pe_per_plane) then
     print *, '*******************************************************'
     print *, 'Invalid PE number and sml_nphi_total'
     print *, '*******************************************************'
     stop
  endif
  if(sml_totalpe < sml_nphi_total* sml_grid_nrho*2 ) then
     print *, '*******************************************************'
     print *, 'Not enough mpi PE: # of PE should >= ', sml_nphi_total*sml_grid_nrho*2
     print *, 'sml_nphi_total=',sml_nphi_total, 'sml_grid_nrho*2=',sml_grid_nrho*2

     print *, '*******************************************************'
     stop
  endif

    ! Make a list of te processes in the new communicator
  call new_communicator
  if(col_mode == 3) then
      call new_mpi_ptl_type
  endif

  ! initialize seeds
  call init_ranx()

  ! set b-field sign
  if(sml_invert_b) then
     sml_bp_sign=-1D0
  else
     sml_bp_sign=1D0
  endif
  !
  if(sml_co_cur_bt) then
     sml_bt_sign=-1D0*sml_bp_sign
  else
     sml_bt_sign=1D0*sml_bp_sign
  endif

  if (.not.sml_use_ts_solver) then
  if(sml_poisson_solver_type>0 .and. sml_iter_solver) then
     sml_iter_solver=.false.
     if(sml_mype==0) print *, '*** sml_poisson_solver_type is non-zero. sml_iter_solver is forced to .false.', sml_poisson_solver_type
  endif
  if(sml_poisson_solver_type>0 .and. sml_use_simple00) then
     sml_use_simple00=.false.
     if(sml_mype==0) print *, '*** sml_poisson_solver_type is non-zero. sml_use_simple00 is forced to .false.', sml_poisson_solver_type
  endif
  if(sml_fsa_solver  .and. sml_poisson_solver_type >0) then
     sml_fsa_solver = .false.
     if(sml_mype==0) print *, ' **** ERROR: sml_fsa_solver is obsolete input. sml_poisson_solver_type is non-zero, so sml_fsa_sover = .true. is ignored'
  endif
  if(sml_fsa_solver .and. sml_poisson_solver_type ==-1) then
     sml_poisson_solver_type = 1
     if(sml_mype==0) print *, 'sml_poisson_solver_type is set to 1,due to sml_fsa_sover = .true.'
  endif
  ! no input for sml_poisson_solver_type
  if(sml_poisson_solver_type == -1 ) then
     sml_poisson_solver_type = 0
  endif
  endif

  sml_nphi_2pi=sml_nphi_total*sml_wedge_n
  sml_2pi_wedge_n=sml_2pi/sml_wedge_n

  ! scale ptl_num by number of threads
  ptl_num = ptl_num * sml_nthreads
  ptl_maxnum = ptl_maxnum * sml_nthreads


  ptl_e_num = ptl_e_num * sml_nthreads
  ptl_e_maxnum = ptl_e_maxnum * sml_nthreads

  spall(1)%num=ptl_num
  spall(1)%min_max_num=ptl_num
  spall(0)%num=ptl_e_num
  spall(0)%min_max_num=ptl_e_num

  ptl_mass(1)=ptl_mass_au*sml_prot_mass
  ptl_charge(1)=ptl_charge_eu*sml_e_charge
  ptl_mass(0)=ptl_e_mass_au*sml_prot_mass
  ptl_charge(0)=ptl_e_charge_eu*sml_e_charge

  ptl_c_m(:)=ptl_charge(:)/ptl_mass(:)
  ptl_c2_2m(:)=0.5D0*ptl_charge(:)**2/ptl_mass(:)

  ptl_deltaf_sp(0)=sml_deltaf_elec
  ptl_deltaf_sp(1)=sml_deltaf

  if(sml_electron_on) then
     ptl_isp=0
  else
     ptl_isp=1
  endif
  ptl_nsp=1

  if(sml_mype==0) print *, 'reading......'
  call add_dir_name1
  call xgc_read



  ! Second Definition part (before memory allocation and interpolation)---


  if(sml_mype==0) print *,'ptl_num for each CPU =', ptl_num

  ! maximum number of particle per processor
  if(ptl_maxnum < ptl_num ) then
     print *, 'too small memory space. Increase ptl_maxnum in module.f90 or decrease ptl_num. '
     stop
  endif

  !append dir to filename.
  call add_dir_name2

  if(sml_f0_grid .and. .not. sml_pol_decomp) then
     if(sml_mype==0) print *, '**** Warning: sml_pol_decomp = .True. enforced due to sml_f0_grid ****'
     sml_pol_decomp=.true.
  endif
  if(sml_f0_grid .and. sml_pol_decomp_simple) then
     if(sml_mype==0) print *, '**** Warning: sml_pol_decomp_simple = .False. enforced due to sml_f0_grid ****'
     sml_pol_decomp_simple=.false.
  endif

  ! memory allocation
  call check_point('memory allocation')
  call ptl_mem_allocation( spall(1), 1, ptl_maxnum  ,ptl_mass(1)  , ptl_charge(1),  sml_nlarmor, ptl_lost_nummax)
  call ptl_mem_allocation( spall(0), 0, ptl_e_maxnum,ptl_mass(0)  , ptl_charge(0),  1          , ptl_lost_nummax)
  call mem_allocation

  call check_point('init_interpolation')
  !if(sml_mype==0) print *, 'init_interpolation'
  call init_interpolation
  call check_point( 'third difinition part')

  ! energy order in MKS unit
  sml_en_order = sml_en_order_kev*1D3 * sml_e_charge
  ! transit time for a particle (sp=1) at the mag axis with pitch1(??) , R=1
  sml_tran = sml_2pi *eq_axis_r/ sqrt(2.D0 * sml_en_order/ptl_mass(1))
  ! delta t of one time step
  sml_dt= sml_dt*sml_tran ! Delta t - step size

  !simulation boundary
  sml_inpsi=sml_inpsi* eq_x_psi  !inner boundary of simulation region
  sml_outpsi=sml_outpsi * eq_x_psi ! outter boundary of simulation region

  sml_inv_f_source_period=1D0/real(sml_f_source_period,8)

!  sml_exb_suppress=0
!  sml_exb_suppress_time=0.5 * sml_tran
!  sml_exb_suppress_width=0.25 * sml_tran

!  sml_f0den_ped_c= sml_f0den_ped_c * eq_x_psi
!  sml_f0den_ped_width=sml_f0den_ped_width * eq_x_psi

  if(sml_deltaf_f0_mode==-1) then
     sml_f0_1_Ln   = sml_f0_1_Ln   /eq_axis_r
     sml_f0_1_Lt   = sml_f0_1_Lt   /eq_axis_r
     sml_f0_1_Lt_e = sml_f0_1_Lt_e /eq_axis_r
  elseif(sml_deltaf_f0_mode==-2) then ! psi
     sml_f0_1_Ln   = sml_f0_1_Ln   /sqrt(eq_x_psi)
     sml_f0_1_Lt   = sml_f0_1_Lt   /sqrt(eq_x_psi)
     sml_f0_1_Lt_e = sml_f0_1_Lt_e /sqrt(eq_x_psi)
  endif

  sml_bd_Te_width=sml_bd_Te_width*eq_x_psi ! radial decay length of electron temperature at the boundary


  !region checking parameters
  eq_x_slope =-(eq_x_r - eq_axis_r)/(eq_x_z - eq_axis_z)
  if(eq_set_x2) then
     call check_point('Put eq_x2_r and eq_x2_z in input file')
  else
     eq_x2_r=eq_axis_r
     eq_x2_z=eq_max_z
  endif
  eq_x2_slope =-(eq_x2_r - eq_axis_r)/(eq_x2_z - eq_axis_z)
  if(sml_mype==0) print *, 'eq_x_slope=', eq_x_slope
  if(sml_mype==0) print *, 'eq_x2_slope=',eq_x2_slope


  !density and temperature : tanh function --------------------------
  ! den - ion(electron) density
  eq_den_x1= eq_den_x1 * eq_x_psi
  eq_den_x2= eq_den_x2 * eq_x_psi
  eq_den_x3= eq_den_x3 * eq_x_psi

  ! tempi - ion temperature
  eq_tempi_x1=eq_tempi_x1 * eq_x_psi
  eq_tempi_x2=eq_tempi_x2 * eq_x_psi
  eq_tempi_x3=eq_tempi_x3 * eq_x_psi

  !electron temperature
  eq_tempe_x1=eq_tempe_x1 * eq_x_psi
  eq_tempe_x2=eq_tempe_x2 * eq_x_psi
  eq_tempe_x3=eq_tempe_x3 * eq_x_psi


  ! ion flow , electron flow
  eq_flowi_x1=eq_flowi_x1 * eq_x_psi
  eq_flowi_x2=eq_flowi_x2 * eq_x_psi
  eq_flowi_x3=eq_flowi_x3 * eq_x_psi

  ! electron flow
  eq_flowe_x1=eq_flowe_x1 * eq_x_psi
  eq_flowe_x2=eq_flowe_x2 * eq_x_psi
  eq_flowe_x3=eq_flowe_x3 * eq_x_psi


  ! z_eff
  eq_zeff_x1=eq_zeff_x1 * eq_x_psi
  eq_zeff_x2=eq_zeff_x2 * eq_x_psi
  eq_zeff_x3=eq_zeff_x3 * eq_x_psi


  ! initial density
  eq_den%shape=eq_den_shape
  eq_den%iny=(/eq_den_v1, eq_den_v2, eq_den_v3/)
  eq_den%inx=(/eq_den_x1, eq_den_x2, eq_den_x3/)
  eq_den%filename=eq_den_file

  call eq_ftn_setup(eq_den)
  ! initial ion temperature
  eq_tempi%shape=eq_tempi_shape
  eq_tempi%iny=(/eq_tempi_v1, eq_tempi_v2, eq_tempi_v3/)
  eq_tempi%inx=(/eq_tempi_x1, eq_tempi_x2, eq_tempi_x3/)
  eq_tempi%filename=eq_tempi_file
  call eq_ftn_setup(eq_tempi)

  ! initial ion flow
  eq_flowi%shape=eq_flowi_shape
  eq_flowi%iny=(/eq_flowi_v1, eq_flowi_v2, eq_flowi_v3/)
  eq_flowi%inx=(/eq_flowi_x1, eq_flowi_x2, eq_flowi_x3/)
  eq_flowi%filename=eq_flowi_file
  call eq_ftn_setup(eq_flowi)


  ! initial electron temperature
  eq_tempe%shape=eq_tempe_shape
  eq_tempe%iny=(/eq_tempe_v1, eq_tempe_v2, eq_tempe_v3/)
  eq_tempe%inx=(/eq_tempe_x1, eq_tempe_x2, eq_tempe_x3/)
  eq_tempe%filename=eq_tempe_file
  call eq_ftn_setup(eq_tempe)

  ! initial electron flow
  eq_flowe%shape=eq_flowe_shape
  eq_flowe%iny=(/eq_flowe_v1, eq_flowe_v2, eq_flowe_v3/)
  eq_flowe%inx=(/eq_flowe_x1, eq_flowe_x2, eq_flowe_x3/)
  eq_flowe%filename=eq_flowe_file
  call eq_ftn_setup(eq_flowe)


  ! initial electron flow
  eq_zeff%shape=eq_zeff_shape
  eq_zeff%iny=(/eq_zeff_v1, eq_zeff_v2, eq_zeff_v3/)
  eq_zeff%inx=(/eq_zeff_x1, eq_zeff_x2, eq_zeff_x3/)
  eq_zeff%filename=eq_zeff_file
  call eq_ftn_setup(eq_zeff)
  



  !-----------------------------------------------------------------------

  !collision input------------------------------------------------------------
  col_pin=col_pin*eq_x_psi
  col_pout=col_pout*eq_x_psi
  col_2_pin=col_2_pin*eq_x_psi
  col_2_pout=col_2_pout*eq_x_psi
  col_2_dp=(col_2_pout-col_2_pin)/real(col_2_m)
  col_2_dtheta= sml_2pi/real(col_2_mtheta)
  col_2_inv_dp=1D0/col_2_dp
  col_2_inv_dtheta=1D0/col_2_dtheta

  col_vb_pin=col_vb_pin*eq_x_psi
  col_vb_pout=col_vb_pout*eq_x_psi
  col_vb_dp=(col_vb_pout-col_vb_pin)/real(col_vb_m)
  col_vb_dtheta= sml_2pi/real(col_vb_mtheta)
  col_vb_inv_dp=1D0/col_vb_dp
  col_vb_inv_dtheta=1D0/col_vb_dtheta

  col_den_imp_ped_c = col_den_imp_ped_c * eq_x_psi
  col_den_imp_ped_width=col_den_imp_ped_width * eq_x_psi

  col_accel_pin1=col_accel_pin1*eq_x_psi
  col_accel_pout1=col_accel_pout1*eq_x_psi
  col_accel_pin2=col_accel_pin2*eq_x_psi
  col_accel_pout2=col_accel_pout2*eq_x_psi

  if (col_mode/=0) call init_col_module(ptl_isp, ptl_nsp, spall(1))


  ! diagnosis ---------------------------------------------------------------------
  diag_1d_pin=diag_1d_pin * eq_x_psi
  diag_1d_pout=diag_1d_pout * eq_x_psi
  !
  diag_neu_pin=diag_neu_pin * eq_x_psi
  diag_neu_pout=diag_neu_pout * eq_x_psi

  diag_f0_period=max(1,diag_f0_period/sml_f_source_period)*sml_f_source_period
  !
  call diag_1d_init
  call diag_neu_init
  if(diag_heat_on) then
     call check_point('diag_heat_init')
     call diag_heat_init
  endif

  ! smooth --------------------------------------------------------------------



  ! source -------------------

  src_pin1  = src_pin1  * eq_x_psi
  src_pout1 = src_pout1 * eq_x_psi
  src_decay_width1 = src_decay_width1 * eq_x_psi
  src_pin2  = src_pin2  * eq_x_psi
  src_pout2 = src_pout2 * eq_x_psi
  src_decay_width2 = src_decay_width2 * eq_x_psi
  src_pin3  = src_pin3  * eq_x_psi
  src_pout3 = src_pout3 * eq_x_psi
  src_decay_width3 = src_decay_width3 * eq_x_psi
  src_pin4  = src_pin4  * eq_x_psi
  src_pout4 = src_pout4 * eq_x_psi
  src_decay_width4 = src_decay_width4 * eq_x_psi

  src_pin1_e  = src_pin1_e  * eq_x_psi
  src_pout1_e = src_pout1_e * eq_x_psi
  src_decay_width1_e = src_decay_width1_e * eq_x_psi
  src_pin2_e  = src_pin2_e  * eq_x_psi
  src_pout2_e = src_pout2_e * eq_x_psi
  src_decay_width2_e = src_decay_width2_e * eq_x_psi
  src_pin3_e  = src_pin3_e  * eq_x_psi
  src_pout3_e = src_pout3_e * eq_x_psi
  src_decay_width3_e = src_decay_width3_e * eq_x_psi
  src_pin4_e  = src_pin4_e  * eq_x_psi
  src_pout4_e = src_pout4_e * eq_x_psi
  src_decay_width4_e = src_decay_width4_e * eq_x_psi

  call src_setup

  ! radiation -----------------------------------

  if(sml_radiation .and. .not. sml_neutral) then
     if(sml_mype==0) then
        print *, '****************************************** Warning **********************************************'
        print *, '*** sml_neutral is .false.'
        print *, '*** sml_radiation requires sml_neutral'
        print *, '*** sml_neutral is turned on and. neu_enforce_no_neutral is .true.'
     endif
     neu_enforce_no_neutral=.true.
     sml_neutral=.true.
  endif

  if(sml_radiation) then
     call check_point('init_radiation')
     call init_radiation
  endif


  ! ---------------------- others --------------------------------------------
     sml_bd_min_r=max(sml_bd_min_r,itp_min_r)
     sml_bd_max_r=min(sml_bd_max_r,itp_max_r)
     sml_bd_min_z=max(sml_bd_min_z,itp_min_z)
     sml_bd_max_z=min(sml_bd_max_z,itp_max_z)

     if(sml_mype==0) print *, 'R range :', sml_bd_min_r,sml_bd_max_r,&
       'Z range :',sml_bd_min_z,sml_bd_max_z
        bnc_min_r=sml_bd_min_r
        bnc_max_r=sml_bd_max_r
        bnc_dr=(bnc_max_r-bnc_min_r)/real(bnc_nr-1)

  ! -------- shift_ie communication and parallel algorithm options -------------
  !! for electrons
  spall(0)%shift_opt(index_shift_ie_max_nthreads)       = sml_elec_shift_max_nthreads
  spall(0)%shift_opt(index_shift_ie_large_limit)        = sml_elec_large_limit

  if (sml_elec_use_alltoall) then
     spall(0)%shift_opt(index_shift_ie_use_alltoall)    = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_alltoall)    = false_shift_ie_opt
  endif

  if (sml_elec_use_hs_barrier0) then
     spall(0)%shift_opt(index_shift_ie_use_hs_barrier0) = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_hs_barrier0) = false_shift_ie_opt
  endif

  if (sml_elec_handshake) then
     spall(0)%shift_opt(index_shift_ie_handshake)       = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_handshake)       = false_shift_ie_opt
  endif

  if (sml_elec_use_hs_barrier1) then
     spall(0)%shift_opt(index_shift_ie_use_hs_barrier1) = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_hs_barrier1) = false_shift_ie_opt
  endif

  if (sml_elec_use_sendrecv) then
     spall(0)%shift_opt(index_shift_ie_use_sendrecv)    = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_sendrecv)    = false_shift_ie_opt
  endif

  if (sml_elec_all_sendrecv) then
     spall(0)%shift_opt(index_shift_ie_all_sendrecv)    = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_all_sendrecv)    = false_shift_ie_opt
  endif

  if (sml_elec_use_isend) then
     spall(0)%shift_opt(index_shift_ie_use_isend)       = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_isend)       = false_shift_ie_opt
  endif

  if (sml_elec_use_rsend) then
     spall(0)%shift_opt(index_shift_ie_use_rsend)       = true_shift_ie_opt
  else
     spall(0)%shift_opt(index_shift_ie_use_rsend)       = false_shift_ie_opt
  endif

  !! for ions
  spall(1)%shift_opt(index_shift_ie_max_nthreads)       = sml_ion_shift_max_nthreads
  spall(1)%shift_opt(index_shift_ie_large_limit)        = sml_ion_large_limit

  if (sml_ion_use_alltoall) then
     spall(1)%shift_opt(index_shift_ie_use_alltoall)    = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_alltoall)    = false_shift_ie_opt
  endif

  if (sml_ion_use_hs_barrier0) then
     spall(1)%shift_opt(index_shift_ie_use_hs_barrier0) = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_hs_barrier0) = false_shift_ie_opt
  endif

  if (sml_ion_handshake) then
     spall(1)%shift_opt(index_shift_ie_handshake)       = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_handshake)       = false_shift_ie_opt
  endif

  if (sml_ion_use_hs_barrier1) then
     spall(1)%shift_opt(index_shift_ie_use_hs_barrier1) = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_hs_barrier1) = false_shift_ie_opt
  endif

  if (sml_ion_use_sendrecv) then
     spall(1)%shift_opt(index_shift_ie_use_sendrecv)    = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_sendrecv)    = false_shift_ie_opt
  endif

  if (sml_ion_all_sendrecv) then
     spall(1)%shift_opt(index_shift_ie_all_sendrecv)    = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_all_sendrecv)    = false_shift_ie_opt
  endif

  if (sml_ion_use_isend) then
     spall(1)%shift_opt(index_shift_ie_use_isend)       = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_isend)       = false_shift_ie_opt
  endif

  if (sml_ion_use_rsend) then
     spall(1)%shift_opt(index_shift_ie_use_rsend)       = true_shift_ie_opt
  else
     spall(1)%shift_opt(index_shift_ie_use_rsend)       = false_shift_ie_opt
  endif

  do i=2,ptl_nsp_max
     spall(i)%shift_opt(:) = spall(1)%shift_opt(:)
  enddo

  !-------------- neutral -------------------------------------------------

  if(sml_neutral) then

     !limiter setup
     if(sml_limiter==1) then
        call limiter_setup
        call check_point('after limiter setup')
     endif

     if(neu_theta_x== -99D0)  neu_theta_x= sml_2pi - acos(  (eq_x_r-eq_axis_r)/dsqrt((eq_x_r-eq_axis_r)**2+(eq_x_z-eq_axis_z)**2)  ) ! X-point theta (radian)
     if(neu_peak_theta== -99D0) neu_peak_theta = neu_theta_x
     if(sml_mype==0)print *,'neu_theta_x==',neu_theta_x,eq_x_r,eq_x_z,eq_axis_r
     neu_ionize2_psi=neu_ionize2_psi * eq_x_psi

     neu_theta_edge(1)=neu_theta_x-neu_delta_theta_lowerx1
     neu_theta_edge(2)=neu_theta_x+neu_delta_theta_lowerx2

     neu_t0=neu_t0*sml_ev2j ! eV to Joule
     neu_base_den=neu_n0
     neu_mfp=neu_mfp0
     if(sml_mype==0 .and. neu_col_mode/=0) then
        print *, 'vth_collision=', sqrt(2D0*neu_t0/ptl_mass(1))
        ! n_e0 should be 'col_denb_edge/2 , otherwise change_neu_mfp routine should be modified.
        if(sml_mype==0) print *,'mfp =' , neu_mfp
        print *, 'theta_x (deg) =' , neu_theta_x*180D0/sml_pi
     endif

     neu_dt=0.002D0*sml_dt/sqrt(neu_temp0*1D-3)
     if(sml_mype==0) print *,'neu_dt:',sml_mype,sml_en_order,eq_axis_r,ptl_mass(1),neu_dt/sml_tran,0.005D0*sml_dt/sqrt(neu_temp0*1D-3)/sml_tran

     neu_psi_edge=neu_psi_edge*eq_x_psi
     neu_grid_max_psi=neu_grid_max_psi*eq_x_psi
     neu_grid_min_psi=neu_grid_min_psi*eq_x_psi
     neu_grid_dtheta=sml_2pi/real(neu_grid_mtheta)
     neu_grid_dpsi=(neu_grid_max_psi-neu_grid_min_psi)/real(neu_grid_mpsi-1)
     neu_grid_dpsi2=(sml_outpsi-sml_inpsi)/real(neu_grid_mpsi-1)
     neu_dtheta=sml_2pi/real(neu_mtheta)
     neu_sep_dtheta=sml_2pi/real(neu_sep_mtheta)

     neu_ion_events=0; neu_cx_events=0; neu_ion_events_weights=0D0; neu_cx_events_weights=0D0
     if(sml_mype==0)print *,'neu_grid_min_psi==',neu_grid_min_psi,eq_x_psi
     if((neu_col_mode == 2) ) then
        call check_point('before neutral2_setup')
        call neutral2_setup
     endif

#ifdef DEGAS2
     if((neu_col_mode == 2) .and. (neu_mode2_period < neu_ion_period)) then
        if(sml_mype==0) then
           print *, 'Should have neu_mode2_period >= neu_ion_period'
           print *, 'Resetting neu_mode2_period = ', neu_ion_period
        endif
        neu_mode2_period=neu_ion_period
     endif
#endif

     allocate(neu_weight_sum_lost(sml_nthreads))
     neu_weight_sum_lost=0D0


     neu_weight_lost_max = neu_lost_rate_max * sml_dt * sml_f_source_period ! neutral3 only
  endif


!-------------------------pre processing ------------------------

  ! initialize


! file description
#ifndef ADIOS
  if(sml_mype==0) then
     open(unit=22, file='fort.tracer', status='replace')
     write(22,'(a400)') '# [1]time [2]r [3]z [4]phi [5]rho_|| [6] E_kin(ev) [7] E_tot(ev) [8]pitch [9]psi [10]canonical P [11]mu [12]weight(df) [13-15] derivs [16] species'
     close(22)
     open(unit=71, file='fort.gam',status='replace')
     write(71,*) '# [1] Time(tran) [2] pot0 [3] dpot [4] idensity0 [5] edensity0'
     !close(71)
     open(unit=72, file='fort.fenergy',status='replace')
     write(72,*) '# [1] Time(tran) [2] Field Energy(J)'
     !close(72)
     write(30,*) '#[1] psi [2] density [3] flow_toro [4] flow_polo [5] flow_para [6] flow_ExB_pol [7] flow_raidal'

     write(50,*)

     write(39,*) '# [1]psi [2] dpdp  [3] E-field at midplane' ! E-field diagnosis
     write(300,*) '# special simulation 3 : [1] energy [2] pitch [3] mass [4] region [5] ptl_index'
     write(301,*) '# specical 3 - reamined particle : [1] K_|| [2] K_perp [3] energy [4] pitch [5] ptl_index'
     write(299,*) '# specical 3 - lost particel : [1] K_|| [2] K_perp [3] energy [4] pitch [5] ptl_index'

     write(103,*) '#bounce setup : '


     open(unit=987,file='fort.pot_zero',status='replace')
     write(987,*) '# r, phi00, den00, grid%psi(grid%itheta0(i))'
     close(987)
  endif
#endif
  if(sml_mype==0) then
     open(unit=10,file='units.txt',status='replace')
     write(10,*) 'sml_totalpe', sml_totalpe
     write(10,*) 'sml_dt(s)=', sml_dt
     write(10,*) 'sml_dt(tau)=',sml_dt/sml_tran
     write(10,*) 'sml_dt(normalized)=',sml_dt
     write(10,*) 'sml_tran(s)=',sml_tran
     write(10,*) 'ptl_num(total)=',ptl_num*sml_totalpe
     write(10,*) 'eq_axis_r(m)=', eq_axis_r
     write(10,*) 'eq_axis_z(m)=', eq_axis_z
     write(10,*) 'eq_axis_B(T)=', eq_axis_b
     write(10,*) 'psi_x(code unit)=', eq_x_psi
     write(10,*) 'psi_x(MKS)=', eq_x_psi
!     write(10,*) 'neu_mean_free_path(m)', neu_mfp
!     write(10,*) 'neu_den (mks)', neu_n0
!     write(10,run_parameters)
     close(10)
  endif
  ! bounce setup
!  ptl_nbounce=0D0
  if(.not. sml_concentric .and. sml_bounce/=0) then
     if(sml_mype==0) print *, 'bounce setup'
     call bounce_setup
  endif

  if(sml_mype==0) print *, 'get_mid_r setup'
  call get_mid_r_setup
  if(sml_mype==0) print *, 'diag_efld setup'


! Smooth initialization
  if(sml_tavg_on==1) then
     smooth00%mode=0
     smooth00%n=1 ! default
     smooth00%type=1 ! default
  else
     smooth00%mode=smooth_mode_in
     smooth00%n=smooth_n_in
     smooth00%type=smooth_type_in
  endif
  call smooth_init(smooth00)

  !need for sml_tavg_on/=0 and diagnosis
  smooth0L%mode=smooth_mode_in
  smooth0L%n=smooth_n_in
  smooth0L%type=smooth_type_in
  call smooth_init(smooth0L)

  ! for high mode smoothing
  smoothH%mode=smooth_H_mode_in
  smoothH%n=smooth_H_n_in
  smoothH%type=smooth_H_type_in
  call smooth_init(smoothH)

  ! for diagnosis smoothing
  smoothdiag%mode=smooth_diag_mode_in
  smoothdiag%n=smooth_diag_n_in
  smoothdiag%type=smooth_diag_type_in
  call smooth_init(smoothdiag)

  ! setup restore temperature module

  !-------------CHECKING INPUT VAILIDITY----------------------------------------------
  if(sml_totalpe*ptl_maxnum > huge(spall(1)%maxgid)) then
     if(sml_mype==0) print *, 'You need 64-bit GID variable for sml_total_pe(',sml_totalpe,') X ptl_maxnum(',ptl_maxnum,').',huge(spall(1)%maxgid)
     stop
  endif

  if(sml_mype==0) then
     call write_units
  endif
contains
  subroutine write_units
    open(unit=15,file='units.m',status='replace')
    write(15,*) 'sml_totalpe =', sml_totalpe, ';'
    write(15,*) 'sml_dt =', sml_dt, ';'
    write(15,*) 'sml_tran=',sml_tran, ';'
    write(15,*) 'vth    =', sqrt(eq_tempi_v1*sml_e_charge/spall(1)%mass), ';'
    write(15,*) 'ptl_num=',ptl_num, ';'
    write(15,*) 'eq_axis_r=', eq_axis_r, ';'
    write(15,*) 'eq_axis_z=', eq_axis_z, ';'
    write(15,*) 'eq_axis_b=', eq_axis_b, ';'
    write(15,*) 'eq_tempi_v1=',eq_tempi_v1, ';'
    write(15,*) 'ptl_ion_mass_au =', ptl_mass_au, ';'
    write(15,*) 'ptl_ion_charge_eu=', ptl_charge_eu, ';'
!    write(15,*) 'ptl_i_num =', ptl_i_num, ';'
!    write(15,*) 'ptl_elec_mass_au =', ptl_elec_mass_au, ';'
!    write(15,*) 'ptl_e_num =',ptl_e_num, ';'
    write(15,*) 'diag_1d_period = ',diag_1d_period, ';'
    write(15,*) 'neu_col_period = ',neu_col_period, ';'
    write(15,*) 'psi_x=', eq_x_psi, ';'
    write(15,*) 'eq_x_z=', eq_x_z, ';'
    write(15,*) 'eq_x_r=', eq_x_r, ';'
    write(15,*) 'sml_wedge_n=', sml_wedge_n, ';'
#ifdef XGC1_EM
    write(15,*) 'v_a=', eq_axis_b/sqrt(sml_mu0*eq_den_v1*ptl_mass(1)),';'
    write(15,*) 'tau_a=', eq_axis_r*sqrt(sml_mu0*eq_den_v1*ptl_mass(1))/eq_axis_b, ';'
    write(15,*) 'beta=', sml_mu0*eq_den_v1*eq_tempi_v1*sml_e_charge/eq_axis_b**2 ,';'
#endif

    close(15)
  end subroutine write_units

end subroutine post_process_input

subroutine read_input
  use sml_module
  use ptl_module
  use itp_module
  use eq_module
  use col_module
  use col_f_module, only : col_f_nthreads
  use diag_module
  use bnc_module
  use grid_class
  use smooth_module
  use pol_decomp_module
  use src_module
  use rad_module
  use perf_monitor
  use neu_module
  use lim_module ! should be merged into neu_module
  use input_module
  use f0_module
  use lbal_module

  implicit none
  include 'mpif.h'

  logical exist
  logical first, comment
  integer :: i,j,ij,ierror,sp
  integer :: input_lth, name_idx
  integer, parameter :: nspace=2
  real (kind=8) :: vol, meanr, r1,r2,z1,z2
  real (kind=8) :: v_the,psi
  character (len=20) :: ptlfile, varname

  namelist /sml_param/ sml_restart, sml_bp_mult, sml_bt_mult,&
       sml_special, sml_machine,sml_mstep, sml_en_order_kev,sml_dt,sml_deltaf,sml_deltaf_elec,&
       sml_concentric,sml_krook_rate,&
       sml_deltaf_f0_mode,sml_marker_temp_factor, sml_marker_min_temp,&
       sml_marker_temp_factor2, sml_marker_temp_factor3, sml_low_mu_fill_population, &
       sml_flat_marker, sml_flat_marker_decay_start1, sml_flat_marker_cutoff1, sml_flat_marker_width1, &
       sml_flat_marker_decay_start2, sml_flat_marker_cutoff2, sml_flat_marker_width2, &
       sml_restart_write_period,sml_monte_num,sml_inpsi,sml_outpsi,&
       sml_bounce,sml_bounce_zero_weight,sml_co_cur_bt, sml_invert_B,&
       sml_read_eq,sml_limiter,&
       sml_bd_min_r,sml_bd_max_r,sml_bd_min_z,sml_bd_max_z,&
       sml_relax_dist, sml_relax_dist_num, &
       sml_f0_grid, sml_f0_grid_lbal_period, &
       sml_f0_grid_min_ptl_imbal, sml_f0_grid_max_ptl_imbal, sml_f0_grid_init_ptl_imbal, &
       sml_f0_grid_alpha, sml_f0_grid_alpha_start, sml_f0_nmu_decomp,& 
       sml_symmetric_f, sml_symmetric_f0g, &
       sml_f_source_period,sml_no_fp_in_f, &
       sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e, &
       sml_node_file, sml_ele_file, sml_bfollow_file,&
       sml_turb_efield, sml_00_efield, sml_turb_poisson, sml_00_poisson, sml_turb_poisson_n0_only, sml_gradpsi_limit, &
       sml_bfollow, &
       sml_tavg_factor, sml_tavg_on,&
       sml_electron_on, sml_electron_hyb, sml_nphi_total, sml_wedge_n, sml_fem_matrix, sml_smtrc_ad_elec, &
       sml_nhybrid,&
       sml_grid_nrho, sml_rhomax, &
       sml_dwdt_exb_only, sml_extra_dwdt, &
       sml_guess_table_size,&
       sml_initial_deltaf_noise,&
       sml_zero_inner_bd,&
       sml_mode_select_on, sml_mode_select_n, sml_mode_initial_n,&
       sml_bd_ext_delta1, sml_bd_ext_delta2, sml_bd_ext_delta3, sml_bd_ext_delta4, &
       sml_bd_ext_delta1H, sml_bd_ext_delta2H, sml_bd_ext_delta3H, sml_bd_ext_delta4H, &
       sml_bd_ext_near_wall, &
       sml_input_file_dir,&
       sml_dwdt_fix_bg, sml_add_pot0, sml_add_pot0_file,&
       sml_replace_pot0, &
       sml_flat_electron_density,&
       sml_suppress_weight_growth, sml_weight_max,&
       sml_restore_temp, sml_restore_temp_period,sml_restore_temp_avg_scale, sml_use_simple00,&
       sml_00_npsi, &
       sml_iter_solver, sml_iter_solver_niter,sml_iter_solver_accel,&
       sml_fsa_solver,sml_poisson_solver_type, &
       sml_simple00_nsmth, &
       sml_max_mat_width, &
       sml_bd_Te_mode, sml_bd_Te_width, sml_zero_out_total_charge, &
       sml_sheath_mode, sml_sheath_init_pot_factor, sml_sheath_adjust, sml_sheath_adjust_factor, &
       sml_exclude_private,sml_rgn1_pot0_only, &
       sml_use_pade, &
       sml_pol_decomp, sml_pol_decomp_simple, sml_max_imbalance, &
       sml_initial_flow, &
       sml_ncycle_half, &
       sml_source, &
       sml_radiation, &
       sml_neutral, sml_neutral_start_step, sml_neutral_use_ion_loss, &
       sml_gpu_ratio, sml_eta, sml_old_grad_perp, &
       sml_perm_gpu_freq, &
       sml_sort_gpu_freq, &
       sml_00_xz_up, &
       sml_ignore_drift_near_wall, sml_ignore_drift_r0, sml_ignore_drift_z0, &
       sml_ignore_drift_slope1, sml_ignore_drift_slope2, &
       sml_elec_shift_max_nthreads, sml_elec_use_alltoall, sml_elec_use_hs_barrier0, &
       sml_elec_large_limit, sml_elec_handshake, sml_elec_use_hs_barrier1, sml_elec_use_sendrecv, &
       sml_elec_all_sendrecv, sml_elec_use_isend, sml_elec_use_rsend, &
       sml_ion_shift_max_nthreads, sml_ion_use_alltoall, sml_ion_use_hs_barrier0, &
       sml_ion_large_limit, sml_ion_handshake, sml_ion_use_hs_barrier1, sml_ion_use_sendrecv, &       
#ifdef XGC1_EM
       sml_use_ts_solver, sml_nrk, &
       sml_hyb_alfven_test, sml_hyb_tearing_test, &
       sml_hyb_linear, sml_hyb_ion_on, sml_hyb_bgraddpe_on, &
       sml_lumped_mass, &
       sml_hyb_debug, sml_use_scaling, sml_n1_diff_coef_perp, &
       sml_n1_diff_coef_para, &
#endif
       sml_ion_all_sendrecv, sml_ion_use_isend, sml_ion_use_rsend
  namelist /ptl_param/  ptl_mass_au,ptl_charge_eu,ptl_num,&
       ptl_e_mass_au, ptl_e_charge_eu, ptl_e_num, ptl_lost_nummax, ptl_maxnum,ptl_e_maxnum, &
       ptl_special_r,ptl_special_z,ptl_special_phi,&
       ptl_special_en_ev,ptl_special_pitch


  namelist /eq_param/ eq_filename,&
       eq_set_x2, eq_x2_r, eq_x2_z, &
       eq_den_shape, &
       eq_den_v1, eq_den_v2, eq_den_v3,&
       eq_den_x1, eq_den_x2, eq_den_x3,&
       eq_tempi_shape, eq_flowi_shape, &
       eq_tempi_v1, eq_tempi_v2, eq_tempi_v3,&  ! electron volt
       eq_tempi_x1, eq_tempi_x2, eq_tempi_x3,&
       eq_flowi_v1, eq_flowi_v2, eq_flowi_v3,&
       eq_flowi_x1, eq_flowi_x2, eq_flowi_x3,&
       eq_tempe_shape, eq_flowe_shape, &
       eq_tempe_v1, eq_tempe_v2, eq_tempe_v3,&  ! electron volt
       eq_tempe_x1, eq_tempe_x2, eq_tempe_x3,&
       eq_flowe_v1, eq_flowe_v2, eq_flowe_v3,&
       eq_flowe_x1, eq_flowe_x2, eq_flowe_x3,&
       eq_zeff_shape,                        &
       eq_zeff_v1,  eq_zeff_v2,  eq_zeff_v3, &
       eq_zeff_x1,  eq_zeff_x2,  eq_zeff_x3, &
       eq_den_file, eq_tempi_file, eq_tempe_file, &
       eq_flowi_file, eq_flowe_file, eq_zeff_file
       

  namelist /col_param/ col_mode, col_period, col_moving_frame,&
       col_pin, col_pout, &
       col_2_pin, col_2_pout,&
       col_2_iteration_method, &
       col_imp_on, col_imp_charge, col_imp_mass,&
       col_den_imp_edge,col_den_imp_out,col_den_imp_ped_c, col_den_imp_ped_width,&
       col_varying_bg, col_vb_period, col_vb_pin, col_vb_pout, &
       col_vb_mtheta, col_vb_m, &
       col_en_col_on, &
       col_accel, col_accel_n, &
       col_accel_pin1, col_accel_pout1, col_accel_factor1, &
       col_accel_pin2, col_accel_pout2, col_accel_factor2, &
       col_3_npsi, col_3_ntheta, col_3_nvr, col_3_nvz, &
       col_3_npsi_sol, col_3_ntheta_sol, &
       col_f_start, col_f_nthreads

  namelist /diag_param/ diag_tracer_period, diag_tracer_n, diag_tracer_sp, &
       diag_particle_mod, diag_particle_period, &
       diag_1d_on, diag_1d_period, diag_1d_npsi, &
       diag_1d_pin, diag_1d_pout, diag_tavg_on, &
       diag_omid_on, &
       diag_3d_on, diag_3d_period, diag_3d_more, &
       diag_f0_period, diag_f0_g_on, diag_f0_df_on, &
       diag_neu_pin, diag_neu_pout, diag_neu_npsi, &
       diag_eflux_on, diag_1d_ne, diag_1d_emin, diag_1d_emax, &
       diag_heat_on, diag_heat_nsection, diag_heat_nr, diag_heat_nz, diag_heat_npsi,&
       diag_heat_rmax1, diag_heat_rmin1, diag_heat_zmax1, diag_heat_zmin1, &
       diag_heat_rmax2, diag_heat_rmin2, diag_heat_zmax2, diag_heat_zmin2, &
       diag_heat_rmax3, diag_heat_rmin3, diag_heat_zmax3, diag_heat_zmin3

  namelist /smooth_param/ smooth_mode_in, smooth_n_in, smooth_type_in, &
       smooth_H_mode_in, smooth_H_n_in, smooth_H_type_in, &
       smooth_diag_mode_in, smooth_diag_n_in, smooth_diag_type_in,&
       smooth_r1_d0_in, smooth_r1_n_in, smooth_r1_type_in, &
       smooth_nearx_d0,smooth_nearx_dr,smooth_nearx_dz, smooth_nearx_nsmth, &
       smooth_pol_efield, smooth_pol_width, smooth_rad_efield, &
       smooth_pol_d0, smooth_grid_rad_width, smooth_grid_rad_sigma

  namelist /src_param/ src_narea, src_narea_e, src_period, src_niter, src_nsubsection, &
       src_pin1, src_pout1, src_decay_width1, src_heat_power1, src_torque1, &
       src_pin2, src_pout2, src_decay_width2, src_heat_power2, src_torque2, &
       src_pin3, src_pout3, src_decay_width3, src_heat_power3, src_torque3, &
       src_pin4, src_pout4, src_decay_width4, src_heat_power4, src_torque4, &
       src_pin1_e, src_pout1_e, src_decay_width1_e, src_heat_power1_e, src_torque1_e, &
       src_pin2_e, src_pout2_e, src_decay_width2_e, src_heat_power2_e, src_torque2_e, &
       src_pin3_e, src_pout3_e, src_decay_width3_e, src_heat_power3_e, src_torque3_e, &
       src_pin4_e, src_pout4_e, src_decay_width4_e, src_heat_power4_e, src_torque4_e

  namelist /rad_param/ rad_filename, rad_impurity_fraction, rad_use_zeff_profile, rad_use_fix_charge, rad_z

  namelist /mon_param/ mon_flush_count, mon_flush_freq

  namelist /neu_param/ neu_sepfile,neu_limfile, neu_col_mode,neu_adjust_n0,&
       neu_ionize_mode,neu_ionize2_psi,&
       neu_cx_period,neu_ion_period,neu_col_period,&
       neu_update_elec,&
       neu_start_time, neu_recycle_rate,&
       neu_varying_mfp,&
       neu_t0,neu_n0,neu_delta_n,neu_delta_theta,neu_mfp0,&
       neu_theta_x, neu_temp_factor,&
       neu_monte_num,neu_mode2_period,neu_mpol,&
       neu_temp0,neu_dt,neu_istep_max,neu_psi_edge,neu_grid_max_psi,neu_peak_theta,neu_grid_mode,&
       neu_lost_rate_max, &
#ifdef DEGAS2
       neu_grid_min_psi,neu_ionize_elec_cooling_ev,neu_flow_on &
       neu_node_file,neu_ele_file,neu_diag_pls,neu_diag_neu,neu_diag_cx
#else
       neu_grid_min_psi,neu_ionize_elec_cooling_ev,neu_flow_on
#endif



  namelist /f0_param/ f0_nmu, f0_nvp, f0_smu_max, f0_vp_max, &
       f0_col_change_weight, f0_f_correction
#ifdef ADIOS
  namelist /adios_param/ adios_stage_restart, adios_stage_restartf0, &
       adios_stage_f0, adios_stage_field3d, adios_stage_particle, &
       adios_stage_dpot, adios_stage_streamer, adios_mfilesys
#endif
  ! 0. reads input file into input file string
#ifndef USE_OLD_READ_INPUT
  call read_input_file
#endif

  call check_point('sml_param')

  ! 1.  set default value to phase1 variables

  sml_nphi_total=16 ! default 16
  sml_wedge_n=1
  sml_input_file_dir='./'
  sml_restart=.false.
  sml_bp_mult=1D0
  sml_bt_mult=1D0
  sml_grid_nrho=6
  sml_rhomax=1D-2  ! 1 cm
  sml_special=0  ! 0 : nomral simulation
  ! 1 : single particel simulation
  ! 2 : undefined
  ! 3 : velocity hole calculation
  sml_machine=0  !                0 : circular
  ! macine number  1: D3D
  !                2: NSTX
  !                3: alcator-CMOD
  !                4: ASDEX-U
  sml_mstep=3000          ! total time step. (simulation time) = dt*mstep
  sml_en_order_kev = 0.2D0  ! energy order in KeV

  sml_dt=0.001D0           !unit of tran --> code unit
  sml_krook_rate=0D0 !unit of 1/sml_dt
  sml_deltaf=.false.             !0: full f, 1: delta f - delta f simulation is not verified yet
  sml_deltaf_elec=.true.
  sml_deltaf_f0_mode=0
  sml_marker_temp_factor=1D0   ! Marker temperature ratio -- (Marker when weight =1 )
  sml_marker_min_temp = 10D0   ! 10 eV minium temperature of marker particle
  sml_flat_marker=.false.
  sml_flat_marker_cutoff1=4D0
  sml_flat_marker_decay_start1=2.5D0
  sml_flat_marker_width1=0.5D0
  sml_restart_write_period=10000000
  sml_inpsi=0.9             ! eq_x_psi unit --> code unit
  sml_outpsi=1.05           ! eq_x_psi unit --> code unit
  ! if current I_T is in -phi direction, sml_minusB=0
  ! if current I_T is in the phi direction, sml_minusB=1
  ! phi is the cylindrical angle
  sml_co_cur_bt=.true.     ! true if cocurrent, false if counter current
  sml_invert_B=.false.        !sml_invert_B=.true. , vec B -> - vec B

  sml_electron_on=.false.
  sml_electron_hyb=.false.

  sml_nhybrid=2
  sml_ncycle_half=50

  sml_relax_dist=0
  sml_relax_dist_num=1000

  sml_turb_efield=.true.
  sml_turb_poisson=.true.
  sml_turb_poisson_n0_only=.false.
  sml_00_efield= .true.
  sml_00_poisson=.true.
  !rh |grad(psi)| threshold for FD (R,Z)-gradient
  sml_gradpsi_limit=1D-3

  sml_00_npsi = 150 ! 150 data point for 00 mode

  sml_tavg_on=0
  sml_tavg_factor= 1D-3

  sml_fem_matrix=1
  sml_smtrc_ad_elec=0
  sml_dwdt_exb_only=.false.
  sml_extra_dwdt=.false.
  sml_dwdt_fix_bg=.false.
  sml_guess_table_size=1000
  sml_initial_deltaf_noise=1D-3
  sml_zero_inner_bd=0
  sml_mode_select_on=0
  sml_mode_select_n=1
  sml_mode_initial_n=0

  sml_f0_grid=.false.            ! turn on grid alpha
  sml_f0_grid_lbal_period=0
  sml_f0_grid_min_ptl_imbal=1.0D0
  sml_f0_grid_max_ptl_imbal=2.0D0
  sml_f0_grid_init_ptl_imbal=-1.0D0
  sml_f0_grid_alpha=0D0
  sml_f0_grid_alpha_start=1
  sml_f0_nmu_decomp=1
  sml_f_source_period=10
  sml_symmetric_f0g=.false.
  sml_symmetric_f=.false.
  sml_no_fp_in_f=.false.
  sml_f0_1_Lt=6.92D0
  sml_f0_1_Lt_e=6.92D0
  sml_f0_1_Ln=sml_f0_1_Lt/3.114D0
  sml_bounce_zero_weight=0

  sml_bd_ext_delta1=0.01 ! inner boundary extension for poisson eq
  sml_bd_ext_delta2=0.02 ! inner boundary extension for charge zeroing out
  sml_bd_ext_near_wall=0D0

  sml_add_pot0=0
  sml_add_pot0_file='pot0.dat'
  sml_replace_pot0=.false.
  sml_flat_electron_density=0
  sml_source=.false.
  sml_radiation =.false.
  sml_neutral=.false.
  sml_neutral_start_step=0
  sml_neutral_use_ion_loss=.false.
  sml_pol_decomp=.false.
  sml_pol_decomp_simple=.false.
  sml_max_imbalance=1.1D0
  sml_rpl_on=.false.
  sml_00_xz_up=.false.

  sml_suppress_weight_growth=.false.  ! if this value is true then weight of deltaf will be between -sml_weight_max and sml_weight_max. default value is .false.
  sml_weight_max=10.   ! maximum weight(= phase(6)) for delta-f simulation. default value is 10.

  ! restore temperature gradient
  sml_restore_temp = .false.  ! restore temperature gradient. Work for delta-f=1, sml_deltaf_f0_mode=-1, and rk2 method only
  sml_restore_temp_period=5 ! how frequently restore temperature
  sml_restore_temp_avg_scale=0.01 ! time average scale - dt base, default=0.01 (means 100 dt scale)

  sml_bd_Te_mode=0
  sml_bd_Te_width=0.01D0

  sml_use_pade=.false.
  sml_initial_flow=.false.
  ! for concentric circular, 00-mode poission equation is solved analytically when we use poisson_two_solvers
  ! it has smoothing routine in it and the number of smoothing operation of charge density is give by sml_simple00_nsmth
  sml_use_simple00=.false.
  sml_simple00_nsmth=0
  sml_iter_solver=.false.
  sml_iter_solver_niter=10
  sml_iter_solver_accel=1D0
  sml_fsa_solver=.false.
  sml_poisson_solver_type=-1  ! 0 for old solver, 1 & 2 for new solver (linear and non-linear)
  sml_zero_out_total_charge=.true.
  sml_gpu_ratio = -1.0
  sml_perm_gpu_freq = -1
  sml_sort_gpu_freq = -1

  sml_ignore_drift_near_wall =.false.
  sml_ignore_drift_r0=1.5575D0
  sml_ignore_drift_z0=-1.1771D0
  sml_ignore_drift_slope1=-0.15D0
  sml_ignore_drift_slope2=0.3D0

#ifdef XGC1_EM
  ! resistivity
  sml_eta=0
  ! fake diffusion coeficients for implicit solver
  sml_n1_diff_coef_perp=0D0
  sml_n1_diff_coef_para=0D0
  ! Use implicit PETSc time stepper
  sml_use_ts_solver = .false.
  sml_nrk=2
  ! use mass solve (true) or 1/area (false)
  !rh sml_mass_solve=.true.
  sml_lumped_mass=.false.
  sml_hyb_linear=.false.
  sml_hyb_alfven_test=.false.
  sml_hyb_tearing_test=.false.
  sml_hyb_bgraddpe_on=.false.
  sml_hyb_ion_on=.true.
  sml_hyb_debug=0
  sml_use_scaling=.true.
#endif
  ! use FE perp. gradient (true) or FD gradient (false)
  sml_old_grad_perp=.true.


  ! shift_ie option defaults; also defined in shift_ie routine
  sml_elec_shift_max_nthreads = -1
  sml_elec_use_alltoall = .false.
  sml_elec_use_hs_barrier0 = .false.
  sml_elec_large_limit = 2048
  sml_elec_handshake = .false.
  sml_elec_use_hs_barrier1 = .false.
  sml_elec_use_sendrecv = .false.
  sml_elec_all_sendrecv = .false.
  sml_elec_use_isend = .false.
  sml_elec_use_rsend = .true.

  sml_ion_shift_max_nthreads = -1
  sml_ion_use_alltoall = .false.
  sml_ion_use_hs_barrier0 = .false.
  sml_ion_large_limit = 2048
  sml_ion_handshake = .false.
  sml_ion_use_hs_barrier1 = .false.
  sml_ion_use_sendrecv = .false.
  sml_ion_all_sendrecv = .false.
  sml_ion_use_isend = .false.
  sml_ion_use_rsend = .true.

  ! read sml_parameters
#ifdef USE_OLD_READ_INPUT
  !USE_OLD_READ_INPUT can hurt performance but the other option needs F2003 capability (path95 cannot handle it)
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"sml_param")
  read(input_string(name_idx:len_trim(input_string)),nml=sml_param)
#endif

#ifdef XGC1_EM
  if (sml_hyb_debug==1) then
    sml_use_scaling=.false.
  endif
  if ( ((.not. sml_use_ts_solver) .or. sml_hyb_ion_on) .and. sml_nrk .ne.2) then
    if (sml_mype==0) print *, 'Your input parameters require sml_nrk=2!!!'
    stop
  endif
#endif

  sml_bd_ext_delta3=sml_bd_ext_delta1 ! outer boundary extension for poisson eq
  sml_bd_ext_delta4=sml_bd_ext_delta2 ! outer boundary extension for charge zeroing out

  ! 2 for core simulation, 1 for edge simulation
  if(sml_outpsi>1D0) then
     sml_bounce=1
  else
     sml_bounce=2
  endif

  if (sml_machine==0) then
     sml_concentric=.false.
     sml_read_eq=.true.
     sml_limiter=0
  else
     sml_concentric=.false.
     sml_read_eq=.true.
     sml_limiter=0
  endif

  if(sml_machine==1) then
     ! simulation boundary for D3D
     sml_bd_min_r=0.7  ! inner boundary
     sml_bd_max_r=2.5  ! outer boundary
     sml_bd_min_z=-1.8 ! lower boundary
     sml_bd_max_z=1.8  ! upper boundary
  elseif(sml_machine==2) then
     ! simulation boundary for NSTX
     sml_bd_min_r=0.15 ! inner boundary
     sml_bd_max_r=1.7  ! outer boundary
     sml_bd_min_z=-1.8 ! lower boundary
     sml_bd_max_z=1.8  ! upper boundary
  elseif(sml_machine==3) then
     ! simulation boundary for alcator-CMOD
     sml_bd_min_r=0.4  ! inner boundary
     sml_bd_max_r=1.   ! outer boundary
     sml_bd_min_z=-.7  ! lower boundary
     sml_bd_max_z=0.65 ! upper boundary
  elseif(sml_machine==4) then
     ! simulation boundary for ASDEX-U
     sml_bd_min_r=1D0
     sml_bd_max_r=2.25D0
     sml_bd_min_z= -1.5D0
     sml_bd_max_z=1.3
  else
     sml_bd_min_r=0.0001
     sml_bd_max_r=1000.
     sml_bd_min_z=-1000.
     sml_bd_max_z=1000.
  endif

  sml_node_file="neo10.1.node"
  sml_ele_file="neo10.1.ele"
  sml_bfollow_file="neo10.bf.dat"
  sml_bfollow=0
  sml_max_mat_width=100

  sml_exclude_private=.true.
  sml_rgn1_pot0_only= .false.

  sml_sheath_mode=0
  sml_sheath_init_pot_factor=2D0
  sml_sheath_adjust=.false.
  sml_sheath_adjust_factor=1D0


  sml_marker_temp_factor2=sml_marker_temp_factor
  sml_marker_temp_factor3=sml_marker_Temp_factor
  sml_low_mu_fill_population=0D0


  sml_flat_marker_cutoff2     =sml_flat_marker_cutoff1
  sml_flat_marker_decay_start2=sml_flat_marker_decay_start1
  sml_flat_marker_width2      =sml_flat_marker_width1
  ! read sml_parameters
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"sml_param")
  read(input_string(name_idx:len_trim(input_string)),nml=sml_param)
#endif

  sml_bd_ext_delta1H=sml_bd_ext_delta1
  sml_bd_ext_delta2H=sml_bd_ext_delta2
  sml_bd_ext_delta3H=sml_bd_ext_delta3
  sml_bd_ext_delta4H=sml_bd_ext_delta4

  ! read sml_parameters
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=sml_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"sml_param")
  read(input_string(name_idx:len_trim(input_string)),nml=sml_param)
#endif

  ! if sml_f0_grid_lbal_period > 0, then it is the period of
  ! f_source calls; turn this into the period in terms of istep, for
  ! convenience
  if (sml_f0_grid_lbal_period > 0) then
     f0_grid_load_balance = .true.
     sml_f0_grid_lbal_period = max(sml_f0_grid_lbal_period*sml_f_source_period,2)
     if (sml_f0_grid_init_ptl_imbal < 1.0D0) then
        sml_f0_grid_init_ptl_imbal = sml_max_imbalance
     endif
  else
     f0_grid_load_balance = .false.
     sml_f0_grid_lbal_period = 1
     sml_f0_grid_init_ptl_imbal = 1.0D0
  endif

#ifdef XGC1_EM
  if ( (sml_hyb_alfven_test .and. sml_hyb_tearing_test) .or. &
       (sml_hyb_alfven_test .and. sml_hyb_linear) .or. &
       (sml_hyb_tearing_test .and. sml_hyb_linear) ) then
    print *, 'Select only one of sml_hyb_alfven_test, sml_hyb_tearing_test, and sml_hyb_linear'
    stop
  endif
#endif

  ! initialize collision load balancing algorithm
  call update_ptl_constraint_init(sml_f0_grid_min_ptl_imbal, &
                                  sml_f0_grid_max_ptl_imbal, &
                                  sml_f0_grid_init_ptl_imbal )

  !---------------------------------------------------------------
  call check_point('ptl_param')
  ptl_mass_au=2D0  !mass ratio to proton
  ptl_e_mass_au=2D-2
  ptl_charge_eu=1D0  ! charge number
  ptl_e_charge_eu=-1D0

  ! special simulation set definition
  if(sml_special==1) then !single particle simulation
     ptl_maxnum=10
     ptl_num=1
     ptl_special_r=1.95D0
     ptl_special_z=0.
     print *, ptl_special_r, ptl_special_z
     ptl_special_phi=0D0
     ptl_special_en_ev=7205
     ptl_special_pitch=0.276D0
  else if(sml_special==3) then
     ptl_num=5000
     ptl_special_n=100
     ptl_special_r=1.95
     ptl_special_z=0.
     ptl_special_phi=0D0
     ptl_special_en_ev=1000D0
  else
     ptl_maxnum=4000
     ptl_num=3000   ! default particle number per process or thread (OpenMP)
  endif

#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=ptl_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"ptl_param")
  read(input_string(name_idx:len_trim(input_string)),nml=ptl_param)
#endif

  ! number of electrons = number of ions
  ptl_e_maxnum=ptl_maxnum
  ptl_e_num=ptl_num


  ptl_lost_nummax=max(ptl_num/10,100)

  if(sml_special==1) then
     if(ptl_num/=1) print *, 'Warning. ptl_num is not 1 for single ptl simulation'
  endif

  sml_monte_num=max(5*ptl_num,10000)

  ! give 20% extra space
  ptl_maxnum=max(ptl_num+ptl_num/5,ptl_maxnum)
  ptl_e_maxnum=max(ptl_e_num+ptl_e_num/5,ptl_e_maxnum)


  !-------------READ Equilbirum file------------------------------
  call check_point('eq_param')

  eq_set_x2=.false.
  eq_x2_r=1.7D0
  eq_x2_z=10D0


  eq_den_shape=1
  eq_den_v1= 5D19   ! m^-3
  eq_den_v2= 0.51D19  ! m^-3
  eq_den_v3= eq_den_v2
  eq_den_x1= 0.5*(sml_inpsi+sml_outpsi) !x_psi unit -> code unit
  eq_den_x2= 0.001  ! x_psi unit -> code unit
  eq_den_x3= 0D0


  !read eq param 1st
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"eq_param")
  read(input_string(name_idx:len_trim(input_string)),nml=eq_param)
#endif

  ! tempi - ion temperature
  eq_tempi_shape=eq_den_shape
  eq_tempi_v1=1D3   !ev
  eq_tempi_v2=0.1D3  !ev
  eq_tempi_v3=eq_tempi_v2
  eq_tempi_x1=eq_den_x1
  eq_tempi_x2=eq_den_x2
  eq_tempi_x3=eq_den_x3


  !read eq param 2nd
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"eq_param")
  read(input_string(name_idx:len_trim(input_string)),nml=eq_param)
#endif

  !electron temperature
  eq_tempe_shape=eq_tempi_shape
  eq_tempe_v1=eq_tempi_v1
  eq_tempe_v2=eq_tempi_v2
  eq_tempe_v3=eq_tempi_v3
  eq_tempe_x1=eq_tempi_x1
  eq_tempe_x2=eq_tempi_x2
  eq_tempe_x3=eq_tempi_x3

  ! initial flow
  eq_flowi_shape=0
  eq_flowi_v1=0D0
  eq_flowi_v2=0D0
  eq_flowi_v3=0D0
  eq_flowi_x1=0D0
  eq_flowi_x2=0D0
  eq_flowi_x3=0D0

  eq_flowe_shape=0
  eq_flowe_v1=0D0
  eq_flowe_v2=0D0
  eq_flowe_v3=0D0
  eq_flowe_x1=0D0
  eq_flowe_x2=0D0
  eq_flowe_x3=0D0

  ! zeffective
  eq_zeff_shape=0
  eq_zeff_v1=1D0
  eq_zeff_v2=1D0
  eq_zeff_v3=1D0
  eq_zeff_x1=1D0
  eq_zeff_x2=0.1D0
  eq_zeff_x3=0D0


#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=eq_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"eq_param")
  read(input_string(name_idx:len_trim(input_string)),nml=eq_param)
#endif

  !-----------------------------------------------------
  ! Collision parameters
  call check_point('col_param')

  col_mode=2  ! 0: off 1 : monte-carlo (non-conserving)
  ! 2: conserving collision (monte-carlo)
  col_period=3  ! collsion each 3 time step
  col_en_col_on=1  ! energy collision on / off
  col_moving_frame=.true.
 
  col_pin= sml_inpsi ! x_psi unit --> will be SI after post_process
  col_pout=sml_outpsi

  col_f_start=1
  col_f_nthreads=sml_nthreads

  col_2_pin=sml_inpsi  ! x_psi unit -> code unit
  col_2_pout=min(sml_outpsi,1.03D0)  ! x_psi unit -> code unit

  col_varying_bg=1  ! change background profile for collision as time goes
  col_vb_period=1  ! period of change
  col_vb_pin=sml_inpsi  ! x_psi unit -> code unit
  col_vb_pout=min(sml_outpsi,1.03D0)  ! x_psi unit -> code unit

  col_imp_on=0        ! impurity collision on/off
  col_imp_charge=4.5  ! impurity charge number
  col_imp_mass=12.    ! impurity mass (AU)
  col_den_imp_edge=eq_den_v1*0.027D0  ! m-3 -> code unit
  col_den_imp_out =eq_den_v2 *0.027D0  ! m^-3 -> code unit
  col_den_imp_ped_c =  eq_den_x1     ! x_psi unit -> code unit
  col_den_imp_ped_width=eq_den_x2 ! x_psi unit -> code unit


  col_accel=.false.    ! artifitial increase of collision
  col_accel_n=2
  col_accel_pin1=sml_inpsi
  col_accel_pout1=sml_inpsi + 0.1*(sml_outpsi-sml_inpsi)
  col_accel_factor1=10
  col_accel_pin2=sml_outpsi - 0.1*(sml_outpsi-sml_inpsi)
  col_accel_pout2=sml_outpsi
  col_accel_factor2=10

  ! Col 3 by E. Yoon
  ! 3: vpic collision (Eulerian NL Fokker-Planck-Landau)
  col_3_npsi   = 64
  col_3_ntheta = 512
  col_3_nvr    = 32
  col_3_nvz    = 64
  col_3_npsi_sol = col_3_npsi/6  !This should be modified as referring to domain ratio
  col_3_ntheta_sol = col_3_ntheta
  !---

  ! read col_parameters
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=col_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"col_param")
  read(input_string(name_idx:len_trim(input_string)),nml=col_param)
#endif

  col_vb_period=max(col_vb_period,col_period)
  if(col_mode==2 .and. col_varying_bg ==1) then
     col_vb_m=col_2_m
     col_vb_mtheta=col_2_mtheta
     col_vb_pin=col_2_pin
     col_vb_pout=col_2_pout
  endif

  ! Col 3 by E. Yoon
  if(col_mode==3) then
     ! In this program, we assume that total number of collision cells are more than total CPU numbers
     ! Check this requirements
     call assert( sml_pe_per_plane .le. col_3_npsi*col_3_ntheta+col_3_npsi_sol*col_3_ntheta_sol, 'Change sml_totalpe > col_3_npsi*col_3_ntheta', ierror)
     call assert( col_3_npsi+col_3_ntheta+col_3_npsi_sol+col_3_ntheta_sol .gt. 0, 'You do NOT want col_mode==3. Please modify input. Program aborted',ierror)
  endif
  !---



  !----------------------------------------------------
  ! diagnosis
  call check_point('diag_param')
  ! tracer : record time history of a given particle
  diag_tracer_period=100
  diag_tracer_n=1  ! index number of traced particle
  diag_tracer_sp=1 ! species index number of traced particle

  ! diag_particle
  diag_particle_mod=0  ! zero for no particle output
  diag_particle_period=1

  ! diag 1d
  diag_1d_on=.true.
  diag_1d_period=10
  diag_1d_npsi=80
  diag_1d_pin=sml_inpsi
  diag_1d_pout=sml_outpsi

  diag_tavg_on=.false.
  diag_omid_on=.false.

  diag_eflux_on=.true.
  diag_1d_ne=100
  diag_1d_emin=0D0

  ! diag 3d
  diag_3d_on=.true.
  diag_3d_more=.false.
  diag_3d_period=-1

  ! diag f0
  diag_f0_period=100
  diag_f0_g_on=.false.
  if(sml_f0_grid) then
    diag_f0_df_on=.true.
  else
    diag_f0_df_on=.false.
  endif

  diag_neu_pin=0.8
  diag_neu_pout=sml_outpsi
  diag_neu_npsi=50

  diag_heat_on=.false.
  diag_heat_nsection=1
  diag_heat_nr=100
  diag_heat_nz=100
  diag_heat_npsi=1000
  diag_heat_rmax1=sml_bd_max_r
  diag_heat_rmin1=sml_bd_min_r
  diag_heat_zmax1=sml_bd_max_z
  diag_heat_zmin1=sml_bd_min_z

  diag_heat_rmax2=sml_bd_max_r
  diag_heat_rmin2=sml_bd_min_r
  diag_heat_zmax2=sml_bd_max_z
  diag_heat_zmin2=sml_bd_min_z

  diag_heat_rmax3=sml_bd_max_r
  diag_heat_rmin3=sml_bd_min_r
  diag_heat_zmax3=sml_bd_max_z
  diag_heat_zmin3=sml_bd_min_z


  ! read diag_parameters
#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=diag_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"diag_param")
  read(input_string(name_idx:len_trim(input_string)),nml=diag_param)
#endif

  if(diag_3d_period==-1) diag_3d_period=diag_1d_period  ! 3d_period is 1d_period

  ! smooth -----------------------------------------------------------------
  call check_point('smooth_param')
  smooth_mode_in=2
  smooth_n_in=30
  smooth_type_in=1  !gaussian type

  smooth_H_mode_in=2
  smooth_H_n_in=3
  smooth_H_type_in=1

  smooth_r1_n_in=0
  smooth_r1_d0_in=0.01
  smooth_r1_type_in=1

  smooth_nearx_nsmth=0
  smooth_nearx_d0=0.01
  smooth_nearx_dr=0.05
  smooth_nearx_dz=0.05

  smooth_pol_efield=.false.
  smooth_pol_width=3  ! Might be a little small for 128 poloidal bins
  smooth_pol_d0=0.1  !! smooth out poloidal variations under 10 cm
  smooth_rad_efield=.false.
  smooth_grid_rad_width=0.75D0
  smooth_grid_rad_sigma=1.7D0

#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=smooth_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"smooth_param")
  read(input_string(name_idx:len_trim(input_string)),nml=smooth_param)
#endif

  smooth_diag_mode_in=smooth_mode_in
  smooth_diag_n_in=smooth_n_in
  smooth_diag_type_in=smooth_type_in

#ifdef USE_OLD_READ_INPUT
  open(unit=14,file='input',action='read')
  read(14,nml=smooth_param)
  close(14)
#else
  name_idx = read_namelist_index(input_string,"smooth_param")
  read(input_string(name_idx:len_trim(input_string)),nml=smooth_param)
#endif

  ! source (heating and torque) ------------------------------------------------------
  if(sml_source ) then
     call check_point('src_param')
     src_narea=1
     src_narea_e=0
     src_period=10
     src_niter=10
     src_nsubsection=10
     ! 1st area
     src_pin1=sml_inpsi
     src_pout1=sml_inpsi+0.01
     src_decay_width1=0.01
     src_heat_power1=0D0  ! Watt
     src_torque1= 0D0     ! N m
     ! 2nd area
     src_pin2=sml_outpsi-0.01
     src_pout2=sml_outpsi
     src_decay_width2=0.01
     src_heat_power2=0D0 ! Watt
     src_torque2=0D0     ! N m
     ! 3rd area
     src_pin3=0D0
     src_pout3=-1D0
     src_decay_width3=0.01
     src_heat_power3=0D0 ! Watt
     src_torque3=0D0     ! N m
     ! 4th area
     src_pin4=0D0
     src_pout4=-1D0
     src_decay_width4=0.01
     src_heat_power4=0D0 ! Watt
     src_torque4=0D0     ! N m

     !electron
     src_pin1_e=sml_inpsi
     src_pout1_e=sml_inpsi+0.01
     src_decay_width1_e=0.01
     src_heat_power1_e=0D0  ! Watt
     src_torque1_e= 0D0     ! N m
     ! 2nd area
     src_pin2_e=sml_outpsi-0.01
     src_pout2_e=sml_outpsi
     src_decay_width2_e=0.01
     src_heat_power2_e=0D0 ! Watt
     src_torque2_e=0D0     ! N m
     ! 3rd area
     src_pin3_e=0D0
     src_pout3_e=-1D0
     src_decay_width3_e=0.01
     src_heat_power3_e=0D0 ! Watt
     src_torque3_e=0D0     ! N m
     ! 4th area
     src_pin4_e=0D0
     src_pout4_e=-1D0
     src_decay_width4_e=0.01
     src_heat_power4_e=0D0 ! Watt
     src_torque4_e=0D0     ! N m



#ifdef USE_OLD_READ_INPUT
     open(unit=14,file='input',action='read')
     read(14,nml=src_param)
     close(14)
#else
     name_idx = read_namelist_index(input_string,"src_param")
     read(input_string(name_idx:len_trim(input_string)),nml=src_param)
#endif
  else
     src_period=10
  endif


  ! impurity radiation
  if(sml_radiation) then
     call check_point('rad_param')


     rad_filename = 'mist_adas_xgc1.dat'
     rad_impurity_fraction = 0.1D0
     rad_use_zeff_profile=.false.
     rad_use_fix_charge=.true.
     rad_z=6D0


#ifdef USE_OLD_READ_INPUT
     open(unit=14,file='input',action='read')
     read(14,nml=rad_param)
     close(14)
#else
     name_idx = read_namelist_index(input_string,"rad_param")
     read(input_string(name_idx:len_trim(input_string)),nml=rad_param)
#endif
  endif


  ! neutral collision ----------------------------------------------------------------
  ! unlike XGC0 all angle is radian in input
  if(sml_neutral .or. sml_radiation) then
     call check_point('neu_param')
     ! neutral collision switch
     ! 1: simple neutral collision - given neutral density, 2: consistent neutral calculation
     neu_col_mode=2
     neu_update_elec=.true.
     neu_adjust_n0=1         ! 0 fixed n0, 1 varying n0 - neu_recycle_rate should be specified.
     neu_ionize_mode=1       ! 1 - limiter hit recycling, 2 - sep flow recycleing
     neu_ionize2_psi=1D0     ! x_psi unit -> code unit

     neu_cx_period=50        ! charge exchange period
     neu_ion_period=50       ! ionization freq. (period)
     neu_col_period=50       ! neutral elastic collision period
     ! Above three will be merged to neu_period

     neu_elastic_col_on=1    ! 1 - turn on neutral elastic collision, 0 -  turn off

     neu_start_time=0*col_vb_period ! Start the neutral effect after some background adjustment
     neu_recycle_rate=0.9D0  ! recycle rate of neutral
     neu_varying_mfp=1       ! 1: varying mean free path  0 : constant mean free path -- for mode 1 only
     neu_t0= eq_tempi_v2     ! ev -> normalized unit, temperature -- for mode 1 only
     neu_n0=5.0D17      ! m^-3 -> nomalized unit , base neutral density
     neu_delta_n=  10.0D0     ! relative neut den over neu_n0 at poloidal peak
     neu_delta_theta=  30D0/180D0*sml_pi  ! poloidal peak witdh (rad)
     ! mean free path = v_0 * tau_ionize = sqrt(2 T_0) / (n_e0 * <sigma V> ionize )
     !neu_mfp0=2D0*sqrt(2D0 * neu_t0) / ( &
     !     eq_den_edge * 0.8D-8 *sqrt(0.75D0*eq_tempe_ev_edge*1D3) &
     !     * (exp(-13.56/eq_tempe_ev_edge/0.75D3)/(1D0+0.01D0*eq_tempe_ev_edge*0.75D3) ) * sml_norm_t &
     !     )
     neu_ionize_elec_cooling_ev=50D0 ! electron cooling energy (eV) per one ionization event
     neu_theta_x = -99D0   ! -99D0 will make neu_theta_x to real x-point angle
     neu_mfp0=0.05         ! cm -> code unit 5cm - initial mean free path

     neu_peak_theta=neu_theta_x  ! poloidal peak theta of neutral birth distribution (rad)
     neu_temp_factor=0.1



     ! neutral 2 only-------------------
     ! neutral 2 gets neutral temp and density consistently with boundary conditions
     neu_monte_num=10*sml_monte_num
     neu_mode2_period=20*col_vb_period
     neu_num=10000
     neu_mpol=10
     neu_temp0=3D0 !neutral temp. boundary condition  eV
     neu_istep_max=20000
     neu_sepfile='sep.dat'
     neu_limfile='lim.dat'

     neu_psi_edge    =sml_outpsi
     neu_grid_max_psi=sml_outpsi  ! minimum psi of neutral grid on which neutral distribution is measured
     neu_grid_min_psi=0.9D0       ! maximum psi of neutral grid on which neutral distribution is measured
     neu_flow_on=0       ! on(1)/off(0) of parallel momentum exchange in plasma-neutral collision process
     neu_grid_mode=0             ! 0: original neutral grid setting, 1: improved neutral grid setting
     neu_delta_theta_lowerx1=20D0/180D0*sml_pi
     neu_delta_theta_lowerx2=20D0/180D0*sml_pi

     !#ifdef DEGAS2
     !
     ! Settings of DEGAS 2 specific variables
     !
     !  neu_node_file="DIII-D_polygon_2.1.node"
     !  neu_ele_file="new_DIII-D_polygon_2.1.ele"
     !  neu_diag_neu=0
     !  neu_diag_pls=0
     !  neu_diag_cx=0
     !#endif


     neu_lost_rate_max=  5D23  ! particle lost per second


#ifdef USE_OLD_READ_INPUT
     open(unit=14,file='input',action='read')
     read(14,nml=neu_param)
     close(14)
#else
     name_idx = read_namelist_index(input_string,"neu_param")
     read(input_string(name_idx:len_trim(input_string)),nml=neu_param)
#endif
  endif



  ! f0_grid settings -----------------------------------------------
  if(sml_f0_grid) then
     call check_point('f0_param')
     f0_col_change_weight=.true.
     f0_f_correction=.false.

     f0_nmu=31
     f0_nvp=15

     !f0_mu0_factor=3D0 ! This does not depend on whether linear mu or sqrt(mu) grid is used
     ! Purpose: if (i_mu==0) then smu=f0_dsmu/f0_mu0_factor else smu=i_mu*f0_dsmu (sqrt(mu))
     !          if (i_mu==0) then mu=f0_dmu/f0_mu0_factor else mu=i_mu*f0_dmu (linear mu grid)

     ! The same for sqrt(mu) and v_perp grid
     f0_smu_max=3.D0
     f0_vp_max=3.D0

  ! read f0_parameters
#ifdef USE_OLD_READ_INPUT
     open(unit=14,file='input',action='read')
     read(14,nml=f0_param)
     close(14)
#else
     name_idx = read_namelist_index(input_string,"f0_param")
     read(input_string(name_idx:len_trim(input_string)),nml=f0_param)
#endif
  endif


  ! performance monitoring -----------------------------------------------
  call check_point('mon_param')
  mon_flush_count = 12
  mon_flush_freq = 0

  if (sml_mype == 0) then
     inquire(FILE="mon_in",EXIST=exist)
     if (exist) then
        open(unit=14,file='mon_in',action='read')
        read(14,nml=mon_param)
        close(14)
     else
#ifdef USE_OLD_READ_INPUT
        open(unit=14,file='input',action='read')
        read(14,nml=mon_param)
        close(14)
#else
        name_idx = read_namelist_index(input_string,"mon_param")
        read(input_string(name_idx:len_trim(input_string)),nml=mon_param)
#endif
     endif

     if (mon_flush_count > 0) then
        mon_flush_freq  = max(1,sml_mstep/mon_flush_count)
     endif
  endif

  call mpi_bcast(mon_flush_count, 1, MPI_INTEGER, 0, sml_comm, ierror)
  call mpi_bcast(mon_flush_freq,  1, MPI_INTEGER, 0, sml_comm, ierror)
#ifdef ADIOS
  ! adios parameters ----------------------------------------------------
  adios_stage_restart = .FALSE.
  adios_stage_restartf0 = .FALSE.
  adios_stage_f0 = .FALSE.
  adios_stage_field3d = .FALSE.
  adios_stage_particle = .FALSE.
  adios_stage_dpot = .FALSE.
  adios_stage_streamer = .FALSE.
  adios_mfilesys = .FALSE.

  inquire(FILE="adios_in",EXIST=exist)
  if (exist) then
     open(unit=14,file='adios_in',action='read')
     read(14,nml=adios_param)
     close(14)
  else
#ifdef USE_OLD_READ_INPUT
     open(unit=14,file='input',action='read')
     read(14,nml=adios_param)
     close(14)
#else
     name_idx = read_namelist_index(input_string,"adios_param")
     read(input_string(name_idx:len_trim(input_string)),nml=adios_param)
#endif
  endif
#endif
  !***************************************************************************


  ! Write out whole parameters
  if(sml_mype==0) then
     open(unit=15,file='fort.input.used',status='replace')
     write(15,nml=sml_param)
     write(15,nml=ptl_param)
     write(15,nml=eq_param)
     write(15,nml=col_param)
     write(15,nml=diag_param)
     write(15,nml=smooth_param)
     write(15,nml=mon_param)
     if(sml_source) write(15,nml=src_param)
     if(sml_radiation) write(15,nml=rad_param)
     if(sml_neutral .or. sml_radiation) write(15,nml=neu_param)
     if(sml_f0_grid) write(15,nml=f0_param)
#ifdef ADIOS  
     write(15,nml=adios_param)
#endif
     close(15)
  endif
  ! 6.  Process  initializing with changing unit


end subroutine read_input

subroutine mem_allocation
  use ptl_module
  use itp_module
  use sml_module
  use eq_module, only : eq_mpsi, eq_mr, eq_mz
  implicit none


  ! itp module allocation
  itp_mpsi=eq_mpsi
  itp_mr=eq_mr
  itp_mz=eq_mz

  allocate(itp_psi_knot(itp_mpsi+itp_korder), itp_I_cscoef(4,itp_mpsi),itp_I_break(itp_mpsi))
  allocate(itp_r_knot(itp_mr+itp_korder_rz),itp_z_knot(itp_mz+itp_korder_rz), &
       itp_psi_bscoef(itp_mr*itp_mz))


end subroutine mem_allocation



subroutine get_mid_r_setup
  use sml_module
  use eq_module
  use itp_module, only : itp_max_r
  implicit none
  real (kind=8) :: r1,r2,r0,rend,dr
  integer :: i,j,find
  integer, parameter :: NN=500
  real(kind=8), external :: psi_interpol
  real (kind=8) :: sqrt_out, sqrt_in, sqrt_p1, sqrt_p2, sqrt_dest

  sqrt_out=sqrt(sml_outpsi)
  sqrt_in=sqrt(sml_inpsi)
  eq_mid_r_dp= (sqrt_out-sqrt_in)/real(eq_mid_r_npsi-1)
  r0=eq_axis_r
  rend=itp_max_r
  dr= (itp_max_r-r0)/real( NN )
  j=1

  do i=1, eq_mid_r_npsi
     sqrt_dest= real(i-1) * eq_mid_r_dp + sqrt_in
     find=0
     ! find proper j ( R-index )
     j=1
     do while( find/=1 )
        r1= real(j-1) * dr + r0
        r2= real(j)* dr + r0
        sqrt_p1=sqrt(psi_interpol(r1,eq_axis_z,0,0))
        sqrt_p2=sqrt(psi_interpol(r2,eq_axis_z,0,0))
!        if( sqrt_p1 < sqrt_dest .and. sqrt_dest <= sqrt_p2 ) then
!           find=1
!        else if ( sqrt_dest <= sqrt_p1 ) then
!           print *, 'Error in get_mid_R_setup'
!           stop
!           find=1
        if( sqrt_dest < sqrt_p2) then
           find=1
        else
           j=j+1
        endif
     enddo
     eq_mid_r_psi(i)= ((sqrt_dest-sqrt_p1)*r2 + (sqrt_p2-sqrt_dest)*r1) / (sqrt_p2-sqrt_p1)
     eq_mid_r_psi(i)=max(eq_mid_r_psi(i),r0)
!     print *, i, eq_mid_r_psi(i), r2, r1, psi1, psi_dest, psi2
  end do

end subroutine get_mid_r_setup

real (kind=8) function get_mid_r(psi)
  use eq_module
  use sml_module
  implicit none
  real  (kind=8),intent(in) :: psi
  real (kind=8) :: aa,bb, sqrt_psi, sqrt_in
  integer :: i
  real , parameter :: epsilon=1e-10

  if( psi < sml_inpsi-epsilon .or. psi > sml_outpsi+epsilon) print * , 'Warning : psi in get_mid_r is out of interpolation range', psi/eq_x_psi,psi

  sqrt_psi=sqrt(psi)
  sqrt_in=sqrt(sml_inpsi)


  i=int((sqrt_psi-sqrt_in)/eq_mid_r_dp)+1
  i=min(eq_mid_r_npsi-1,max(1,i))
  bb= (sqrt_psi-sqrt_in)/eq_mid_r_dp + 1D0 - real(i)
  aa=1D0-bb

  get_mid_r = eq_mid_r_psi(i)*aa + eq_mid_r_psi(i+1) * bb

end function get_mid_r

!not used anymore
!subroutine reset_nphi_total
!  use sml_module
!
!  if(.not. sml_turb_poisson) then
!     sml_nphi_total=1
!  endif
!
!end subroutine reset_nphi_total

subroutine add_dir_name1
  use sml_module, only : sml_input_file_dir
  use eq_module, only : eq_filename

  implicit none
  integer :: endpos

  endpos=len_trim(sml_input_file_dir)
  eq_filename       = sml_input_file_dir(1:endpos) // '/' // eq_filename

!!$  print *, eq_filename

end subroutine add_dir_name1
subroutine add_dir_name2
  use sml_module, only : sml_node_file, sml_ele_file, sml_bfollow_file,&
       sml_input_file_dir, sml_add_pot0_file
  use eq_module
  use neu_module
  use rad_module
  implicit none
  integer :: endpos

  endpos=len_trim(sml_input_file_dir)
  sml_node_file     =  sml_input_file_dir(1:endpos) // '/' // sml_node_file
  sml_ele_file      =  sml_input_file_dir(1:endpos) // '/' // sml_ele_file
  sml_bfollow_file  =  sml_input_file_dir(1:endpos) // '/' // sml_bfollow_file
  sml_add_pot0_file =  sml_input_file_dir(1:endpos) // '/' // sml_add_pot0_file
  eq_den_file       =  sml_input_file_dir(1:endpos) // '/' //eq_den_file
  eq_tempi_file     =  sml_input_file_dir(1:endpos) // '/' //eq_tempi_file
  eq_tempe_file     =  sml_input_file_dir(1:endpos) // '/' //eq_tempe_file
  eq_flowi_file     =  sml_input_file_dir(1:endpos) // '/' //eq_flowi_file
  eq_flowe_file     =  sml_input_file_dir(1:endpos) // '/' //eq_flowe_file
  eq_zeff_file      =  sml_input_file_dir(1:endpos) // '/' //eq_zeff_file

  
  rad_filename     =  sml_input_file_dir(1:endpos) // '/' //rad_filename

  neu_sepfile      =  sml_input_file_dir(1:endpos) // '/' //neu_sepfile
  neu_limfile      =  sml_input_file_dir(1:endpos) // '/' //neu_limfile


!!$  print *, 'node',sml_node_file
!!$  print *, 'ele',sml_ele_file
!!$  print *, 'bf',sml_bfollow_file
!!$  print *, 'neu',neu_sepfile
!!$  print *, 'lim',lim_filename

end subroutine add_dir_name2

subroutine check_point(str)
  use sml_module
  implicit none
  character (len=*) :: str
  integer :: n,ierr
  include 'mpif.h'

  call mpi_barrier(sml_comm,ierr)
  if(sml_mype==0) print *, str
end subroutine check_point

subroutine check_point0(str)
  use sml_module
  implicit none
  character (len=*) :: str
  integer :: n,ierr
  include 'mpif.h'

  if(sml_mype==0) print *, str
end subroutine check_point0

subroutine new_communicator
   use sml_module
   implicit none
   include 'mpif.h'
   integer :: i, j, k, ierr, plane_0_pe
   integer :: mpi_comm_world_group, new_sml_comm_group, sml_comm_group
   integer :: plane_ranks(0:sml_pe_per_plane-1)
   integer :: intpl_ranks(0:(sml_totalpe/sml_pe_per_plane)-1)
#ifdef ADIOS
   integer :: adios_ranks(0:(sml_totalpe/sml_pe_per_plane)-1)
#endif
#ifndef PLANE_MAJOR
   integer :: new_sml_comm_ranks(0:sml_totalpe-1)
#endif

#ifdef PLANE_MAJOR

   ! get the group underlying sml_comm
   call mpi_comm_group(sml_comm, sml_comm_group, ierr)

#else

   call check_point('Interplane-Major ordering for sml_comm communicator')

   ! redefine sml_comm pe ordering from consecutive within planes
   ! (Plane_Major) to consecutive across planes (Interplane-Major)
   k = 0
   do j=0,sml_pe_per_plane-1
      do i=0,(sml_totalpe/sml_pe_per_plane)-1
         new_sml_comm_ranks(k)=sml_pe_per_plane*i + j
         k=k+1
      enddo
   enddo

   ! get the group underlying sml_comm (== mpi_comm_world)
   call mpi_comm_group(MPI_COMM_WORLD, mpi_comm_world_group,ierr)

   ! create new group permuting sml_comm pe ordering to Interplane-Major
   call mpi_group_incl(mpi_comm_world_group, sml_totalpe, &
                       new_sml_comm_ranks, new_sml_comm_group, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_group_incl - new_sml_comm',ierr

   ! Free original sml_comm communicator
   Call mpi_comm_free(sml_comm, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_comm_free - sml_comm',ierr

   ! Create the new communicator
   call mpi_comm_create(MPI_COMM_WORLD, new_sml_comm_group, sml_comm,ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_comm_create - new_sml_comm',ierr

   ! reset sml_mype for the new communicator ordering
   call mpi_comm_rank(sml_comm,sml_mype,ierr)

   ! "get" the group underlying sml_comm
   sml_comm_group = new_sml_comm_group

#endif

   ! PLANE MPI COMMUNICATOR
   call check_point('plane mpi communication')
   plane_0_pe=int(sml_mype/sml_pe_per_plane)*sml_pe_per_plane
   do i=0, sml_pe_per_plane-1
      plane_ranks(i)=plane_0_pe + i
   enddo

   ! Create the new plane group
   call mpi_group_incl(sml_comm_group, sml_pe_per_plane, plane_ranks, &
                       sml_plane_group, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_group_incl - plane',ierr

   ! Create the new plane communicator
   call mpi_comm_create(sml_comm, sml_plane_group, sml_plane_comm, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_comm_create - plane',ierr

   call mpi_comm_size(sml_plane_comm, sml_plane_totalpe, ierr)
   call mpi_comm_rank(sml_plane_comm, sml_plane_mype, ierr)
   sml_plane_index=sml_mype/sml_plane_totalpe

   ! INTER-PLANE MPI COMMUNICATOR
   call check_point('inter-plane mpi communication')
   do i=0, (sml_totalpe/sml_pe_per_plane)-1
      intpl_ranks(i)=sml_plane_mype + i*sml_pe_per_plane
   enddo

   ! Create the new inter-plane group
   call mpi_group_incl(sml_comm_group, (sml_totalpe/sml_pe_per_plane), &
                       intpl_ranks, sml_intpl_group, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_group_incl - intpl',ierr

   ! Create the new inter-plane communicator
   call mpi_comm_create(sml_comm, sml_intpl_group, sml_intpl_comm, ierr)
   if(ierr/=MPI_SUCCESS) print *, 'error in mpi_comm_create - intpl',ierr

   call mpi_comm_size(sml_intpl_comm, sml_intpl_totalpe, ierr)
   call mpi_comm_rank(sml_intpl_comm, sml_intpl_mype, ierr)

#ifdef ADIOS
   call check_point('adios mpi communication?')
   do i=0, sml_intpl_totalpe-1
      adios_ranks(i)= i*sml_pe_per_plane
   enddo
   call mpi_group_incl(sml_comm_group, sml_intpl_totalpe, adios_ranks, sml_adios_group, ierr)
   call mpi_comm_create(sml_comm, sml_adios_group, sml_adios_comm, ierr)

   if (adios_mfilesys) then
      sml_pe_per_filesys = ceiling(real(sml_totalpe)/2.0)
      sml_mfilesys_index = sml_mype/sml_pe_per_filesys
      sml_mfilesys_mype = mod(sml_mype, sml_pe_per_filesys)
      call MPI_Comm_split(sml_comm, sml_mfilesys_index, sml_mfilesys_mype, sml_mfilesys_comm, ierr)
   endif
#endif

   call check_point('end of commnunicator setup')
end subroutine

! new_mpi_ptl_type -- for collision 3 shift - by E. Yoon
! sp list to be sent
! 1.  sendright(iphs:iphe,m)=sp%ptl(ptl_rshift(m))%ph : real, ptl_nphase=6
! 2.  sendright(icts:icte,m)=sp%ptl(ptl_rshift(m))%ct : real, ptl_nconst=3
! 3.  sendright(iph0s:iph0e,m)=sp%phase0(:,ptl_rshift(m)) : real, ptl_nphase
! 4.  sendrightid(m)=sp%ptl(ptl_rshift(m))%gid : integer, 1
! The front number implies packing order
subroutine new_mpi_ptl_type
   use sml_module
   use ptl_module, only : ptl_nphase, ptl_nconst
   use col_module, only : mpi_ptl_type
   implicit none
   include 'mpif.h'
   integer :: blockcounts(0:3), offsets(0:3), oldtypes(0:3), extent
   integer :: ierr, st

   blockcounts(0) = ptl_nphase
   blockcounts(1) = ptl_nconst
   blockcounts(2) = ptl_nphase
   blockcounts(3) = 1

   oldtypes(0) = MPI_REAL8
   oldtypes(1) = MPI_REAL8
   oldtypes(2) = MPI_REAL8
   oldtypes(3) = MPI_INTEGER

   call MPI_TYPE_EXTENT(MPI_REAL8, extent, ierr)
   offsets(0) = 0
   offsets(1) = ptl_nphase*extent
   offsets(2) = offsets(1)+ptl_nconst*extent
   offsets(3) = offsets(2)+ptl_nphase*extent

   ! Define structured type and commit it
   call MPI_TYPE_STRUCT(4,blockcounts, offsets, oldtypes, mpi_ptl_type, ierr)
   call MPI_TYPE_COMMIT(mpi_ptl_type, ierr)
end subroutine


subroutine setup_grid(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module, only : eq_x_psi, eq_x_z, eq_rgrid, eq_zgrid, eq_mr, eq_mz
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: ierr,i, npsi, inode, k
  real (kind=8) :: psi0, epsi
  integer, parameter :: wall_num=100
  real (8), allocatable :: ptmp(:)

  grid%guess_n=(/sml_guess_table_size,sml_guess_table_size/)
  grid%guess_min=(/sml_bd_min_r,sml_bd_min_z/)
  grid%guess_max=(/sml_bd_max_r,sml_bd_max_z/)
  if(sml_mype==0)  then
     print *,'Guess min', grid%guess_min
     print *,'Guess max', grid%guess_max
  endif

  if(sml_mype==0) print *, 'init_grid'
  call init_grid(grid,sml_node_file,sml_ele_file)
  !min-max of grid
  grid%npsi00=sml_00_npsi
  grid%psi00min=sml_inpsi
  grid%psi00max=sml_outpsi
  grid%dpsi00=(grid%psi00max-grid%psi00min)/real(grid%npsi00-1)
  !rho
  grid%nrho=sml_grid_nrho  ! rho indexing starts from 0
  grid%rhomax=sml_rhomax   ! rho min is 0.
  grid%drho = grid%rhomax/real(grid%nrho)

  !rh This is shifted to init_grid
  !rh if(sml_mype==0) print *, 'init_guess_table'
  !rh call init_guess_table(grid)

  call psn_mem_alloc(psn,grid%nnode,grid%ntriangle,grid%npsi00,grid%nrho,sml_nhybrid)
!!$  if(sml_sheath_mode/=0) then
!!$     !find wall point
!!$     do i=1, grid%nnode
!!$        if(grid%p(i)==wall_num) exit
!!$     enddo
!!$     psn%wall_start=i
!!$     psn%nwall=grid%nnode-i+1
!!$     call psn_sheath_alloc(psn)
!!$  endif
#ifdef OLD_INIT_FF
  call check_point('init bfollow, init_ff')
  if(sml_turb_poisson) then
     call init_bfollow(grid,psn)
     call init_ff(grid,psn,.true.)
     call init_ff(grid,psn,.false.)
  endif
#else
  call check_point('calling init_ff')
  if(sml_turb_poisson) then
     call init_ff(grid,psn)
  endif
#endif

  call check_point("after init_ff")
  !call extend_boundary2(grid,psn)  ! moved after init_sheath

  !  call get_matrix00(grid,psn)
#ifdef XGC_DEBUG1
  if(sml_mype==0) call output_matrix(grid,psn)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Identify flux surfaces and calculate safety factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call check_point("Some simple operations 1")

  !1. find the flux surface number: npsi
  epsi = 2D-4
  npsi=1
  psi0 = grid%psi(1)
  do inode =2, grid%nnode
     if(abs(grid%psi(inode)-grid%psi(inode-1))>=epsi) then
       npsi = npsi + 1
     endif
  enddo
  grid%nsurf=npsi

  call check_point("Some simple operations 2")

  !2. set itheta0 and ntheta0 for each flux surface
  allocate(grid%itheta0(npsi), grid%ntheta0(npsi))
  k=1
  grid%itheta0(1) = 1
  grid%ntheta0(1) = 1
  do inode = 2, grid%nnode
     if(abs(grid%psi(inode)-grid%psi(inode-1)).lt.epsi) then
        grid%ntheta0(k) = grid%ntheta0(k)+1
     else
        k = k+1
        grid%itheta0(k) = inode
        grid%ntheta0(k) = 1
     endif
  enddo

  call check_point("q evaluation")

  allocate(ptmp(npsi))  !tmp variable to avoid warning message of ifort -C
  ptmp=grid%psi(grid%itheta0)
  call q_evaluation(grid,eq_rgrid,eq_zgrid,eq_mr,eq_mz,grid%nsurf,ptmp)
  deallocate(ptmp)

  call check_point("After q evaluation")

#ifdef XGC1_EM
  !  SK-HELP B_cross_gradB_cal is in EM
  if(sml_gstep==1 .and. sml_mype==0) then ! for debug
    do i=1, npsi
       write(1234,*) grid%psi(grid%itheta0(i)), grid%qsafety2(i), grid%j0_par(grid%itheta0(i))
    enddo
  endif
#endif

  call check_point("setup end")

end subroutine setup_grid

subroutine extend_boundary2(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  real (kind=8) :: in_bd_psi, out_bd_psi
  logical :: is_in, is_out, is_private, is_near_wall
  real (kind=8), parameter :: epsil=0.005


  !is_near_wall : include nodes which is connected to wall node

  ! nonzero n-mode,  potential boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta1H*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta3H*eq_x_psi
  is_in=.true.
  is_out=.true.
!  is_out=(sml_sheath_mode<=0)
  is_private=.false.
  is_near_wall=.false.
  if(sml_mype==0) print *, 'nonzero n-mode, potential boundary'
  if(sml_mype==0) print *, 'eq_x_psi',eq_x_psi
  if(sml_mype==0) print *, 'in_bd_psi',in_bd_psi
  if(sml_mype==0) print *, 'sml_inpsi',sml_inpsi
  if(sml_mype==0) print *, 'sml_bd_ext_delta1H',sml_bd_ext_delta1H
  if(sml_mype==0) print *, 'out_bd_psi',out_bd_psi
  if(sml_mype==0) print *, 'sml_outpsi',sml_outpsi 
  if(sml_mype==0) print *, 'sml_bd_ext_delta3H',sml_bd_ext_delta3H
  call extend_bd2_single(grid,psn,psn%pbdH_2,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)

  ! nonzero n-mode, charge boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2H*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta4H*eq_x_psi
  is_in=.true.
  is_out=.true.
  is_private=.false.
  is_near_wall=.true.
  if(sml_mype==0) print *, 'nonzero n-mode, charge boundary'
  if(sml_mype==0) print *, 'eq_x_psi',eq_x_psi
  if(sml_mype==0) print *, 'in_bd_psi',in_bd_psi
  if(sml_mype==0) print *, 'sml_inpsi',sml_inpsi
  if(sml_mype==0) print *, 'sml_bd_ext_delta2H',sml_bd_ext_delta2H
  if(sml_mype==0) print *, 'out_bd_psi',out_bd_psi
  if(sml_mype==0) print *, 'sml_outpsi',sml_outpsi 
  if(sml_mype==0) print *, 'sml_bd_ext_delta4H',sml_bd_ext_delta4H
  call extend_bd2_single(grid,psn,psn%cbdH_2,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)

#ifdef XGC1_EM
 ! fluid electron boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2H*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta4H*eq_x_psi
  is_in=.true.
  is_out=.true.
  is_private=.false.
  is_near_wall=.true.
  if(sml_mype==0) print *, 'nonzero n-mode, charge boundary'
  call extend_bd2_single(grid,psn,psn%febdH_2,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)
#endif

  ! n=0-mode, potential boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta1*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta3*eq_x_psi
  if(sml_rgn1_pot0_only) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta3), out_bd_psi)
  is_in=(sml_zero_inner_bd==1)
  is_out=.true. !(sml_sheath_mode<=0)
  is_private=.NOT. sml_exclude_private
  is_near_wall=.false.
  if(sml_mype==0) print *, 'n=0-mode, potential boundary'
  if(sml_mype==0) print *, 'eq_x_psi',eq_x_psi
  if(sml_mype==0) print *, 'in_bd_psi',in_bd_psi
  if(sml_mype==0) print *, 'sml_inpsi',sml_inpsi
  if(sml_mype==0) print *, 'sml_bd_ext_delta1',sml_bd_ext_delta1
  if(sml_mype==0) print *, 'out_bd_psi',out_bd_psi
  if(sml_mype==0) print *, 'sml_outpsi',sml_outpsi 
  if(sml_mype==0) print *, 'sml_bd_ext_delta3',sml_bd_ext_delta3
  call extend_bd2_single(grid,psn,psn%pbd0_2,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)

  ! n=0-mode, charge boundary
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2*eq_x_psi
  out_bd_psi=sml_outpsi - sml_bd_ext_delta4*eq_x_psi
  if(sml_rgn1_pot0_only) out_bd_psi=min(eq_x_psi*(1D0-sml_bd_ext_delta4),out_bd_psi)
  is_in=.true.
  is_out=.true.
  is_private=.NOT. sml_exclude_private
  is_near_wall=.true.
  if(sml_mype==0) print *, 'n=0-mode, charge boundary'
  if(sml_mype==0) print *, 'eq_x_psi',eq_x_psi
  if(sml_mype==0) print *, 'in_bd_psi',in_bd_psi
  if(sml_mype==0) print *, 'sml_inpsi',sml_inpsi
  if(sml_mype==0) print *, 'sml_bd_ext_delta2',sml_bd_ext_delta2
  if(sml_mype==0) print *, 'out_bd_psi',out_bd_psi
  if(sml_mype==0) print *, 'sml_outpsi',sml_outpsi 
  if(sml_mype==0) print *, 'sml_bd_ext_delta4',sml_bd_ext_delta4
  call extend_bd2_single(grid,psn,psn%cbd0_2,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)

  ! additional boundary to remove axis density 
  in_bd_psi=sml_inpsi + sml_bd_ext_delta2H*eq_x_psi
  call extend_bd2_single(grid,psn,psn%cbd0_tmp,in_bd_psi,out_bd_psi,is_in, is_out, is_private, is_near_wall)
end subroutine


subroutine extend_bd2_single(grid,psn,bd,in_bd_psi,out_bd_psi,is_in,is_out,is_private, is_near_wall)
  use grid_class
  use psn_class
  use sml_module
  use boundary_class
  use eq_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(boundary2_type) :: bd
  real (kind=8) :: in_bd_psi, out_bd_psi
  logical, intent(in) :: is_in, is_out, is_private, is_near_wall
  !
  integer :: i
  real (kind=8) , parameter :: epsilon=1D-8
  logical :: lbd
  integer :: start, end, nseg

  ! count segments
  nseg=0
  start=0

  do i=1, grid%nnode
     call is_bd(i,is_near_wall,lbd)
     if(start==0) then
        if(lbd) then
           nseg=nseg+1
           start=i
           end=i
        endif
     else ! start/=0 --> previous node (i-1) is boundary node
        if(lbd) then
           end=end+1
        else
           start=0
        endif
     endif
  enddo


  ! allocated memory
  bd%nseg=nseg
  allocate(bd%iseg(nseg))

  nseg=0
  start=0

  do i=1, grid%nnode
     call is_bd(i,is_near_wall,lbd)
     if(start==0) then
        if(lbd) then
           nseg=nseg+1
           start=i
           end=i
        endif
     else ! start/=0 --> previous node (i-1) is boundary node
        if(lbd) then
           end=end+1
        else
           start=0
        endif
     endif
     if(start/=0) then
        bd%iseg(nseg)%start=start
        bd%iseg(nseg)%end=end
     endif
  enddo

  if(sml_mype==0) then
     print *, 'extended bd nseg', bd%nseg
     do i=1, bd%nseg
        print *, '     iseg(',i,')=' , bd%iseg(i)
     enddo
  endif

contains
  subroutine is_bd(i,is_near_wall,is_bd1)
    implicit none
    integer :: i
    logical :: is_near_wall, is_bd1
    ! when is_near_wall is on
    integer :: j,tr, k, nd
    real (8) :: r1, z1, r2, z2, dist, min_dist

    if(is_in .and. grid%psi(i) - epsilon < in_bd_psi) then
       is_bd1=.true.
       return
    elseif(is_out .and. grid%psi(i)> out_bd_psi) then
       is_bd1=.true.
       return
    elseif(.not. is_private .and. grid%rgn(i)/=1 .and. grid%rgn(i)/=2 ) then
       is_bd1=.true.
       return
    else if(is_near_wall) then
       !check if i-node is connected to wall
       do j=1, grid%num_t_node(i)
          tr=grid%tr_node(j,i)
          do k=1,3
             nd=grid%nd(k,tr)
             if(grid%rgn(nd)==grid_rgn_wall) then
                is_bd1=.true.
                return
             endif
          enddo
       enddo
       ! check if i-node is near wall (sml_bd_ext_near_wall)
       if(sml_bd_ext_near_wall>0D0) then
          min_dist=1D99
          r1=grid%x(1,i)
          z1=grid%x(2,i)
          do j=1, psn%nwall
            r2=grid%x(1,psn%wall_nodes(j))
            z2=grid%x(2,psn%wall_nodes(j))
            dist=(r1-r2)*(r1-r2)+(z1-z2)*(z1-z2)
            min_dist=min(min_dist,dist)
            if(min_dist < sml_bd_ext_near_wall*sml_bd_ext_near_wall) then
              is_bd1=.true.
              return
            endif
          enddo
       endif


       ! wall node not found
       is_bd1=.false.
    else
       is_bd1=.false.
    endif

  end subroutine is_bd
end subroutine extend_bd2_single


 !! read grid information from file
subroutine init_grid(grid,filename1,filename2)
  use grid_class
  use sml_module
  use eq_module
  use fld_module
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  character (len=65), intent(in) :: filename1,filename2 !! (1) node file, (2) triangle element file
!  real (8), intent(in) :: xpsi, xz
  integer, parameter :: nodefile=2001,trifile=2002
  integer :: i,j,k,dum,n_n,n_t,ierr
  real (kind=8) , external :: psi_interpol
  real (kind=8) :: psi_val
#ifdef XGC1_EM
  real (kind=8), allocatable :: j0_par_eb(:), grad_j0_par_eb(:,:)
#endif
  ! For magnetic field info
  type(fld_type) :: fld
  real (kind=8) :: over_B, over_B2, inv_r, dbdr, dbdz, dbdphi, grad_psi(2), over_abs_grad_psi
  logical :: rz_outside


  if(sml_mype==0) then
     !read node -----------------------------------------------------------------
     !open file
     print *, 'open node file (*.node) :', filename1
     open(nodefile, file=filename1,status='old',form='formatted')
     print *, 'open connectivity file (*.ele) :', filename2
     open(trifile, file=filename2,status='old',form='formatted')

     !read
     ! read number of nodes
     read(nodefile,*) grid%nnode, dum,dum,dum
  endif

  call MPI_BCAST(grid%nnode, 1,  MPI_INTEGER, 0, sml_comm,ierr)

  n_n=grid%nnode
  ! read nodes
  allocate(grid%x(2,n_n),grid%rgn(n_n),grid%psi(n_n))
  allocate(grid%node_area(n_n),grid%inv_node_vol(n_n),grid%node_vol(n_n), grid%node_vol_ff(n_n,0:1))
  if(sml_f0_grid) allocate(grid%node_vol_nearest(n_n))
  allocate(grid%bfield(4,n_n))
  allocate(grid%v_curv(3,n_n),grid%v_gradb(3,n_n),grid%nb_curl_nb(n_n))

#ifdef XGC1_EM
  allocate(grid%gradpsi(n_n,2),grid%absgradpsi(n_n))
  allocate(grid%dene(n_n), grid%tempe(n_n))
  allocate(grid%B_cross_gradB(3,n_n), grid%j0_par(n_n),grid%curl_nb(2,n_n))
  allocate(grid%tearing_drive(n_n,3),grid%tearing_drive2(n_n))
#endif

  !    allocate(grid%sort(n_n),grid%inv_sort(n_n))
  if(sml_mype==0) then
     do i=1, n_n
        !          read(nodefile,*) dum, grid%x(1,i),grid%x(2,i),dum, dum, dum, &
        !               dum, dum, dum, dum, dum
        read(nodefile,*) dum, grid%x(1,i),grid%x(2,i),grid%rgn(i)
        !    grid%nn(2,i),grid%rgn(i),grid%it(i),old_node_num(i),bd_flag(i)
     enddo
  endif
  !Broad Casting data
  call MPI_BCAST(grid%x, n_n*2, MPI_DOUBLE_PRECISION, 0, sml_comm,ierr)
  ! end of reading and broadcasting -------------------------------------------------


  if(sml_mype==0) print *, 'grid%nnode =', grid%nnode

  allocate(grid%rtmp1(n_n),grid%rtmp2(n_n))

  !set psi and rgn
  if(sml_mype==0) then
     do i=1, n_n
        grid%psi(i)=psi_interpol(grid%x(1,i),grid%x(2,i),0,0)
        if(grid%rgn(i)==1) then ! rgn was boundary flag
           grid%rgn(i)=grid_rgn_wall  ! region number for wall
        else
           if(grid%psi(i) > eq_x_psi -epsil_psi ) then
              grid%rgn(i)=2
           elseif( is_rgn1(grid%x(1,i),grid%x(2,i),grid%psi(i)) ) then
              grid%rgn(i)=1
           else
              grid%rgn(i)=3
           endif
        endif
     enddo
  endif
  call MPI_BCAST(grid%rgn, n_n, MPI_INTEGER, 0, sml_comm,ierr)
  call MPI_BCAST(grid%psi, n_n, MPI_REAL8,   0, sml_comm,ierr)

  ! read number of triangles - element file
  if(sml_mype==0) then
     read(trifile,*) grid%ntriangle,dum,dum
  endif
  !Broad Casting data
  call MPI_BCAST(grid%ntriangle, 1,  MPI_INTEGER, 0, sml_comm,ierr)

  n_t=grid%ntriangle
  ! number of triangle
  !    --------------------------------------------------------
  !    use the column grid%mapping(1:2,3,itr) in mapping array to hold the
  !    coordinate for 3rd vertex in triangle itr
  !    --------------------------------------------------------
  allocate( grid%nd(3,n_t), grid%adj(3,n_t),grid%mapping(2,3,n_t))
  allocate( grid%tr_vol(n_t),grid%tr_area(n_t) )
  !read from file
  if(sml_mype==0) then
     do i=1, grid%ntriangle
        read(trifile,*) dum,grid%nd(1,i),grid%nd(2,i),grid%nd(3,i)
     enddo
  endif

  call MPI_BCAST(grid%nd, n_t*3, MPI_INTEGER, 0, sml_comm,ierr)
  ! end of reading and broadcasting ----------------------------------------------

  call init_triangle(grid)

  !rh We have to call this before init_gradient_mat, because we already need the triangle search
  if(sml_mype==0) print *, 'init_guess_table'
  call init_guess_table(grid)

  call get_nphi(grid%nphi)
  grid%delta_phi=atan(1.D0)*8D0/real(sml_nphi_2pi)
  grid%iphi_offset=grid%nphi*sml_intpl_mype
  grid%phimin=real(grid%iphi_offset)*grid%delta_phi
  grid%phimax=grid%phimin + grid%delta_phi*real(grid%nphi)
  !    print *, mype, grid%iphi_offset, grid%phimin/atan(1D0), grid%phimax/atan(1D0)

  do i=1, n_n
     call bvec_interpol(grid%x(1,i),grid%x(2,i),0D0,grid%bfield(1,i),grid%bfield(2,i),grid%bfield(3,i))
     grid%bfield(4,i)=sqrt(sum(grid%bfield(1:3,i)**2))

     ! Prepare data for evaluation of v-grid drift-velocities
     fld%r=grid%x(1,i)
     fld%z=grid%x(2,i)
     call field(fld,0.D0,rz_outside)
     over_B=1.D0/grid%bfield(4,i)
     over_B2=over_B*over_B
     if (sml_cylindrical) then
       inv_r=1.D0/eq_axis_r
     else
       inv_r=1.D0/fld%r
     endif

     grad_psi(1)=psi_interpol(fld%r,fld%z,1,0)
     grad_psi(2)=psi_interpol(fld%r,fld%z,0,1)
     over_abs_grad_psi=1.D0/sqrt(sum(grad_psi(:)**2))  ! --> What to do with X-point???

     if (sml_cylindrical) then
       ! Cylindrical limit
       grid%nb_curl_nb(i)= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - fld%dbpdr*fld%bz )
     else
       ! Toroidal geometry
       grid%nb_curl_nb(i)= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/fld%r  + fld%dbpdr)*fld%bz )
     endif

     dbdr   = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
     dbdz   = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
     dbdphi = 0.D0

     ! To get the real velocities, one has to multiply D and some mu and rho factors
     ! v_gradb = Bxgrad(|B|)
     ! v_curv = curl(B)
     if (sml_cylindrical) then
       ! Radial curv drifts
       grid%v_curv(1,i)  = grad_psi(1) * (fld%dbzdp*inv_r - fld%dbpdz) &
                          +grad_psi(2) * (fld%dbpdr-fld%dbrdp*inv_r)
       ! Poloidal curv drift
       grid%v_curv(2,i)  = (-grad_psi(2) * (fld%dbzdp*inv_r - fld%dbpdz) &
                          +grad_psi(1) * (fld%dbpdr-fld%dbrdp*inv_r) ) &
                         * over_abs_grad_psi
     else
       ! Radial curv drifts
       grid%v_curv(1,i)  = grad_psi(1) * (fld%dbzdp*inv_r - fld%dbpdz) &
                          +grad_psi(2) * (fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r)
       ! Poloidal curv drift
       grid%v_curv(2,i)  = (-grad_psi(2) * (fld%dbzdp*inv_r - fld%dbpdz) &
                            +grad_psi(1) * (fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r) ) &
                           * over_abs_grad_psi
     endif
     ! Radial grad(B) drifts
     grid%v_gradb(1,i) = grad_psi(1) * (fld%bphi * dbdz - fld%bz * dbdphi/fld%r) &
                        +grad_psi(2) * (fld%br * dbdphi/fld%r - fld%bphi * dbdr)
     ! Poloidal grad(B) drifts
     grid%v_gradb(2,i) = (-grad_psi(2) * (fld%bphi * dbdz - fld%bz * dbdphi/fld%r) &
                          +grad_psi(1) * (fld%br * dbdphi/fld%r - fld%bphi * dbdr) ) &
                         * over_abs_grad_psi
     ! Toroidal magnetic drifts
     grid%v_curv(3,i)  = fld%dbrdz - fld%dbzdr
     grid%v_gradb(3,i) = fld%bz * dbdr - fld%br * dbdz

#ifdef XGC1_EM
     grid%gradpsi(i,1)=grad_psi(1)
     grid%gradpsi(i,2)=grad_psi(2)
     grid%absgradpsi(i)=sqrt(sum(grad_psi**2))

     grid%B_cross_gradB(1,i) = fld%bphi*dbdz         ! R
     grid%B_cross_gradB(2,i) = -fld%bphi*dbdr           ! Z
     grid%B_cross_gradB(3,i) = -(fld%br*dbdz-fld%bz*dbdr) ! phi
  
     grid%j0_par(i) = (fld%dbrdz - fld%dbzdr)/sml_mu0
     if (sml_cylindrical) then
       grid%curl_nb(1,i) = fld%bphi*dbdz*over_B2  - fld%dbpdz*over_B !R
       grid%curl_nb(2,i) = -fld%bphi*dbdr*over_B2 + fld%dbpdr*over_B !Z
     else
       grid%curl_nb(1,i) = fld%bphi*dbdz*over_B2  - fld%dbpdz*over_B !R
       grid%curl_nb(2,i) = -fld%bphi*dbdr*over_B2 + (fld%bphi*inv_r + fld%dbpdr)*over_B !Z
     endif

     grid%dene(i)=eq_ftn(grid%psi(i),fld%r,fld%z,eq_den)
     grid%tempe(i)=eq_ftn(grid%psi(i),fld%r,fld%z,eq_tempi)*sml_e_charge ! equillibrium electron temperature
#endif

  enddo

  call get_node_vol(grid)

  !rh FD gradient operator in (R,Z) plane
  call init_gradient_mat(grid)
  !rh call init_pol_smooth_mat(grid)

  do i=1, n_n
     call bvec_interpol(grid%x(1,i),grid%x(2,i),0D0,grid%bfield(1,i),grid%bfield(2,i),grid%bfield(3,i))
     grid%bfield(4,i)=sqrt(sum(grid%bfield(:,i)**2))
  enddo

#ifdef XGC1_EM
     !rh Calculate the tearing drive term b x grad(j0/(eB))
     allocate(grad_j0_par_eb(grid%nnode,2),j0_par_eb(grid%nnode))
     ! -j0/(eB) ---> We use "-" because it is an electron current
     ! and the kink drive contains grad(n0*u_e0/B)!
     ! Chen and Parker define e<0, we define e>0!
     do i=1,grid%nnode
       j0_par_eb(i)=grid%j0_par(i)/(-sml_e_charge*grid%bfield(4,i))
     enddo
     call grid_deriv(grid,j0_par_eb,grad_j0_par_eb(:,1),grad_j0_par_eb(:,2))

     do i=1,grid%nnode
       grid%tearing_drive(i,1) = grid%bfield(3,i)*grad_j0_par_eb(i,2)/grid%bfield(4,i)
       grid%tearing_drive(i,2) =-grid%bfield(3,i)*grad_j0_par_eb(i,1)/grid%bfield(4,i)
       grid%tearing_drive(i,3) = ( grid%bfield(2,i)*grad_j0_epar_b(i,1) - grid%bfield(1,i)*grad_j0_par_eb(i,2) ) &
                                / grid%bfield(4,i)
       grid%tearing_drive2(i) = ( grid%curl_nb(1,i)*grad_j0_par_eb(i,1) &
                                 +grid%curl_nb(2,i)*grad_j0_par_eb(i,2) )
     enddo

     ! Generate matrix for (b x grad(j0/eB)).grad
     call make_v_dot_grad_mat(grid,grid%kink_mat,grid%tearing_drive)

     deallocate(grad_j0_par_eb,j0_par_eb)
#endif

end subroutine init_grid

subroutine init_sheath(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  !
!  integer, parameter :: rgn_wall=100
  integer :: i, count

  ! count nwall
  count=0
  do i=1, grid%nnode
     if(grid%rgn(i)==grid_rgn_wall) count=count+1
  enddo
  psn%nwall=count

  call psn_sheath_alloc(psn) ! allocate psn%sheath_pot and psn%wall_nodes
  !if(sml_sheath_adjust) then
  allocate(psn%node_to_wall(grid%nnode)) ! merge to psn_sheath_alloc later
  !endif

  psn%node_to_wall=0
  count=0
  do i=1, grid%nnode

     if(grid%rgn(i)==grid_rgn_wall) then
        count=count+1
        psn%wall_nodes(count)=i

        !if(sml_sheath_adjust) psn%node_to_wall(i)=count
        psn%node_to_wall(i)=count
     endif
  enddo


end subroutine init_sheath
