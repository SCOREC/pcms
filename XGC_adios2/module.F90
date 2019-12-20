#ifndef USE_INQUIRE_DIRECTORY
#define DIR_INQ FILE
#else
#define DIR_INQ DIRECTORY
#endif
!! Simulation parameter and control module
module sml_module
  integer :: sml_comm !! communicator
  integer :: sml_comm_null !! MPI_COMM_NULL communicator for adios posix method
  integer :: sml_mype !! processor index
  integer :: sml_totalpe !! total number of processors
  integer :: sml_nthreads !! total number of OpenMP threads per process
  integer :: sml_mstep  !! total number of time steps for simulation
  integer :: sml_istep  !! current time step number
  integer :: sml_gstep  !! global time step
  integer :: sml_ipc    !! 2nd order Runge-Kutta index, or predictor-corrector index
  integer :: sml_epc    !! 2nd order Runge-Kutta index, or predictor-corrector index for electrons
  real (8) :: sml_time !! simulation time
  real (8) :: sml_dt   !! time step size
  real (8) :: sml_krook_rate !Krook rate term

  logical :: sml_deltaf, sml_deltaf_elec  !! delta-f switch. 0 for off, 1 for on. delta-f simulation is not verified yet, especially for output file.
  logical :: sml_extra_dwdt  ! extra dwdt

  logical :: sml_dwdt_exb_only !! delta-f weight evolution switch. false for whole (grad_b, curl_b,and ExB), true for exb only
  logical :: sml_dwdt_fix_bg  !! delta-f weight evolution switch. false for (1-w) factor true for (1) factor. default is .false.

  real (8) :: sml_marker_temp_factor !! Apply virtual temperature for initial loading. The virtual temperature will be (temp_factor)*(real temperature).
  real (8) :: sml_marker_min_temp ! minimum temperature of marker
  real (8) :: sml_marker_temp_factor2 ! perp temperature factor
  real (8) :: sml_marker_temp_factor3 ! low mu fill temperature
  real (8) :: sml_low_mu_fill_population  !

  !flat marker distribution
  logical :: sml_flat_marker  ! use flat marker distribution function when it is true
  real (8) :: sml_flat_marker_decay_start1, sml_flat_marker_cutoff1, sml_flat_marker_width1
  real (8) :: sml_flat_marker_decay_start2, sml_flat_marker_cutoff2, sml_flat_marker_width2

  logical :: sml_initial_flow ! initial parallel flow applied.

  logical :: sml_concentric !! Simulation for concentric circular geometry
  logical :: sml_read_eq       !! Read B-field profile from a file
  !   real (8), parameter :: sml_2pi=6.28318530717959D0  !! 2 Pi
  !  real (8), parameter :: sml_pi =3.14159265358979D0 !! Pi
  real (8) :: sml_2pi !! 2 Pi
  real (8) :: sml_pi  !!  Pi
  real (8) :: sml_sqrtpi  !! Sqrt(Pi)
  real (8) :: sml_sqrt2pi !! Sqrt(2Pi)
  real (8) :: sml_2pi_wedge_n

  real (8), parameter :: sml_e_charge=1.6022D-19  !! electron charge (MKS)
  real (8), parameter :: sml_epsilon0=8.8542D-12  !! permittivity of free space (MKS)
  real (8), parameter :: sml_prot_mass=1.6720D-27 !! proton mass (MKS)
  real (8), parameter :: sml_elec_mass=9.1094D-31 !! electron mass (MKS)
!  real (8), parameter :: sml_elec_mass=0.01D0*sml_prot_mass !! electron mass (MKS)

  real (8), parameter :: sml_ev2j=sml_e_charge, sml_j2ev=1D0/sml_e_charge ! e, 1/e

#ifdef XGC1_EM
  integer :: sml_nrk
#else
#ifdef PURE_RK4
  integer, parameter :: sml_nrk = 4 ! RK4 step
#else
  integer, parameter :: sml_nrk = 2 ! RK2 step
#endif
#endif

  integer :: sml_bounce         !! Bounce routine switch  0 for off, 1 for inner boundary, 2 for both boundaries
  logical :: sml_restart        !! Restarting simulation switch.  0 for original run, 1 for restart
  integer :: sml_run_count      !! couting number of runs : 1 for fresh start and 2, 3, 4, ... for restart
  logical :: sml_invert_B         !! Change the direction of whole B-field
  logical :: sml_co_cur_bt        !! Change the direction of toroidal B-field only
  real (8) :: sml_bp_sign       !! -1 for inverted B   1 for normal B <-- given by sml_invert_B
  real (8) :: sml_bt_sign       !! -bp_sign for co-current bt, bp_sign for couter-current bt <- given by sml_co_curr_bt

  integer :: sml_monte_num      !! Number of sample for Monte-Carlo volume calculation. It should be greater than particle number for correct simulation.
  real (8) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z !! simulation boundary in R-Z space
  real (8) :: sml_marker_den !! Marker density. Used for loading and initial weight calculation.
  real (8) :: sml_inpsi, sml_outpsi  !! Inner and outer boundary for initial loading.
  integer :: sml_restart_write_period !! Number of steps between particle data dumps for restart.

  real (8) :: sml_en_order_kev  !! Characteristic ion energy (keV)
  real (8) :: sml_en_order      !! Characteristic ion energy in normalized unit
  real (8) :: sml_tran !! Torodial transit time of an ion with the characteristic energy.

  integer :: sml_machine  !! Machine type
  integer :: sml_limiter  !! Switch for limiter
  integer :: sml_special  !! Switch for special simulation. eg. single particle simulation.

  real (8) :: sml_bp_mult !! Artificial multiplication of poloidal B-field
  real (8) :: sml_bt_mult !! Artificial multiplication of toroidal B-field

  !for gyrokinetic
  integer, parameter  :: sml_nlarmor=8 !4
  integer :: sml_nphi_total, sml_nphi_2pi
  integer :: sml_wedge_n

  logical :: sml_turb_efield    ! set zero turbulece efield
  logical :: sml_turb_poisson   ! don't solve poisson equation
  logical :: sml_turb_poisson_n0_only ! this will set dpot=0 right after turb possion solver
  logical :: sml_00_efield      ! set zero 00 efield
  logical :: sml_00_poisson     ! don't solve 00 poisson equation
  real (8) :: sml_gradpsi_limit

  integer :: sml_grid_nrho
  real (8) :: sml_rhomax

  !for multi speces
!  real (8) ::   sml_c2_2m_rel(2), sml_c_m_rel(2),sml_mass_rel(2), sml_charge_rel(2)

  !for ExB suppress
!  integer :: sml_exb_suppress
  real (8) :: sml_exb_on(2)=(/1D0, 0D0/), sml_exb_suppress_time, sml_exb_suppress_width


  ! for electron hybrid scheme
  integer :: sml_nhybrid
  integer :: sml_ncycle_half


  ! for initial distribution relaxation
  integer :: sml_relax_dist
  integer :: sml_relax_dist_num

  ! for delta-f simulation f0 function
!  real (8) :: sml_f0den_edge, sml_f0den_out, sml_f0den_ped_c,sml_f0den_ped_width
  logical :: sml_f0_grid
  integer :: sml_f0_grid_lbal_period
  real (8) :: sml_f0_grid_max_ptl_imbal
  real (8) :: sml_f0_grid_min_ptl_imbal
  real (8) :: sml_f0_grid_init_ptl_imbal
  integer :: sml_f0_nmu_decomp
  real (8) :: sml_f0_grid_alpha
  integer :: sml_f0_grid_alpha_start
  integer :: sml_deltaf_f0_mode
  real (8), allocatable :: sml_f0_n(:), sml_f0_gradn(:,:)
  real (8) :: sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e  ! for sml_deltaf_f0_mode==-1  1/Ln , 1/Lt - artificial f0
  logical :: sml_symmetric_f0g
  logical :: sml_symmetric_f
  integer :: sml_f_source_period
  real (8) :: sml_inv_f_source_period
  logical :: sml_no_fp_in_f  ! f used for f_source will not have particle information --> wrong f but smooth


  ! grid file input
  character (len=65) :: sml_node_file, sml_ele_file, sml_bfollow_file
  character (len=65) :: sml_diag_node_file, sml_diag_ele_file

  ! time averaging potential
  real (8) :: sml_tavg_factor
  integer :: sml_tavg_on

  !multi PE for one plane
  integer :: sml_plane_per_pe, sml_pe_per_plane
#ifdef ADIOS
  integer :: sml_adios_group, sml_adios_comm
  ! set by adios_param namelist
  logical :: adios_stage_restart = .FALSE.
  logical :: adios_stage_restartf0 = .FALSE.
  logical :: adios_stage_f0 = .FALSE.
  logical :: adios_stage_field3d = .FALSE.
  logical :: adios_stage_particle = .FALSE.
  logical :: adios_stage_dpot = .FALSE.
  logical :: adios_stage_streamer = .FALSE.
  logical :: adios_mfilesys = .FALSE.
  ! parameters for using multiple filesystem
  integer :: sml_pe_per_filesys, sml_mfilesys_index, sml_mfilesys_mype, sml_mfilesys_comm
#endif
  integer :: sml_plane_group, sml_plane_comm, sml_plane_totalpe, &
             sml_plane_mype, sml_plane_index
  integer :: sml_intpl_group, sml_intpl_comm, sml_intpl_totalpe, &
             sml_intpl_mype
  ! electron on/off switch
  logical :: sml_electron_on

  ! fluid hybrid electron
  logical :: sml_electron_hyb
  logical :: sml_use_ts_solver
  real (8) :: sml_eta  !! resistivity setup for fluid electrons tearing mode
  logical :: sml_old_grad_perp  ! true: FE grad_perp, false: FD grad_perp
  !rh logical :: sml_mass_solve  ! true: mass matrix solver, false: 1/node_area
  logical :: sml_hyb_tearing_test  ! Only terms for linear tearing mode
  logical :: sml_hyb_alfven_test   ! Only terms for simple Alfven dynamics
  logical :: sml_hyb_linear        ! Only linear terms in electron fluid equations
  logical :: sml_hyb_ion_on        ! False: no ion terms in el. fluid eq, true: ion terms activated
  logical :: sml_lumped_mass       ! Use lumped mass or acc. mass matrix (only with sml_mass_solve=true)
  logical :: sml_use_scaling       ! Use scaling factors in PETSc implicit code
  integer :: sml_hyb_debug         ! Debug level for electron equations
  logical :: sml_hyb_bgraddpe_on   ! Switch on/off b.grad(dpe) term

  integer :: sml_fem_matrix, sml_smtrc_ad_elec

  !for p-c.
  real (8), allocatable :: sml_pc_coeff(:,:)

  ! grid guess table size
  integer :: sml_guess_table_size

  ! delta-f loading parameter
  real (8) :: sml_initial_deltaf_noise

  ! zero inner boundary of 00 solver
  integer :: sml_zero_inner_bd

  ! toroidal mode select for linear GK simulation
  integer :: sml_mode_select_on, sml_mode_select_n
  integer :: sml_mode_initial_n

  !
  real (8) :: sml_bd_ext_delta1, sml_bd_ext_delta2, sml_bd_ext_delta3, sml_bd_ext_delta4
  real (8) :: sml_bd_ext_delta1H, sml_bd_ext_delta2H, sml_bd_ext_delta3H, sml_bd_ext_delta4H
  real (8) :: sml_bd_ext_near_wall

  integer :: sml_bfollow
  ! input file directory
  character (len=65) :: sml_input_file_dir

  integer :: sml_bounce_zero_weight

  logical :: sml_use_pade=.true.
  logical :: sml_use_simple00=.false.
  integer :: sml_simple00_nsmth=0
  logical :: sml_iter_solver
  integer :: sml_iter_solver_niter
  real (8) :: sml_iter_solver_accel
  logical :: sml_fsa_solver
  integer :: sml_poisson_solver_type

  integer :: sml_add_pot0  ! additional potential 0 for no pot, 1 for read-in, 2 for neoclassical pot0 (simple)
  logical :: sml_replace_pot0 ! replace add_pot0 to pot0, which means initial 00 mode potential doesn't change
  character (len=256) :: sml_add_pot0_file
  real (8) :: sml_dpot_min_factor, sml_dpot_max_factor
  integer :: sml_flat_electron_density ! smoothing of background electron density for ion only simulation.   0 for no smoothing (default) 1 for smoothing.
  logical :: sml_suppress_weight_growth  ! if this value is true then weight of deltaf will be between -sml_weight_max and sml_weight_max. default value is .false.
  real (8) :: sml_weight_max   ! maximum weight(= phase(6)) for delta-f simulation. default value is 10.

  ! for restore temperature gradient
  logical :: sml_restore_temp  ! restore temperature gradient. Work for delta-f=1, sml_deltaf_f0_mode=-1, and rk2 method only
  integer :: sml_restore_temp_period ! how frequently restore temperature
  real (8) :: sml_restore_temp_avg_scale ! time average scale - dt base, default=0.01 (means 100 dt scale)
  integer :: sml_max_mat_width

  ! for cold electron at the boundary
  integer :: sml_bd_Te_mode    ! 0 for no effect, 1 for innner only, 2 for both boundary, 3 for outer only
  real (8) :: sml_bd_Te_width !radial decay length of electron temperature at the boundary

  logical :: sml_zero_out_total_charge ! total charge is zero always

  real (8),parameter :: sml_boundary_diagonal=1D-6  ! diagonal value of boundary node point

  logical :: sml_exclude_private

  integer :: sml_sheath_mode
  real (8) :: sml_sheath_init_pot_factor
  logical :: sml_sheath_adjust
  real (8) :: sml_sheath_adjust_factor

  logical :: sml_rgn1_pot0_only

  logical :: sml_source

  logical :: sml_radiation

  logical :: sml_neutral
  integer :: sml_neutral_start_step
  logical :: sml_neutral_use_ion_loss

  logical :: sml_pol_decomp
  logical :: sml_pol_decomp_simple
  real (8):: sml_max_imbalance
  logical :: sml_rpl_on

  logical :: sml_00_xz_up ! use upper regionof xz when obtaining 00 mode
  !for diagnosis
  integer, parameter :: sml_n_vf_diag=8

  integer :: sml_00_npsi

  real (8), parameter :: sml_load_maxe = 12D0
  real (8) :: sml_gpu_ratio
  integer :: sml_perm_gpu_freq
  integer :: sml_sort_gpu_freq

  logical :: sml_ignore_drift_near_wall
  real (8):: sml_ignore_drift_r0
  real (8):: sml_ignore_drift_z0
  real (8):: sml_ignore_drift_slope1
  real (8):: sml_ignore_drift_slope2

  ! electron and ion shift_ie communication and parallel algorithm options
  integer :: sml_elec_shift_max_nthreads
  logical :: sml_elec_use_alltoall
  logical :: sml_elec_use_hs_barrier0
  integer :: sml_elec_large_limit
  logical :: sml_elec_handshake
  logical :: sml_elec_use_hs_barrier1
  logical :: sml_elec_use_sendrecv
  logical :: sml_elec_all_sendrecv
  logical :: sml_elec_use_isend
  logical :: sml_elec_use_rsend

  integer :: sml_ion_shift_max_nthreads
  logical :: sml_ion_use_alltoall
  logical :: sml_ion_use_hs_barrier0
  integer :: sml_ion_large_limit
  logical :: sml_ion_handshake
  logical :: sml_ion_use_hs_barrier1
  logical :: sml_ion_use_sendrecv
  logical :: sml_ion_all_sendrecv
  logical :: sml_ion_use_isend
  logical :: sml_ion_use_rsend
#ifdef XGC1_EM
  real (8), parameter :: sml_mu0=1.2566370614D-6 ! vacuum permeability, magnetic constant
  real (8) :: sml_n1_diff_coef_perp !< artificial perpendicular diffusion coefficient
  real (8) :: sml_n1_diff_coef_para !< artificial parallel diffusion coefficient
#endif
#ifdef CYLINDRICAL
  logical, parameter :: sml_cylindrical=.true.
#else
  logical, parameter :: sml_cylindrical=.false.
#endif
end module sml_module

!! remainining time estimation module
module rem_module
  logical  :: rem_estimate  = .false. !! enable (disable) remaining execution time logic
  real (8) :: rem_starttime = -1.0    !! first recorded wallclock timestamp

  integer  :: rem_walltime          !! estimated time remaining at beginning of job (seconds)
  integer  :: rem_walltime_src = -1 !! where estimated time came from
                                    !! -1: not yet determined; 0: no estimate;
                                    !! 1: from input; 2: from pre_aprun;
                                    !! 3: from xgc_syslog (nonconservative estimate)
                                    !! 4: from xgc_syslog (conservative estimate)

  integer  :: rem_final_step_est      !! estimate of last step can compute before time expires
  logical  :: rem_final_step_upd = .false. !! flag indicating that rem_final_step_est was broadcast

  integer  :: rem_step_loop_cnt = -1  !! number of loop iterations in current running sum
  real (8) :: rem_step_loop_start     !! most recent timestamp for start of step loop
  real (8) :: rem_step_loop_sum = 0.0 !! current step loop cost running sum
  real (8) :: rem_step_loop_max = 0.0 !! max observed step loop cost

  logical  :: rem_restart_any = .false. !! flag indicating whether a restart has ever occurred
  logical  :: rem_restart = .false.     !! flag indicating whether restart occurred this loop count
  real (8) :: rem_restart_write_start   !! most recent timestamp for start of restart write
  real (8) :: rem_restart_write         !! most recent cost for restart write
  real (8) :: rem_restart_write_max = 0.0 !! max observed restart write cost
end module

!! Equilibrium module
module eq_module
  use EZspline_obj
  use EZspline
  character (len=80):: eq_header, eq_filename
  real (8) :: eq_min_r, eq_max_r, eq_min_z, eq_max_z
  real (8) :: eq_axis_r, eq_axis_z, eq_axis_b
  real (8) :: eq_x_psi, eq_x_z, eq_x_r, eq_x_slope
  logical :: eq_set_x2
  real (8) :: eq_x2_z, eq_x2_r, eq_x2_slope
  integer :: eq_mr, eq_mz, eq_mpsi,eq_sep
  real (8), allocatable :: eq_I(:), eq_psi_grid(:),eq_rgrid(:), eq_zgrid(:),eq_psi_rz(:,:)


  ! input variable for profile
  integer :: eq_den_shape, eq_tempi_shape, eq_tempe_shape, eq_flowi_shape, eq_flowe_shape
  real (8) :: eq_den_v1,   eq_den_v2,   eq_den_v3,   eq_den_x1,   eq_den_x2,   eq_den_x3
  real (8) :: eq_tempi_v1, eq_tempi_v2, eq_tempi_v3, eq_tempi_x1, eq_tempi_x2, eq_tempi_x3
  real (8) :: eq_flowi_v1, eq_flowi_v2, eq_flowi_v3, eq_flowi_x1, eq_flowi_x2, eq_flowi_x3
  real (8) :: eq_tempe_v1, eq_tempe_v2, eq_tempe_v3, eq_tempe_x1, eq_tempe_x2, eq_tempe_x3
  real (8) :: eq_flowe_v1, eq_flowe_v2, eq_flowe_v3, eq_flowe_x1, eq_flowe_x2, eq_flowe_x3
  ! Z_eff profile for impurity radiation
  integer :: eq_zeff_shape
  real (8) :: eq_zeff_v1, eq_zeff_v2, eq_zeff_v3, eq_zeff_x1, eq_zeff_x2, eq_zeff_x3




  ! profile interpolation
  character (len=100) :: eq_den_file, eq_tempi_file, eq_tempe_file, eq_flowi_file, eq_flowe_file, eq_zeff_file

  integer, parameter :: eq_mid_r_npsi=50  !!  For poloidal flux -> midplane R function. Number of evaluation points
  real (8) :: eq_mid_r_psi(eq_mid_R_npsi) !! For poloidal flux -> midplane R function. Function data
  real (8) :: eq_mid_r_dp  !! For poloidal flux -> midplane R function. Delta psi

  real (8) , parameter :: epsil_psi =  1D-5
  ! data structure for equilbirum profile
  type eq_ftn_type
     integer :: shape
     real (8):: inx(3), iny(3)
     real (8):: sv(6)
     ! for arbitrary profile function - use pspline
     character (len=100) :: filename
     type (EZspline1_r8) :: spl
     real (8) :: min, max
  end type eq_ftn_type

  type(eq_ftn_type) :: eq_tempi, eq_tempe, eq_den, eq_flowi, eq_flowe, eq_zeff
contains
  ! check if region 1 or 2
  logical function is_rgn12(r,z,psi)
    implicit none
    real (8) :: r,z,psi


    if(psi > eq_x_psi -epsil_psi .or. &
         -(r-eq_x_r)*eq_x_slope + (z-eq_x_z) > 0D0 .and. -(r-eq_x2_r)*eq_x2_slope + (z-eq_x2_z) < 0D0   ) then
       is_rgn12=.true.
    else
       is_rgn12=.false.
    endif
  end function is_rgn12

 ! check if region 1
  logical function is_rgn1(r,z,psi)
    implicit none
    real (8) :: r,z,psi

    if(psi <= eq_x_psi-epsil_psi .and. &
         -(r-eq_x_r)*eq_x_slope + (z-eq_x_z) > 0D0 .and. -(r-eq_x2_r)*eq_x2_slope + (z-eq_x2_z) < 0D0   ) then
       is_rgn1=.true.
    else
       is_rgn1=.false.
    endif
  end function is_rgn1


  !setup eq_ftn functions
  subroutine eq_ftn_setup(ftn)
    implicit none
    type(eq_ftn_type) :: ftn

    select case(ftn%shape)
    case(0) ! constant
       continue  ! do nothing
    case(1) ! hyperbolic tanh
       ftn%sv(1)=0.5D0*(ftn%iny(1)+ftn%iny(2))
       ftn%sv(2)=0.5D0*(ftn%iny(1)-ftn%iny(2))
       ftn%sv(3)=2D0/ftn%inx(2)
    case(2) ! linear
       ftn%sv(1)=(ftn%iny(2)-ftn%iny(1))/ ftn%inx(2) ! slope
       ftn%sv(2)=ftn%iny(1)-ftn%sv(1)*ftn%inx(1)     ! offset
    case(3) ! A exp(B*w tanh( (r-r0)/w )
       ftn%sv(1)=ftn%inx(2)/ftn%inx(3)
       ftn%sv(2)=1D0/ftn%inx(2)
    case(4)  ! modified hyperbolic tanh : sqrt linear + tanh --> with if statement
       ftn%sv(1)=0.5D0*(ftn%iny(1)+ftn%iny(2))
       ftn%sv(2)=0.5D0*(ftn%iny(1)-ftn%iny(2))
       ftn%sv(3)=2D0/ftn%inx(2)
       ftn%sv(4)=(ftn%iny(3)-ftn%iny(1)) / (sqrt(ftn%inx(3)) - sqrt(ftn%inx(1)-0.5D0*ftn%inx(2)))
       ftn%sv(5)= - ftn%sv(4)*sqrt(ftn%inx(1))
    case(5)   ! modified hyperbolic tanh : sqrt linear + tanh --> smooth function
       ! A ( (1+z slope) e^z - e^-z )/(e^z + e^-z) + B
       ! z= 4 ( center - sqrt(center psi) )/ width
       ! iny(1) - edge_val : ped_top
       ! iny(2) - out_val : ped bottom
       ! inx(1) - center ,  inx(2) - width
       ! iny(3) - core(axis)_val , inx(3) - core(axis)_psi
       ftn%sv(1)=0.5D0*(ftn%iny(1)+ftn%iny(2))    !
       ftn%sv(2)=0.5D0*(ftn%iny(1)-ftn%iny(2))    ! height
       ftn%sv(3)=4D0/ftn%inx(2)                   ! width inverse
       ftn%sv(4)=(ftn%iny(3)-ftn%iny(1)) / (sqrt(ftn%inx(3)) - sqrt(ftn%inx(1))) &
            /( - (ftn%sv(2)+1D-99) * ftn%sv(3) * sqrt(ftn%inx(1)) )  ! slope
    case(-1) ! arbitrary profile - file input
#ifdef PSPLINE
       call init_ftn_spline(ftn)
#else
       print *, 'PSPLINE library is required for -1 shape option . Use -DPSPLINE in compiler option'
#endif
    case default
       print *, 'Invalid shape number in eq_ftn_setup', ftn%shape
       stop
    end select
  end subroutine eq_ftn_setup

  ! function evaluation
  real(8) function eq_ftn(p,r,z,ftn)
    use sml_module
    implicit none
    type(eq_ftn_type) :: ftn
    real(8) :: p, r, z, tmp, tmp2, tmp3, tmp4

!   if(z>eq_x_z .or. p > eq_x_psi) then
    if( is_rgn12(r,z,p) ) then
       tmp=p
    else
       tmp=eq_x_psi
    endif

    !tmp=min(max(p,sml_pmin), sml_pmax)

    select case(ftn%shape)
    case(0) ! constant
       eq_ftn=ftn%iny(1)
    case(1) ! hyperbolic tanh
       eq_ftn = ftn%sv(2)*tanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
    case(2) ! linear
       eq_ftn = ftn%sv(1)*tmp + ftn%sv(2)
    case(3) ! a exp(-B w tanh( (r-r0)/w ))
       eq_ftn = ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2)))
    case(4)
       eq_ftn = ftn%sv(2)*dtanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
       if(tmp < ftn%inx(1)-0.5D0*ftn%inx(2) ) eq_ftn = eq_ftn + ftn%sv(4)*sqrt(tmp) + ftn%sv(5)
    case(5)
       tmp2=ftn%sv(3)*( ftn%inx(1)-sqrt(ftn%inx(1)*tmp) )  ! z
       tmp3=exp(tmp2)                           ! expz
       tmp4=exp(-tmp2)                          ! expmz
       ! A * ( (1+z*slope)*expz - expmz )/(expz + expmz) + B
       eq_ftn = ftn%sv(2)*( (1+tmp2*ftn%sv(4))*tmp3 - tmp4 )/(tmp3+tmp4) + ftn%sv(1)
    case(-1)
       tmp=min(max(tmp,ftn%min),ftn%max)
       call ftn_evaluation(ftn%spl,tmp,eq_ftn)
    case default
       print *, 'Invalid shape number in eq_ftn', ftn%shape
       stop
    end select

  end function eq_ftn

  ! function evaluation
  real(8) function eq_dftn(p,r,z,ftn)
    use sml_module
    implicit none
    type(eq_ftn_type) :: ftn
    real(8) :: p, r, z, tmp


    ! region searching and enforcing region 3 value to x-point value
!    if(z>eq_x_z .or. p > eq_x_psi) then
    if(is_rgn12(r,z,p)) then
       tmp=p
    else
       tmp=eq_x_psi
       eq_dftn=0D0
       return
    endif


    select case(ftn%shape)
    case(0) ! constant
       eq_dftn=0
    case(1) ! hyperbolic tanh
       eq_dftn = -ftn%sv(3)*ftn%sv(2)*(1D0/cosh((ftn%inx(1)-tmp)*ftn%sv(3)))**2
    case(2) ! linear
       eq_dftn = ftn%sv(1)
    case(3)
       eq_dftn = - ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2))) &
            /((cosh((ftn%inx(1)-tmp)*ftn%sv(2)))**2*ftn%inx(3))


    case(-1)
       if(ftn%min < tmp .and. tmp < ftn%max) then
          call ftn_derivative(ftn%spl,tmp,1,eq_dftn)
       else
          eq_dftn=0D0
       endif
    case default
       print *, 'Invalid shape number in eq_dftn', ftn%shape
       stop
    end select

    ! case  4 and case 5 does not have derivative function, yet. --> case 5 is required

  end function eq_dftn

  ! using pspline
  subroutine init_ftn_spline(ftn)
    use sml_module
    implicit none
    type(eq_ftn_type) :: ftn
    integer, parameter :: funit=16
    real (8), allocatable :: psi(:), var(:)
    integer :: num, i, flag
    integer :: ierr, BCS1(2)

    ! open file
    open(unit=funit,file=ftn%filename,action='read')
    read(funit,*) num
    if(num<=0) then
       print *, 'error in profile ftn init : invalid number of data', num, ftn%filename
       stop
    endif

    !mem allocation
    allocate(psi(num),var(num))

    ! read data
    do i=1, num
       read(funit,*) psi(i),var(i)
    enddo

    read(funit,*) flag
    if(flag/=-1) then
       print *, 'error in profile ftn init : invalid number of data or ending flag -1 is not set correctly', flag, ftn%filename
       stop
    endif
    close(16)

    if(sml_mype==0) then
       open(unit=funit,file=trim(ftn%filename)//'.used',status='replace')
       do i=1, num
          write(funit,*) psi(i),var(i)
       enddo
    endif

    ! convert normalization
    psi=psi*eq_x_psi

    ! save psi range
    ftn%min=psi(1)
    ftn%max=psi(num)

    ! set size and boundary condition
    BCS1=0
    call ezspline_init(ftn%spl,num,BCS1,ierr)
    call EZspline_error(ierr)

    ! get coefficient
    ftn%spl%x1=psi
    call ezspline_setup(ftn%spl,var,ierr)
    call EZspline_error(ierr)

    deallocate(psi,var)
  end subroutine init_ftn_spline

  ! pspline funtion value
  subroutine ftn_evaluation(spl,psi,val)
    use EZspline_obj
    use EZspline
    implicit none
    type (EZspline1_r8) :: spl
    real (8) :: psi, val
    integer :: ierr

    call ezspline_interp(spl,psi,val,ierr)

  end subroutine ftn_evaluation

  ! pspline derivative
  subroutine ftn_derivative(spl,psi,der,val)
    use EZspline_obj
    use EZspline
    implicit none
    type (EZspline1_r8) :: spl
    real (8) :: psi, val
    integer :: der
    integer :: ierr

    call ezspline_derivative(spl,der,psi,val,ierr)
  end subroutine ftn_derivative
end module eq_module


!! Particle module
module ptl_module
  integer, parameter ::  ptl_nphase=6, ptl_nconst=3, ptl_nphase2=12
  integer, parameter :: ptl_nsp_max=1   ! 1 for single ion species
  integer :: ptl_isp, ptl_nsp

  ! shift_ie communication options. See subroutine shift_ie in
  ! pol_decomp.F90 for description
  integer, parameter :: num_shift_ie_opts              = 10
  integer, parameter :: index_shift_ie_max_nthreads    =  1
  integer, parameter :: index_shift_ie_use_alltoall    =  2
  integer, parameter :: index_shift_ie_use_hs_barrier0 =  3
  integer, parameter :: index_shift_ie_large_limit     =  4
  integer, parameter :: index_shift_ie_handshake       =  5
  integer, parameter :: index_shift_ie_use_hs_barrier1 =  6
  integer, parameter :: index_shift_ie_use_sendrecv    =  7
  integer, parameter :: index_shift_ie_all_sendrecv    =  8
  integer, parameter :: index_shift_ie_use_isend       =  9
  integer, parameter :: index_shift_ie_use_rsend       = 10
  integer, parameter :: use_def_shift_ie_opt           = -2
  integer, parameter :: true_shift_ie_opt              =  1
  integer, parameter :: false_shift_ie_opt             =  0

  real (8) :: ptl_mass(0:ptl_nsp_max), ptl_charge(0:ptl_nsp_max),  ptl_c_m(0:ptl_nsp_max), ptl_c2_2m(0:ptl_nsp_max)
  logical ::  ptl_deltaf_sp(0:ptl_nsp_max)
  type ptl_type
     real (8) :: ph(ptl_nphase) ! 1-r, 2-z, 3-zeta, 4-rho_parallel
     real (8) :: ct(ptl_nconst) ! 1?- mu, 2-w0, 3-f0
     integer (8) :: gid
#ifdef PURE_RK4
     real (8) :: dph(ptl_nphase),dpht(ptl_nphase)
#endif
  end type ptl_type

  type species_type
     integer :: num  !! number of particle for each processor
     integer :: min_max_num !! maximum number of particles per process after most recent load balancing
     integer :: maxnum !! upper bound on number of particles per process
     integer :: type !! species type number 1 for ion, 2 for electron
     real (8) :: mass   !! Ion mass (Atomic Unit)
     real (8) :: charge !! Ion charge (electron charge unit)
     real (8) :: c_m
     integer (8) ::  maxgid !! Largest global ID number used.

     ! shift_ie communication options for this particle type
     integer :: shift_opt(num_shift_ie_opts)

     ! phase variables -- make one set of array for cache performance
     type(ptl_type), allocatable :: ptl(:)

     ! rk4 variables
     real (8), allocatable :: phase0(:,:) !! Ion phase variable - storage for old phse variable

     ! lost particle save
     integer :: lost_num, lost_nummax
     real (8), allocatable :: lost_index(:)

     ! particle-grid position save
     integer, allocatable :: tr_save(:)  !size of electron and ion are different
     real (8), allocatable :: p_save(:,:) ! size of electron and ion are different


     ! save rhoi
     real (8), allocatable :: rhoi(:) ! rho_perp of each particle
     real (8), allocatable :: rhon(:,:) !Normalized rho on grid
  end type species_type


  ! electron subcycling
  integer :: ptl_enum_save, ptl_tracer_n_save
  type(ptl_type), allocatable :: ptl_ephase_save(:)

  !special simulation purpose : such as single particle simulation
  real (8) :: ptl_special_r, ptl_special_z,  ptl_special_phi, ptl_special_en_ev, ptl_special_pitch
  integer :: ptl_special_n
  integer :: ptl_maxnum, ptl_e_maxnum,ptl_num, ptl_e_num, ptl_lost_nummax
  real (8) :: ptl_mass_au, ptl_charge_eu, ptl_e_mass_au, ptl_e_charge_eu


  integer, parameter :: pir=1, piz=2, pip=3, pirho=4, pim=1, piw1=5, piw2=6, piw0=2, pif0=3
  contains
    subroutine ptl_mem_allocation( sp, sp_type, maxnum,mass, charge,  nlarmor, lost_nummax)
      implicit none
      integer, intent(in) :: sp_type, maxnum, nlarmor, lost_nummax
      real (8), intent(in) :: mass, charge
      type(species_type) :: sp

      integer :: np

      sp%maxnum=maxnum
      np=ptl_nphase
      sp%mass=mass
      sp%charge=charge
      sp%c_m=charge/mass

      sp%lost_nummax=lost_nummax
      sp%lost_num=0

      allocate( sp%ptl(maxnum), sp%phase0(np,maxnum), &
           sp%lost_index(lost_nummax) )

      allocate( sp%tr_save(maxnum),sp%p_save(3,maxnum) )

      if(sp_type==1) then !for ion
         sp%type=1
         allocate(sp%rhoi(maxnum))
         allocate(sp%rhon(2,maxnum))
      else  ! for electron
         sp%type=0
         ! for subcycling
         allocate(ptl_ephase_save(maxnum))
      endif

    end subroutine ptl_mem_allocation

    subroutine ptl_mem_reallocation( sp, maxnum, lost_nummax )
      implicit none
      integer, intent(in) :: maxnum, lost_nummax
      type(species_type) :: sp

      ! local variables
      integer :: i, ierror
      type(ptl_type), allocatable :: tmp_ptl(:)
      type(ptl_type), allocatable :: tmp_ptl_ephase_save(:)
      real (8),       allocatable :: tmp_phase0(:,:)
      real (8),       allocatable :: tmp_lost_index(:)
      real (8),       allocatable :: tmp_p_save(:,:)
      integer,        allocatable :: tmp_tr_save(:)
      real (8),       allocatable :: tmp_rhoi(:)

      if (maxnum > sp%maxnum) then
         ! ptl
         allocate  ( tmp_ptl(sp%maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            tmp_ptl(i)%ph(:) = sp%ptl(i)%ph(:)
            tmp_ptl(i)%ct(:) = sp%ptl(i)%ct(:)
            tmp_ptl(i)%gid   = sp%ptl(i)%gid
         enddo
         deallocate( sp%ptl )

         allocate  ( sp%ptl(maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            sp%ptl(i)%ph(:) = tmp_ptl(i)%ph(:)
            sp%ptl(i)%ct(:) = tmp_ptl(i)%ct(:)
            sp%ptl(i)%gid   = tmp_ptl(i)%gid
         enddo
         deallocate( tmp_ptl )

         ! phase0
         allocate  ( tmp_phase0(ptl_nphase,sp%maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            tmp_phase0(:,i) = sp%phase0(:,i)
         enddo
         deallocate( sp%phase0 )

         allocate  ( sp%phase0(ptl_nphase,maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            sp%phase0(:,i) = tmp_phase0(:,i)
         enddo
         deallocate( tmp_phase0 )

         ! tr_save
         allocate  ( tmp_tr_save(sp%maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            tmp_tr_save(i) = sp%tr_save(i)
         enddo
         deallocate( sp%tr_save )

         allocate  ( sp%tr_save(maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            sp%tr_save(i) = tmp_tr_save(i)
         enddo
         deallocate( tmp_tr_save )

         ! p_save
         allocate  ( tmp_p_save(3,sp%maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            tmp_p_save(:,i) = sp%p_save(:,i)
         enddo
         deallocate( sp%p_save )

         allocate  ( sp%p_save(3,maxnum), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%maxnum
            sp%p_save(:,i) = tmp_p_save(:,i)
         enddo
         deallocate( tmp_p_save )

         if(sp%type==1) then
            ! rhoi (for ion)
            allocate  ( tmp_rhoi(sp%maxnum), stat=ierror )
            if (ierror /= 0) then
               write(6,*) "allocation error in ptl_mem_reallocation"
!pw            call flush(6)
            endif
            do i=1,sp%maxnum
               tmp_rhoi(i) = sp%rhoi(i)
            enddo
            deallocate( sp%rhoi )

            allocate  ( sp%rhoi(maxnum), stat=ierror )
            if (ierror /= 0) then
               write(6,*) "allocation error in ptl_mem_reallocation"
!pw            call flush(6)
            endif
            do i=1,sp%maxnum
               sp%rhoi(i) = tmp_rhoi(i)
            enddo
            deallocate( tmp_rhoi )
         else
            ! ptl_ephase_save (for electron subcycling)
            allocate  ( tmp_ptl_ephase_save(sp%maxnum), stat=ierror )
            if (ierror /= 0) then
               write(6,*) "allocation error in ptl_mem_reallocation"
!pw            call flush(6)
            endif
            do i=1,sp%maxnum
               tmp_ptl_ephase_save(i)%ph(:) = ptl_ephase_save(i)%ph(:)
               tmp_ptl_ephase_save(i)%ct(:) = ptl_ephase_save(i)%ct(:)
               tmp_ptl_ephase_save(i)%gid   = ptl_ephase_save(i)%gid
            enddo
            deallocate( ptl_ephase_save )

            allocate  ( ptl_ephase_save(maxnum), stat=ierror )
            if (ierror /= 0) then
               write(6,*) "allocation error in ptl_mem_reallocation"
!pw            call flush(6)
            endif
            do i=1,sp%maxnum
               ptl_ephase_save(i)%ph(:) = tmp_ptl_ephase_save(i)%ph(:)
               ptl_ephase_save(i)%ct(:) = tmp_ptl_ephase_save(i)%ct(:)
               ptl_ephase_save(i)%gid   = tmp_ptl_ephase_save(i)%gid
            enddo
            deallocate( tmp_ptl_ephase_save )
         endif

         sp%maxnum=maxnum
      endif

      if (lost_nummax > sp%lost_nummax) then
         ! lost_index
         allocate  ( tmp_lost_index(sp%lost_nummax), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%lost_nummax
            tmp_lost_index(i) = sp%lost_index(i)
         enddo
         deallocate( sp%lost_index )

         allocate  ( sp%lost_index(lost_nummax), stat=ierror )
         if (ierror /= 0) then
            write(6,*) "allocation error in ptl_mem_reallocation"
!pw         call flush(6)
         endif
         do i=1,sp%lost_nummax
            sp%lost_index(i) = tmp_lost_index(i)
         enddo
         deallocate( tmp_lost_index )

         sp%lost_nummax=lost_nummax
      endif

    end subroutine ptl_mem_reallocation

end module ptl_module


!! Module for Field value
Module fld_module
  type fld_type
    real (8) :: r, z, phi  !! NOT USED. R,Z coordinate variable
    real (8) :: I         !! R*Bt
    real (8) :: epot      !! Electric potential
    real (8) :: br, bz, bphi  ! B_r, B_z, B_phi
    real (8) :: dbrdr, dbrdz, dbrdp, dbzdr, dbzdz, &
         dbzdp, dbpdr, dbpdz, dbpdp !! dB_x / dx , x=R,Z,phi
    real (8)   :: Er, Ez, Ephi  !! E_r, E_z, E_phi
    real (8) :: Er00, Ez00 ! only for electron

    !electron weight calculation
    real (8) :: ddpotdt   ! dPhi/dt -- for adiabatic term

    !for weight calculation
    real (8) :: dIdpsi, dpsidr,dpsidz !! dI/dpsi, dpsi/dr, dpsi/dz for weight calculation
    real (8) :: psi  !! Poloidal Flux
    !for weight calculation : f0 information
    !real (8) :: f0_den,f0_gradn(2),f0_temp,f0_gradt(2)
 end type fld_type
end module


!! Interpolation module
module itp_module
#if defined(PSPLINE)
  use EZspline_obj
  type(EZspline2_r8) :: spl(0:2,0:2)
  type(EZspline1_r8) :: spl_psi
#endif


  integer :: itp_mpsi, itp_mr, itp_mz
  integer, parameter :: itp_mr2=100, itp_mz2=100
  integer,parameter :: itp_korder=3, itp_korder_rz=5, itp_korder_rz2=5
  real (8) :: itp_min_r, itp_max_r, itp_min_z, itp_max_z
  real (8) :: itp_min_psi, itp_max_psi

  real (8), allocatable :: itp_I_cscoef(:,:), itp_I_break(:)
  real (8), allocatable :: itp_psi_knot(:), itp_r_knot(:), itp_z_knot(:)
  real (8), allocatable :: itp_psi_bscoef(:)

  real (8) :: itp_rgrid2(itp_mr2), itp_zgrid2(itp_mz2)
  real (8) :: itp_r_knot2(itp_mr2+itp_korder_rz2), itp_z_knot2(itp_mz2+itp_korder_rz2)

  real (8) :: itp_psi_bscoef00(itp_mr2*itp_mz2), itp_psi_bscoef01(itp_mr2*itp_mz2), &
       itp_psi_bscoef10(itp_mr2*itp_mz2), itp_psi_bscoef02(itp_mr2*itp_mz2), &
       itp_psi_bscoef11(itp_mr2*itp_mz2), itp_psi_bscoef20(itp_mr2*itp_mz2)
end module

!! Diagnosis module
module diag_module
  integer, parameter :: diag_max_sp_num=2  ! two species

  ! tracer
  integer :: diag_tracer_period !! Period for tracer
  integer :: diag_tracer_n     !! Particle index for tracer routine
  integer :: diag_tracer_sp     ! Species index for tracer

  ! diag_particle
  integer :: diag_particle_mod       ! # of particle to be saved = total particle / diag_particle_mo
  integer :: diag_particle_period  ! how often write out

  ! for 1D
  integer, parameter :: diag_1d_npv1=15

  logical  :: diag_1d_on
  integer  :: diag_1d_period
  integer  :: diag_1d_npsi  !
  real (8) :: diag_1d_pin, diag_1d_pout
  real (8) :: diag_1d_dp, diag_1d_dp_inv
  real (8),allocatable :: diag_1d_vol(:)
  integer :: diag_1d_isp,diag_1d_nsp  ! same as ptl_isp/nsp
  logical  :: diag_tavg_on

  real (8), allocatable :: diag_1d_f_pv1(:,:,:,:)
  real (8), allocatable :: diag_1d_df_pv1(:,:,:,:)
  real (8), allocatable :: diag_1d_f0_pv1(:,:,:,:)
  real (8), allocatable :: diag_1d_tavg_f_pv1(:,:,:,:)
  real (8), allocatable :: diag_1d_tavg_df_pv1(:,:,:,:)
  real (8), allocatable :: diag_2d_ef0_pv(:,:,:,:,:)

  !outer midplane
  logical :: diag_omid_on
  real (8), allocatable :: diag_1d_omid_f_pv1(:,:,:,:)


  ! for 3D
  logical :: diag_3d_on
  integer :: diag_3d_period
  logical :: diag_3d_more
  integer, parameter :: diag_3d_nvar=3
  real (8), allocatable :: diag_3d_add(:,:,:)

  ! for f0
  integer :: diag_f0_period
  logical :: diag_f0_g_on

  logical :: diag_f0_df_on
  integer, parameter :: diag_f0_df_npv1=5  ! number of variables to be stored.
                                        ! (1) volume, (2) density, (3) density change, (4) energy change, (5) momentum change
  integer, parameter :: diag_f0_df_nsource=4   ! number of type of sources
                                               ! (1) collision, (2) heat_torqe, (3) neutral, (4) radiation
  real (8), allocatable :: diag_f0_df_pv1(:,:,:,:) ! 


  ! for neutral
  integer :: diag_neu_npsi
  real (8) :: diag_neu_pin, diag_neu_pout
  real (8) :: diag_neu_dp_inv
  integer, parameter  :: diag_neu_nv=4
  real (8), allocatable :: diag_neu_port(:,:,:), diag_neu_vol(:)

  ! for ion flux
  integer :: diag_1d_ne
  logical :: diag_eflux_on
  real (8) :: diag_1d_emin, diag_1d_emax
  real (8), allocatable :: diag_1d_eflux_pv(:,:,:,:,:)
  real (8), allocatable :: diag_2d_dflux_pv(:,:,:,:,:)

  !for heat
  logical :: diag_heat_on
  integer, parameter :: diag_heat_nvar=5
  integer :: diag_heat_nsection
  real (8) :: diag_heat_rmax1, diag_heat_rmin1, diag_heat_zmax1, diag_heat_zmin1
  real (8) :: diag_heat_rmax2, diag_heat_rmin2, diag_heat_zmax2, diag_heat_zmin2
  real (8) :: diag_heat_rmax3, diag_heat_rmin3, diag_heat_zmax3, diag_heat_zmin3
  real (8) :: diag_heat_rmax(3), diag_heat_rmin(3), diag_heat_zmax(3), diag_heat_zmin(3), diag_heat_dr(3), diag_heat_dz(3)
  real (8) :: diag_heat_pmax(3), diag_heat_pmin(3), diag_heat_dp(3)
  integer :: diag_heat_nr, diag_heat_nz, diag_heat_npsi
  real (8), allocatable :: diag_heat_pv(:,:,:,:,:,:),diag_heat_pv_psi(:,:,:,:,:)
end module diag_module

!****************************************************************************
! simulation parameter for collision module definition
!
! first created : 2000/10/26
! last modified : 2011/06/01
!****************************************************************************
!! Collision module
module col_module
  integer :: col_mode
  integer :: col_period
  integer, parameter :: col_imp=0 ! imp=1 impurity species, 0=none
  logical :: col_moving_frame

  real (8) :: col_pin, col_pout ! now those are used for col_f
  integer :: col_f_start ! start time of collisions

  !2002/05/15 -- col mode 2 variables
  integer,parameter :: col_2_m=50 !! Number of radial slice for conserving collision
  integer,parameter :: col_2_mtheta=8 !! Number of poloidal slice for conserving collision
  real (kind=8) :: col_2_pin, col_2_pout !! Inner and outer boundary
  real (kind=8) :: col_2_dp      !! Delta psi of a cell for conserving collision
  real (kind=8) :: col_2_dtheta  !! Delta theta of a cell for conserving collision
  real (kind=8) :: col_2_inv_dp, col_2_inv_dtheta

  real (kind=8), allocatable  :: col_2_vol(:,:) !! Volume of each cell for conserving collision
  real (kind=8), allocatable :: col_2_delta_vsum_diag(:,:,:,:) !!for conservation accuracy check
  real (kind=8), allocatable :: col_2_w0(:,:) !! for delta-f collision

  integer :: col_2_iteration_method=2  !! 1. iteration method    2. matrix method

  !2004/01/16 -- impurity collision
  integer :: col_imp_on=0  !!Switch for impurity collision
  real (kind=8) :: col_imp_charge, col_imp_mass  !! Charge(electorn) and Mass(AU) of impurity
  real (kind=8) :: col_den_imp_edge,col_den_imp_out,col_den_imp_ped_c, col_den_imp_ped_width  !! Profile information of impurity ions. Tanh model.

  !2002/08/30 -- varying background
  integer :: col_varying_bg !!Switch for background plasma update
  integer :: col_vb_m=50, col_vb_mtheta=8       !! Number of background update points
  integer :: col_vb_period  !! Period of background updating
  real (kind=8) :: col_vb_pin, col_vb_pout !! Inner and outer boundary and delta psi for background update
  real (kind=8) :: col_vb_dp, col_vb_dtheta, col_vb_inv_dp, col_vb_inv_dtheta
  real (kind=8) , allocatable :: col_vb(:,:,:,:), col_vb_vol(:,:) !! Density, Temperature and velocity of each shell for backgournd update

  !2002/10/09
  integer :: col_en_col_on   !! Switch for energy collision

  !2008/10/31
  logical  :: col_accel  ! artificial collision amplifying
  real (8) :: col_accel_pin1, col_accel_pout1, col_accel_factor1
  real (8) :: col_accel_pin2, col_accel_pout2, col_accel_factor2
  integer :: col_accel_n

  !conserving collision
  integer, allocatable :: col_org_j(:),col_org_l(:)
  real (kind=8),allocatable :: col_org_b(:),col_coeff_a(:), col_coeff_b(:), col_coeff_d(:),&
                               col_org_nuE(:),col_org_nuD(:),col_org_weight(:),&
                               col_org_vp_m(:), col_org_v2(:)
  logical, allocatable :: col_gid(:)

  !2013-02-01 -- col mode 3 : vpic
  integer :: mpi_ptl_type
  integer :: col_3_npsi, col_3_ntheta         ! The number of volumes in psi and theta, respectively
  integer :: col_3_npsi_sol, col_3_ntheta_sol ! The number of volumes in psi and theta in "SOL"
  integer :: col_3_nvr, col_3_nvz      ! The number of "GRIDS" in vr and vz, respectively
  integer :: col_3_ntotal_r, col_3_ntotal_v
  integer :: col_3_nset, col_3_nset_mod, col_3_nset_max, col_3_nset_allocated
  real (kind=8) :: col_3_dpsi_core, col_3_dtheta_core, col_3_dpsi_sol, col_3_dtheta_sol
  real (kind=8) :: col_3_dt, col_3_theta_offset
  integer ,allocatable :: col_3_tag(:)
  !2013-02-06 --- col_3_decomp
  integer, allocatable :: col_3d_scounts(:), col_3d_rcounts(:), col_3d_sdispls(:), col_3d_rdispls(:)
  real (kind=8), allocatable  :: col_3_vol_r(:)  !! Volume of each cell for conserving collision
  !2013-02-16 --- col_3_subgrid
  integer, allocatable, dimension(:) :: col_3_subg_num, col_3_refg_indx
  real (kind=8), allocatable, dimension(:) :: col_3_subg_dth
  integer :: col_3_total_subg_num
  logical :: col_3_sol_solve
  ! required minimum number of particles
  integer :: col_3_min_popul=100
  !2013-02-21 --- col_3_subgrid_correction due to limiter
  integer, allocatable, dimension(:) :: col_3_grid_valid_list
  !2013-02-23 --- SUPER LU
  integer, allocatable, dimension(:,:) :: index_map_LU
  integer, allocatable, dimension(:) :: LU_cvalues, LU_rowindx, LU_colptr
  integer :: LU_n, LU_nnz, LU_nrhs, LU_ldb

  !2013-02-07 --- col_3_type
  type col_3_type
      integer, allocatable, dimension(:) :: num
      real (kind=8),allocatable, dimension(:,:,:) :: col_3_vmap, col_3_vmap_prev !col_3_vmap_prev : for debug (save vmap before collision)
      !col_3_nset_allocated size
      real (kind=8), allocatable, dimension(:) :: lz, lr, Vel_avg_z, deltax, deltay
      real (kind=8), allocatable, dimension(:,:) :: Evolume_inV_proc, vol
      real (kind=8), allocatable, dimension(:) :: den, t_ev
      real (kind=8), allocatable, dimension(:) :: org_weight_sum, org_v2_avg
      real (kind=8), allocatable, dimension(:) :: negative_weight_sum, negative_v2_sum, negative_vp_sum
      integer :: negative_count
      integer, allocatable, dimension(:) :: negative_tag

      !ptl num size
      integer, allocatable, dimension (:) :: ipsth, ivel_r, ivel_z
      real (kind=8), allocatable, dimension(:) :: vpic_ptl_v_parallel, vpic_ptl_v_perp
      real (kind=8), allocatable, dimension(:,:) :: p2m_fc
  end type col_3_type
  type col_3_core_type
      real (kind=8) :: numeric_T, numeric_Teq, numeric_vth2
      real (kind=8) :: mass, mesh_dr, mesh_dz, dens, ens
      real (kind=8), allocatable, dimension(:) :: mesh_r, mesh_r_half, mesh_z, mesh_z_half, local_center_volume, vol
      real (kind=8), allocatable, dimension(:,:) :: delta_r, delta_z
  end type col_3_core_type

#ifdef COL_NO_EDGE_CONSERV
  !col routine for no edge collision conservation part
  integer :: col_no_conserv_m_num=5
#endif

contains
  subroutine col_accel_factor(psi,factor)
    implicit none
    real (8) :: psi, factor

    factor=1D0
    if(col_accel_n>=1) then
       if(col_accel_pin1 < psi .and. psi < col_accel_pout1 ) then
          factor=col_accel_factor1*factor
       endif
    endif
    if(col_accel_n>=2) then
       if(col_accel_pin2 < psi .and. psi < col_accel_pout2 ) then
          factor=col_accel_factor2*factor
       endif
    endif
    return
  end subroutine col_accel_factor

  subroutine init_col_module(isp,nsp,sp)
    use sml_module, only : sml_2pi, sml_mype, sml_plane_mype, sml_pe_per_plane
    use ptl_module
    implicit none
    integer, intent(in) :: isp, nsp
    type(species_type):: sp
    integer :: i, j, k, sep, itr
    real (kind=8) :: tri_c(2), tri_c_r, tri_c_theta
    integer, allocatable :: node_counter(:)

    if(col_varying_bg==1) then
        allocate(col_vb_vol(col_vb_mtheta,0:col_vb_m))
        allocate(col_vb(3,col_vb_mtheta,0:col_vb_m,isp:nsp))
    endif

    if(col_mode==2) then
        allocate(col_2_vol(col_2_mtheta, 0:col_2_m))
        allocate(col_org_j(sp%maxnum), col_org_l(sp%maxnum), col_org_b(sp%maxnum),&
                 col_coeff_a(sp%maxnum), col_coeff_b(sp%maxnum), col_coeff_d(sp%maxnum), &
                 col_org_nuE(sp%maxnum), col_org_nuD(sp%maxnum), &
                 col_org_weight(sp%maxnum), col_org_vp_m(sp%maxnum), &
                 col_org_v2(sp%maxnum), col_gid(sp%maxnum))
    endif

    if(col_mode==3) then
        allocate(col_3_tag(0:col_3_nset_allocated-1))
        allocate(col_3_vol_r(0:col_3_ntotal_r-1))

        ! Allocation : collision cells =>> CPU_destination
        sep = col_3_nset_max*col_3_nset_mod
        k=0
        do i=0, col_3_ntotal_r-1
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

        !! Allocation : nodes =>> collision cells
        !! 1. count for dynamic array
        !!  parallelization can be possible
        !allocate(node_counter(0:col_3_ntotal_r-1))  !!For local calculation
        !node_counter = 0
        !do itr=1,grid%ntriangle
        !    tri_c = (grid%x(:,grid%nd(1,itr))+grid%x(:,grid%nd(2,itr))+grid%x(:,grid%nd(3,itr)))/3D0
        !    tri_c_r = sqrt((tri_c(1)-eq_axis_r)*(tri_c(1)-eq_axis_r) + (tri_c(2)-eq_axis_z)*(tri_c(2)-eq_axis_z))
        !    tri_c_theta= acos((tri_c(1)-eq_axis_r)/tri_c_r)
        !    if(tri_c(2) .lt. eq_axis_z) then
        !        tri_c_theta= sml_2pi-tri_c_theta
        !    endif
        !
        !    if(grid%rgn(itr) .eq. 1) then
        !        i=floor(grid%psi(itr)/col_3_dpsi_core)*col_3_ntheta+tri_c_theta/col_3_dtheta_core
        !        call assert(i .lt. col_3_npsi*col_3_ntheta, 'logical error for indexing of col_3',j)
        !    else !SOL including private regions
        !        i=col_3_npsi*col_3_ntheta+ floor((grid%psi(itr)-eq_x_psi)/col_3_dpsi_sol)*col_3_ntheta+tri_c_theta/col_3_dtheta_sol
        !        call assert(i .ge. col_3_npsi*col_3_ntheta .and. i .lt. col_3_ntotal_r, 'logical error for indexing of col_3',j)
        !    endif
        !    node_counter(i) = node_counter(i)+1
        !enddo

        !! 2. save node index for me
        !! parallelization can be poissible
        !allocate(col_3_mynode(max(node_counter(col_3_tag)),col_3_nset_allocated))  !global
        !node_counter = 0
        !do itr=1,grid%ntriangle
        !    tri_c = (grid%x(:,grid%nd(1,itr))+grid%x(:,grid%nd(2,itr))+grid%x(:,grid%nd(3,itr)))/3D0
        !    tri_c_r = sqrt((tri_c(1)-eq_axis_r)*(tri_c(1)-eq_axis_r) + (tri_c(2)-eq_axis_z)*(tri_c(2)-eq_axis_z))
        !    tri_c_theta= acos((tri_c(1)-eq_axis_r)/tri_c_r)
        !    if(tri_c(2) .lt. eq_axis_z) then
        !        tri_c_theta= sml_2pi-tri_c_theta
        !    endif
        !
        !    if(grid%rgn(itr) .eq. 1) then
        !        i=floor(grid%psi(itr)/col_3_dpsi_core)*col_3_ntheta+tri_c_theta/col_3_dtheta_core
        !        call assert(i .lt. col_3_npsi*col_3_ntheta, 'logical error for indexing of col_3',j)
        !    else !SOL including private regions
        !        i=col_3_npsi*col_3_ntheta+ floor((grid%psi(itr)-eq_x_psi)/col_3_dpsi_sol)*col_3_ntheta+tri_c_theta/col_3_dtheta_sol
        !        call assert(i .ge. col_3_npsi*col_3_ntheta .and. i .lt. col_3_ntotal_r, 'logical error for indexing of col_3',j)
        !    endif
        !
        !    do k=0:col_3_nset_allocated-1
        !        if(i .eq. col_3_tag(k)) then
        !            node_counter(k) = node_counter(k)+1
        !            col_3_mynode(node_counter(k),k) = i
        !            exit
        !        endif
        !    enddo
        !enddo
        !deallocate(node_counter)

    endif

  end subroutine init_col_module

  subroutine col2_mem_reallocation(old_maxnum, new_maxnum)
    implicit none
    integer, intent(in) :: old_maxnum, new_maxnum

    ! local
    integer :: i, ierror
    integer, allocatable      :: tmp_int(:)
    real (kind=8),allocatable :: tmp_real8(:)
    logical, allocatable      :: tmp_logical(:)

    if ((new_maxnum > old_maxnum) .and. (col_mode==2)) then

       ! integer fields
       allocate  ( tmp_int(old_maxnum), stat=ierror )

       ! col_org_j
       do i=1,old_maxnum
          tmp_int(i) = col_org_j(i)
       enddo
       deallocate( col_org_j )
       allocate  ( col_org_j(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_j(i) = tmp_int(i)
       enddo

       ! col_org_l
       do i=1,old_maxnum
          tmp_int(i) = col_org_l(i)
       enddo
       deallocate( col_org_l )
       allocate  ( col_org_l(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_l(i) = tmp_int(i)
       enddo

       deallocate( tmp_int )

       ! real8 fields
       allocate  ( tmp_real8(old_maxnum), stat=ierror )

       ! col_org_b
       do i=1,old_maxnum
          tmp_real8(i) = col_org_b(i)
       enddo
       deallocate( col_org_b )
       allocate  ( col_org_b(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_b(i) = tmp_real8(i)
       enddo

       ! col_coeff_a
       do i=1,old_maxnum
          tmp_real8(i) = col_coeff_a(i)
       enddo
       deallocate( col_coeff_a )
       allocate  ( col_coeff_a(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_coeff_a(i) = tmp_real8(i)
       enddo

       ! col_coeff_b
       do i=1,old_maxnum
          tmp_real8(i) = col_coeff_b(i)
       enddo
       deallocate( col_coeff_b )
       allocate  ( col_coeff_b(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_coeff_b(i) = tmp_real8(i)
       enddo

       ! col_coeff_d
       do i=1,old_maxnum
          tmp_real8(i) = col_coeff_d(i)
       enddo
       deallocate( col_coeff_d )
       allocate  ( col_coeff_d(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_coeff_d(i) = tmp_real8(i)
       enddo

       ! col_org_nuE
       do i=1,old_maxnum
          tmp_real8(i) = col_org_nuE(i)
       enddo
       deallocate( col_org_nuE )
       allocate  ( col_org_nuE(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_nuE(i) = tmp_real8(i)
       enddo

       ! col_org_nuD
       do i=1,old_maxnum
          tmp_real8(i) = col_org_nuD(i)
       enddo
       deallocate( col_org_nuD )
       allocate  ( col_org_nuD(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_nuD(i) = tmp_real8(i)
       enddo

       ! col_org_weight
       do i=1,old_maxnum
          tmp_real8(i) = col_org_weight(i)
       enddo
       deallocate( col_org_weight )
       allocate  ( col_org_weight(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_weight(i) = tmp_real8(i)
       enddo

       ! col_org_vp_m
       do i=1,old_maxnum
          tmp_real8(i) = col_org_vp_m(i)
       enddo
       deallocate( col_org_vp_m )
       allocate  ( col_org_vp_m(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_vp_m(i) = tmp_real8(i)
       enddo

       ! col_org_v2
       do i=1,old_maxnum
          tmp_real8(i) = col_org_v2(i)
       enddo
       deallocate( col_org_v2 )
       allocate  ( col_org_v2(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_org_v2(i) = tmp_real8(i)
       enddo

       deallocate( tmp_real8 )

       ! logical fields
       allocate  ( tmp_logical(old_maxnum), stat=ierror )

       ! col_gid
       do i=1,old_maxnum
          tmp_logical(i) = col_gid(i)
       enddo
       deallocate( col_gid )
       allocate  ( col_gid(new_maxnum), stat=ierror )
       do i=1,old_maxnum
          col_gid(i) = tmp_logical(i)
       enddo

       deallocate( tmp_logical )
    endif

  end subroutine col2_mem_reallocation

end module col_module

module bnc_module
  integer, parameter :: bnc_nr=100
  real (8) :: bnc_min_r, bnc_max_r, bnc_dr
  real (8) :: bnc_z_psi_min(bnc_nr)
  real (8) :: bnc_arg_pass_r
end module bnc_module

!----------------------PERFORMANCE STUFF------------------------------------------

module perf_monitor
#if defined(CAM_TIMERS)
  use perf_mod, only: t_startf, t_stopf, t_adj_detailf
#endif
  integer, parameter :: mon_NX=44
  integer, parameter :: mon_N=19
  integer, parameter :: mon_N2=44

!PETSc log events (all from main.F90)
!level 0
  integer, parameter :: TOTAL_           =  1
!level 1
  integer, parameter :: INIT_            =  2
  integer, parameter :: FIRST_           =  3
  integer, parameter :: MAIN_LOOP_       =  4
!level 2
  integer, parameter :: COLLISION_       =  5
  integer, parameter :: SOURCE_          =  6
  integer, parameter :: RESTART_WRITE_   =  7
!level 3
  integer, parameter :: PUSH_I_          =  8
  integer, parameter :: DIAGNOSIS_       =  9
  integer, parameter :: SHIFT_I_         = 10
  integer, parameter :: NEUTRAL_         = 11
  integer, parameter :: CHARGEI_         = 12
  integer, parameter :: POISSON_         = 13
!level 4
  integer, parameter :: POISSON_E1_      = 14
!level 5
  integer, parameter :: PUSHE_           = 15
  integer, parameter :: SHIFT_E_         = 16
  integer, parameter :: CHARGEE_         = 17
  integer, parameter :: POISSON_E2_      = 18
  integer, parameter :: F0_GRID_         = 19
!
!Other main.F90 log events
!level 1 (FINALIZE = TOTAL- (INIT+FIRST+MAIN_LOOP))
  integer, parameter :: FINALIZE_        = 24
!level 2 (IPC_LOOP = PUSH_I+DIAGNOSIS+SHIFT_I+CHARGEI+ELECTRONS+POISSON+GET_POT_GRAD)
  integer, parameter :: PETSC_INIT_      = 25
  integer, parameter :: ADIOS_INIT_      = 26
  integer, parameter :: SETUP_           = 27
  integer, parameter :: INIT_PUSHMOD_GPU_= 28
  integer, parameter :: SET_WEIGHTS_F_   = 29
  integer, parameter :: SHIFT_I_F_       = 30
  integer, parameter :: SHIFT_E_F_       = 31
  integer, parameter :: POISSON_F_       = 32
  integer, parameter :: GET_POT_GRAD_F_  = 33
  integer, parameter :: MAIN_LOOP_RED_   = 34
  integer, parameter :: SET_WEIGHTS_     = 35
  integer, parameter :: IPC_LOOP_        = 36
!level 3 (ELECTRONS = POISSON_E1+HYBRID_LOOP)
  integer, parameter :: ELECTRONS_       = 37
  integer, parameter :: GET_POT_GRAD_    = 38
  integer, parameter :: MEM_CLEAN_RSTRT_ = 39
!level 4 (HYBRID_LOOP = DDPOTDT_UPD+SOL_ELEC_PHASE+ELECTRON_LOOP+CHARGEE+POISSON_E2)
  integer, parameter :: MEM_CLEAN_I_     = 40
  integer, parameter :: HYBRID_LOOP_     = 41
!level 5 (ELECTRON_LOOP = CHARGEE_SEARCH+PUSH_E+SHIFT_E)
  integer, parameter :: DDPOTDT_UPD_     = 42
  integer, parameter :: SOL_ELEC_PHASE_  = 43
  integer, parameter :: GAT_FIELD_INFO_  = 44

#if !defined(NO_PETSC)
  integer :: event(mon_N)
#else
  real (8) :: mon_time(mon_N)
  real (8) :: mon_sum(mon_N)
#endif
  character (len=21) :: mon_str (2*mon_NX) ! must have 2*mon_NX > mon_N (yuck)
  logical            :: mon_sync(2*mon_NX)
  integer            :: mon_adj(2*mon_NX)
  integer            :: mon_flush_count
  integer            :: mon_flush_freq

contains
  subroutine init_perf_monitor()
    use sml_module
#if defined(CAM_TIMERS)
    use perf_mod, only: t_barrier_onf, t_initf
#endif
    implicit none
    integer i
    logical masterproc, exist
#if ( defined TIMING_BARRIERS )
    logical :: timing_barrier = .true.
#else
    logical :: timing_barrier = .false.
#endif
    character (len=21) :: ctemp1, ctemp2
    include 'mpif.h'

#if defined(CAM_TIMERS)
    if (sml_mype == 0) then
       masterproc = .true.
    else
       masterproc = .false.
    endif
    inquire(FILE="perf_in",EXIST=exist)
    if (exist) then
       call t_initf("perf_in", LogPrint=masterproc, &
                    Mpicom=sml_comm, &
                    MasterTask=masterproc)
    else
       call t_initf("input", LogPrint=masterproc, &
                    Mpicom=sml_comm, &
                    MasterTask=masterproc)
    endif
    timing_barrier = t_barrier_onf()
#endif

! PETSc events
    mon_str(TOTAL_)            = 'TOTAL                '
    mon_str(INIT_)             = 'INIT                 '
    mon_str(FIRST_)            = 'FIRST                '
    mon_str(MAIN_LOOP_)        = 'MAIN_LOOP            '
    mon_str(COLLISION_)        = 'COLLISION            '
    mon_str(SOURCE_)           = 'SOURCE               '
    mon_str(RESTART_WRITE_)    = 'RESTART_WRITE        '
    mon_str(PUSH_I_)           = 'PUSH_I               '
    mon_str(DIAGNOSIS_)        = 'DIAGNOSIS            '
    mon_str(SHIFT_I_)          = 'SHIFT_I              '
    mon_str(NEUTRAL_)          = 'NEUTRAL              '
    mon_str(CHARGEI_)          = 'CHARGEI              '
    mon_str(POISSON_)          = 'POISSON              '
    mon_str(POISSON_E1_)       = 'POISSON_E1           '
    mon_str(PUSHE_)            = 'PUSHE                '
    mon_str(SHIFT_E_)          = 'SHIFT_E              '
    mon_str(CHARGEE_)          = 'CHARGEE              '
    mon_str(POISSON_E2_)       = 'POISSON_E2           '
    mon_str(F0_GRID_)          = 'F0_GRID              '
! Non-PETSc events
    mon_str(FINALIZE_)         = 'FINALIZE             '
    mon_str(PETSC_INIT_)       = 'PETSC_INIT           '
    mon_str(ADIOS_INIT_)       = 'ADIOS_INIT           '
    mon_str(SETUP_)            = 'SETUP                '
    mon_str(INIT_PUSHMOD_GPU_) = 'INIT_PUSHMOD_GPU     '
    mon_str(SET_WEIGHTS_F_)    = 'SET_WEIGHTS_F        '
    mon_str(SHIFT_I_F_)        = 'SHIFT_I_F            '
    mon_str(SHIFT_E_F_)        = 'SHIFT_E_F            '
    mon_str(POISSON_F_)        = 'POISSON_F            '
    mon_str(GET_POT_GRAD_F_)   = 'GET_POT_GRAD_F       '
    mon_str(MAIN_LOOP_RED_)    = 'MAIN_LOOP_RED        '
    mon_str(SET_WEIGHTS_)      = 'SET_WEIGHTS          '
    mon_str(IPC_LOOP_)         = 'IPC_LOOP             '
    mon_str(ELECTRONS_)        = 'ELECTRONS            '
    mon_str(GET_POT_GRAD_)     = 'GET_POT_GRAD         '
    mon_str(MEM_CLEAN_RSTRT_)  = 'MEM_CLEAN_RSTRT      '
    mon_str(MEM_CLEAN_I_)      = 'MEM_CLEAN_I          '
    mon_str(HYBRID_LOOP_)      = 'HYBRID_LOOP          '
    mon_str(DDPOTDT_UPD_)      = 'DDPOTDR_UPD          '
    mon_str(SOL_ELEC_PHASE_)   = 'SOL_ELEC_PHASE       '
    mon_str(GAT_FIELD_INFO_)   = 'GAT_FIELD_INFO       '

    mon_sync(:) = .false.
    mon_sync(INIT_)            = timing_barrier
    mon_sync(FIRST_)           = timing_barrier
    mon_sync(SETUP_)           = timing_barrier
    mon_sync(INIT_PUSHMOD_GPU_)= timing_barrier
    mon_sync(MAIN_LOOP_RED_)   = timing_barrier
    mon_sync(PUSH_I_)          = timing_barrier
    mon_sync(DIAGNOSIS_)       = timing_barrier
    mon_sync(SHIFT_I_)         = timing_barrier
    mon_sync(NEUTRAL_)         = timing_barrier
    mon_sync(CHARGEI_)         = timing_barrier
    mon_sync(POISSON_E1_)      = timing_barrier
    mon_sync(DDPOTDT_UPD_)     = timing_barrier
    mon_sync(SOL_ELEC_PHASE_)  = timing_barrier
    mon_sync(GAT_FIELD_INFO_)  = timing_barrier
    mon_sync(PUSHE_)           = timing_barrier
    mon_sync(SHIFT_E_)         = timing_barrier
    mon_sync(CHARGEE_)         = timing_barrier
    mon_sync(POISSON_E2_)      = timing_barrier
    mon_sync(F0_GRID_)         = timing_barrier
    mon_sync(POISSON_)         = timing_barrier
    mon_sync(GET_POT_GRAD_)    = timing_barrier
    mon_sync(COLLISION_)       = timing_barrier
    mon_sync(SOURCE_)          = timing_barrier
    mon_sync(RESTART_WRITE_)   = timing_barrier
    mon_sync(FINALIZE_)        = timing_barrier

    mon_adj(:) = 0
    mon_adj(TOTAL_)            = 1
    mon_adj(SHIFT_I_F_)        = 1
    mon_adj(SHIFT_E_F_)        = 1
    mon_adj(POISSON_F_)        = 1
    mon_adj(GET_POT_GRAD_F_)   = 1

    do i=1,mon_NX
       ctemp1 = mon_str(i)
       ctemp2(1: 5) = 'sync_'
       ctemp2(6:21) =  ctemp1(1:16)
       mon_str(i+mon_NX) = ctemp2
    enddo

#if defined(NO_PETSC)
    mon_time(:)=0D0
    mon_sum(:)=0D0
#endif

  end subroutine init_perf_monitor

  subroutine petsc_perf_init( ierror )
  integer ierror
  !
  integer i, ierr
  do i=1,mon_N
     call PetscLogEventRegister( mon_str(i), 0, event(i), ierr )
  enddo
  ierror=0
  end subroutine petsc_perf_init

  subroutine flush_perf_monitor( istep )
#if defined(CAM_TIMERS)
    use sml_module
    use rem_module
    use perf_mod, only: t_prf, t_startf, t_stopf
#endif
    implicit none
    integer istep
#if defined (CAM_TIMERS)
    integer ierr
    integer i, str_length
    logical exist
    character (len=40):: filename
    character (len=9):: cstep
    include 'mpif.h'

    call t_startf("sync1_t_prf")
    call mpi_barrier(sml_comm,ierr)
    call t_stopf("sync1_t_prf")

    do i=1,40
      filename(i:i) = " "
    enddo

    str_length = 0
    if (sml_mype == 0) then
       inquire(DIR_INQ="timing",EXIST=exist)
       if (exist) then
          filename(1:7) = "timing/"
          str_length = 7
          inquire(DIR_INQ="timing/checkpoints",EXIST=exist)
          if (exist) then
             filename(str_length+1:str_length+12) = "checkpoints/"
             str_length = str_length + 12
          endif
       endif
    endif

    if (istep < 0) then
      filename(str_length+1:str_length+20) = "timing_initializ.txt"
    else
      write(cstep,'(i9.9)') istep
      filename(str_length+1:str_length+7) = "timing_"
      filename(str_length+8:str_length+16) = cstep(1:9)
      filename(str_length+17:str_length+20) = ".txt"
    endif
    call t_prf(filename=trim(filename), mpicom=sml_comm, &
               num_outpe=1, global_stats=.true.)

!   'hijack' flush_perf_monitor to also update remaining time estimation
    if ((rem_restart_any) .and. (rem_estimate)) then
       call mpi_bcast(rem_final_step_est, 1, mpi_integer, 0, sml_comm, ierr)
       rem_final_step_upd = .true.
    endif
!
    call t_startf("sync2_t_prf")
    call mpi_barrier(sml_comm,ierr)
    call t_stopf("sync2_t_prf")
#endif

  end subroutine flush_perf_monitor

  subroutine finish_perf_monitor( ierr )
    use sml_module
#if defined(CAM_TIMERS)
    use perf_mod, only: t_finalizef, t_prf
#endif
    implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h>
#else
#include <petsc/finclude/petsc.h>
#endif
    PetscErrorCode::i,ierr
    integer str_length
    logical exist
    character (len=40):: filename
    call petsc_end(ierr);CHKERRQ(ierr)

#if defined (CAM_TIMERS)
    do i=1,40
      filename(i:i) = " "
    enddo

    str_length = 0
    exist = .false.
    if (sml_mype == 0) then
       inquire(DIR_INQ="timing",EXIST=exist)
    endif
    call mpi_bcast(exist,1,mpi_logical,0,sml_comm,ierr)
    if (exist) then
      filename(1:7) = "timing/"
      str_length = 7
    endif

    filename(str_length+1:str_length+14) = "timing_all.txt"
    call t_prf(filename=trim(filename), mpicom=sml_comm)
    call t_finalizef ()
#endif

  end subroutine finish_perf_monitor

#if !defined(CAM_TIMERS)
  subroutine t_startf( eventname )
  character(len=*), intent(in) :: eventname
  return
  end subroutine t_startf

  subroutine t_stopf( eventname )
  character(len=*), intent(in) :: eventname
  return
  end subroutine t_stopf

  subroutine t_adj_detailf( adjustment )
  integer, intent(in) :: adjustment
  return
  end subroutine t_adj_detailf
#endif

end module perf_monitor

#ifndef OPTIM_GYRO_AVG_MAT
module mat_class
  type mat_type
     integer :: n,m,width
     real (8), allocatable :: value(:,:)
     integer, allocatable :: eindex(:,:),nelement(:)
  end type mat_type
contains
  ! create a empty matrix
  subroutine new_mat(mat,n,w)
    use sml_module, only : sml_mype
    implicit none
    type(mat_type) :: mat
    integer,intent(in) :: n,w
    integer :: st
    mat%n=n
    mat%m=n ! default, square matrix
    mat%width=w
    allocate(mat%value(w,n), mat%eindex(w,n), mat%nelement(n),STAT=st)
    if(st>0) then
       print *, 'insufficient memory for matrix allocation', n,w
       stop
    endif
    mat%nelement=0
    mat%value=0D0
  end subroutine new_mat

  subroutine set_non_square_mat(mat,m)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: m
    mat%m = m
  end subroutine set_non_square_mat

  !delete a matrix
  subroutine del_mat(mat)
    implicit none
    type(mat_type) :: mat

    if(allocated(mat%value)) then
       deallocate(mat%value,mat%eindex,mat%nelement)
    endif
    mat%n=-1
    mat%width=-1
  end subroutine del_mat

    ! set value
  subroutine set_value(mat,i,j,value,flag)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: i,j,flag
    real (8) :: value
    integer :: l

    !check same node
    do l=1,mat%nelement(i)
       ! redundant point
       if(j==mat%eindex(l,i))then
          if(flag==0) then
             mat%value(l,i)=value
          else
             mat%value(l,i)=mat%value(l,i)+value
          endif
          return
       endif
    enddo

    ! new point
    if(mat%width <= mat%nelement(i) ) then
       print *, 'Error : not enough memory space for matrix.'
       print *, '        Increase width ', mat%width, mat%nelement(i),i
       stop
       return
    endif
    mat%nelement(i)=mat%nelement(i)+1
    l=mat%nelement(i)
    mat%eindex(l,i)=j
    mat%value(l,i)=value
  end subroutine set_value
  subroutine get_max_width(mat,mwidth)
    implicit none
    type(mat_type) :: mat
    integer, intent(out) :: mwidth
    integer :: i

    mwidth=-1
    do i=1, mat%n
       if(mat%nelement(i)> mwidth) then
          mwidth=mat%nelement(i)
       endif
    enddo
  end subroutine get_max_width
  subroutine output_matrix(mat,filename)
    implicit none
    integer :: i,j
    type(mat_type) :: mat
    character (len=20) :: filename
    open(unit=201,file=filename,status='replace')
    do i=1, mat%n
     do j=1, mat%nelement(i)
        write(201,*) i,mat%eindex(j,i),mat%value(j,i)
     enddo
    enddo
    close(201)

  end subroutine output_matrix

  subroutine mat_mult(mat, x, y)
    implicit none
    type(mat_type) :: mat
    real (8), intent(in) :: x(mat%m)
    real (8), intent(out) :: y(mat%n)
    integer :: i, j

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       y(i)=0D0
       do j=1, mat%nelement(i)
          y(i) = y(i) + mat%value(j,i) * x(mat%eindex(j,i))
       enddo
    enddo
  end subroutine mat_mult
  ! matrix multiplication of transpose of M : y= Mat^T x
  subroutine mat_transpose_mult(mat,x,y)
    implicit none
    type(mat_type) :: mat
    real (8), intent(in) :: x(mat%n)
    real (8), intent(out) :: y(mat%m)
    integer :: i, j, k

    y=0D0

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       do j=1, mat%nelement(i)
          k = mat%eindex(j,i)
          y(k) = y(k) + mat%value(j,i) * x(i)
       enddo
    enddo
  end subroutine mat_transpose_mult

  subroutine mat_mult_tensor(mat, x, nv, y)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: nv
    real (8), intent(in) :: x(nv,mat%n)
    real (8), intent(out) :: y(nv,mat%m)
    integer :: i, j

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       y(:,i)=0D0
       do j=1, mat%nelement(i)
          y(:,i) = y(:,i) + mat%value(j,i) * x(:,mat%eindex(j,i))
       enddo
    enddo
  end subroutine mat_mult_tensor


end module mat_class
#else
! ======================================
! New non-regular matrices
! @author: J. Dominski
module mat_class
  type mat_type
     integer :: n,m,ntot  ! n == #nodes 
     real (8), allocatable :: value(:)
     integer, allocatable :: width(:),nelement(:),eindex(:),a(:) !width is not absolutely necessary, can implement other algo
  end type mat_type
contains
  ! create a empty matrix
  subroutine new_mat(mat,n,w)
    !use sml_module, only : sml_mype
    implicit none
    type(mat_type) :: mat
    integer,intent(in) :: n,w
    integer :: st
    mat%n=n
    mat%m=n
    allocate(mat%width(n),mat%nelement(n),mat%a(n),STAT=st)
    if(st>0) then
       print *, 'insufficient memory for matrix allocation', n
       stop
    endif
    mat%nelement=0
    mat%width=w
    if(w.gt.0)then
      call new_value(mat)
    endif
  end subroutine new_mat

  ! allocate values 
  subroutine new_value(mat)
    implicit none
    type(mat_type) :: mat
    integer :: st,i
    if(allocated(mat%value)) deallocate(mat%value,mat%eindex)
    mat%ntot=sum(mat%width)
    !mat%ntot=sum(mat%nelement)
    !mat%nelement=0
    mat%a(1)=0
    do i=2,mat%n
      mat%a(i)=sum(mat%width(1:(i-1)))
    enddo
    allocate(mat%value(mat%ntot),mat%eindex(mat%ntot),STAT=st)
    if(st>0) then
       print *, 'insufficient memory for matrix allocation', mat%ntot
       stop
    endif
    mat%value=0D0
  end subroutine new_value

  subroutine set_non_square_mat(mat,m)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: m
    mat%m = m
  end subroutine set_non_square_mat

  !delete a matrix
  subroutine del_mat(mat)
    implicit none
    type(mat_type) :: mat

    if(allocated(mat%value)) then
       deallocate(mat%value,mat%eindex)
    endif

    if(allocated(mat%nelement))then
       deallocate(mat%nelement,mat%width)
    endif

  end subroutine del_mat
!  function get_value(mat,i,j)
!    implicit none
!    type(mat_type) :: mat
!    integer, intent(in) :: i,j
!    real (8) :: get_value
!    integer :: l,a
!    a=i*(mat%n-1)
!    get_value=0D0
!    !check same node
!    do l=1,mat%nelement(i)
!       ! redundant point
!       if(j==mat%eindex(a+l))then
!         get_value=mat%value(a+l)
!         return
!       endif
!    enddo
!  end subroutine
    ! set value
  subroutine set_value(mat,i,j,value,flag)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: i,j,flag
    real (8) :: value
    integer :: l,a
    a=mat%a(i)
    !check same node
    do l=1,mat%nelement(i)
       ! redundant point
       if(j==mat%eindex(a+l))then
          if(flag==0) then
             mat%value(a+l)=value
          else
             mat%value(a+l)=mat%value(a+l)+value
          endif
          return
       endif
    enddo

    ! new point
    if(mat%width(i) <= mat%nelement(i) ) then
       print *, 'Error : not enough memory space for matrix.'
       print *, '        Increase width ', mat%width(i), mat%nelement(i),i,j
       stop
       return
    endif
    mat%nelement(i)=mat%nelement(i)+1
    l=mat%nelement(i)
    mat%eindex(a+l)=j
    mat%value(a+l)=value
  end subroutine set_value
  subroutine set_value_width(mat,i,j)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: i,j
    integer :: l
    do l=1,mat%width(i)
       if(j==mat%eindex(mat%a(i)+l)) return
    enddo
    mat%width(i)=mat%width(i)+1
    l=mat%nelement(i)
  end subroutine set_value_width
  subroutine get_max_width(mat,mwidth)
    implicit none
    type(mat_type) :: mat
    integer, intent(out) :: mwidth
    integer :: i

    mwidth=-1
    do i=1, mat%n
       if(mat%nelement(i)> mwidth) then
          mwidth=mat%nelement(i)
       endif
    enddo
  end subroutine get_max_width
  subroutine output_matrix(mat,filename)
    implicit none
    integer :: i,j
    type(mat_type) :: mat
    character (len=20) :: filename
    open(unit=201,file=filename,status='replace')
    do i=1, mat%n
     do j=1, mat%nelement(i)
        write(201,*) i,mat%eindex(mat%a(i)+j),mat%value(mat%a(i)+j)
     enddo
    enddo
    close(201)

  end subroutine output_matrix

  subroutine mat_mult(mat, x, y)
    implicit none
    type(mat_type) :: mat
    real (8), intent(in) :: x(mat%m)
    real (8), intent(out) :: y(mat%n)
    integer :: i, j

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       y(i)=0D0
       do j=1, mat%nelement(i)
          y(i) = y(i) + mat%value(mat%a(i)+j) * x(mat%eindex(mat%a(i)+j))
       enddo
    enddo
  end subroutine mat_mult
  ! matrix multiplication of transpose of M : y= Mat^T x
  subroutine mat_transpose_mult(mat,x,y)
    implicit none
    type(mat_type) :: mat
    real (8), intent(in) :: x(mat%n)
    real (8), intent(out) :: y(mat%m)
    integer :: i, j, k

    y=0D0

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       do j=1, mat%nelement(i)
          k = mat%eindex(mat%a(i)+j)
          y(k) = y(k) + mat%value(mat%a(i)+j) * x(i)
       enddo
    enddo
  end subroutine mat_transpose_mult

  subroutine mat_mult_tensor(mat, x, nv, y)
    implicit none
    type(mat_type) :: mat
    integer, intent(in) :: nv
    real (8), intent(in) :: x(nv,mat%n)
    real (8), intent(out) :: y(nv,mat%m)
    integer :: i, j

!!$###OMP PARALLEL DO PRIVATE(I,J)
    do i=1, mat%n
       y(:,i)=0D0
       do j=1, mat%nelement(i)
          y(:,i) = y(:,i) + mat%value(mat%a(i)+j) * x(:,mat%eindex(mat%a(i)+j))
       enddo
    enddo
  end subroutine mat_mult_tensor


end module mat_class
#endif

module boundary_class
  type range_type
     integer :: start,end
  end type range_type
  type boundary2_type
     integer :: nseg
     type(range_type), allocatable :: iseg(:)
  end type boundary2_type
contains
  logical function is_inside(i,bd)
    ! check if inside of boundary
    implicit none
    type(boundary2_type) :: bd
    integer , intent(in) ::i
    integer :: is

    is_inside=.true.  ! - the other way would be faster -- # of inside nodes > # of boundary nodes

    do is=1, bd%nseg
       if( bd%iseg(is)%start<= i .and. i<= bd%iseg(is)%end ) then
          is_inside=.false.
          return
       endif
    enddo
  end function is_inside

  subroutine set_boundary2_values(arr,value,bd)
    implicit none
    type(boundary2_type) :: bd
    real (8) :: arr(bd%iseg(bd%nseg)%end), value
    integer :: i

    do i=1, bd%nseg
       arr(bd%iseg(i)%start:bd%iseg(i)%end)=value
    enddo

  end subroutine set_boundary2_values

end module boundary_class


module smooth_module
  use mat_class
  type smooth_type
     real (8), allocatable :: weight(:)
     integer :: n, mode, type
  end type smooth_type
  type smooth_r_type
     type(mat_type) :: mat
     real(8) :: d0 ! basic distance
     integer :: n ! number of point
     integer :: type
  end type smooth_r_type
  type smooth_nearx_type
     type(mat_type) :: mat
     real(8) :: d0 ! distance
     integer :: n
     integer :: nsmth
     integer ,allocatable :: nodes(:)
     real (8),allocatable :: tmp(:)
     real (8) :: dr, dz

  end type smooth_nearx_type
  type(smooth_type) :: smooth00,smooth0L,smoothH,smoothdiag
  type(smooth_r_type) :: smooth_r1
  type(smooth_nearx_type) :: smooth_nearx1
  integer :: smooth_n_in, smooth_mode_in,smooth_type_in
  integer :: smooth_H_n_in, smooth_H_mode_in, smooth_H_type_in
  integer :: smooth_diag_n_in, smooth_diag_mode_in, smooth_diag_type_in
  ! input parameters for smooth_r
  real (8) :: smooth_r1_d0_in
  integer :: smooth_r1_n_in, smooth_r1_type_in
  ! input parameters for smooth_nearx
  real (8) :: smooth_nearx_d0,smooth_nearx_dr,smooth_nearx_dz
  integer :: smooth_nearx_nsmth

  ! new poloidal smoothing operator
  type(mat_type) :: smooth_pol_mat

  ! Smoothing of poloidal potential before computing the poloidal field
  logical :: smooth_pol_efield, smooth_rad_efield
  integer :: smooth_pol_width
  real (kind=8) :: smooth_pol_d0
  real (kind=8) :: smooth_grid_rad_width, smooth_grid_rad_sigma


  contains
    ! initialize smooth_init
    subroutine smooth_init(smooth)
      implicit none
      type(smooth_type) :: smooth
      integer :: i

      if(smooth%mode==-1) return  ! No smooth
      if(smooth%mode== 0) return  ! flux average

      allocate(smooth%weight(smooth%n))
      do i=1, smooth%n
         if(smooth%type==1) then
            smooth%weight(i)=exp(- real(i-1)**2/real(smooth%n)**2*4D0 )
         else
            smooth%weight(i)=1D0
         endif
      enddo
      smooth%weight(1)=0.5D0*smooth%weight(1)
    end subroutine smooth_init
    !
    ! delete smooth_pol object
    subroutine smooth_pol_delete(smooth)
      implicit none
      type(smooth_type) :: smooth
      if(allocated(smooth%weight)) then
         deallocate(smooth%weight)
      endif
    end subroutine smooth_pol_delete
    !
    ! initialize smooth_r
    subroutine smooth_r_init1(smooth_r,d0,n,type,nnode)
      implicit none
      type(smooth_r_type) :: smooth_r
      real (8) :: d0
      integer :: n, type,nnode

      smooth_r%n=n
      smooth_r%d0=d0
      smooth_r%type=type

      call new_mat(smooth_r%mat,nnode,(2*n+1)*3)

      ! call smooth_r_init2 with grid information

    end subroutine smooth_r_init1
    !
    ! delete smooth_r object
    subroutine smooth_r_delete(smooth_r)
      implicit none
      type(smooth_r_type) :: smooth_r
!      if(allocated(smooth_r%mat)) then
         call del_mat(smooth_r%mat)
!      endif
    end subroutine smooth_r_delete
end module

!! OpenMP module
module omp_module
contains
  subroutine split_indices(total,num_pieces,ibeg,iend)
    implicit none

    integer :: total
    integer :: num_pieces
    integer :: ibeg(num_pieces), iend(num_pieces)
    integer :: itmp1, itmp2, ioffset, i

    if (num_pieces > 0) then
       itmp1 = total/num_pieces
       itmp2 = mod(total,num_pieces)
       ioffset = 0
       do i=1,itmp2
          ibeg(i) = ioffset + 1
          iend(i) = ioffset + (itmp1+1)
          ioffset = iend(i)
       enddo
       do i=itmp2+1,num_pieces
          ibeg(i) = ioffset + 1
          if (ibeg(i) > total) then
             iend(i) = ibeg(i) - 1
          else
             iend(i) = ioffset + itmp1
             ioffset = iend(i)
          endif
       enddo
    endif

  end subroutine split_indices
end module omp_module

module src_module
  ! module
  integer, parameter :: src_nmax=4
  integer :: src_narea, src_narea_e, src_niter
  real (8), dimension(src_nmax) :: src_pin, src_pout, src_decay_width, src_heat_power, src_torque
  real (8), dimension(src_nmax) :: src_pin_e, src_pout_e, src_decay_width_e, src_heat_power_e, src_torque_e
  integer :: src_period, src_nsubsection
  real (8) :: src_pin1, src_pout1, src_decay_width1, src_heat_power1, src_torque1
  real (8) :: src_pin2, src_pout2, src_decay_width2, src_heat_power2, src_torque2
  real (8) :: src_pin3, src_pout3, src_decay_width3, src_heat_power3, src_torque3
  real (8) :: src_pin4, src_pout4, src_decay_width4, src_heat_power4, src_torque4

  real (8) :: src_pin1_e, src_pout1_e, src_decay_width1_e, src_heat_power1_e, src_torque1_e
  real (8) :: src_pin2_e, src_pout2_e, src_decay_width2_e, src_heat_power2_e, src_torque2_e
  real (8) :: src_pin3_e, src_pout3_e, src_decay_width3_e, src_heat_power3_e, src_torque3_e
  real (8) :: src_pin4_e, src_pout4_e, src_decay_width4_e, src_heat_power4_e, src_torque4_e
contains
  subroutine src_setup
    src_pin         = (/ src_pin1,         src_pin2,        src_pin3,        src_pin4   /)
    src_pout        = (/ src_pout1,        src_pout2,       src_pout3,       src_pout4  /)
    src_decay_width = (/ src_decay_width1, src_decay_width2,src_decay_width3,src_decay_width4 /)
    src_heat_power  = (/ src_heat_power1,  src_heat_power2, src_heat_power3, src_heat_power4 /)
    src_torque      = (/ src_torque1,      src_torque2,     src_torque3,     src_torque4/)

    src_pin_e         = (/ src_pin1_e,         src_pin2_e,        src_pin3_e,        src_pin4_e   /)
    src_pout_e        = (/ src_pout1_e,        src_pout2_e,       src_pout3_e,       src_pout4_e  /)
    src_decay_width_e = (/ src_decay_width1_e, src_decay_width2_e,src_decay_width3_e,src_decay_width4_e /)
    src_heat_power_e  = (/ src_heat_power1_e,  src_heat_power2_e, src_heat_power3_e, src_heat_power4_e /)
    src_torque_e      = (/ src_torque1_e,      src_torque2_e,     src_torque3_e,     src_torque4_e/)

  end subroutine src_setup
end module src_module


module rad_module
  ! impurity radiation module
  character (len=80):: rad_filename
  
  integer :: rad_nT, rad_nn0ne

  real (8), allocatable :: rad_Te_ev(:), rad_n0_ne(:)
  real (8), allocatable :: rad_Lz(:,:), rad_avgZ(:,:), rad_avgZ2(:,:)

  real (8) :: rad_impurity_fraction ! simplest model : n_imp / n_e is constant
  
  ! for impurity fraction using z_eff
  logical :: rad_use_zeff_profile
  logical :: rad_use_fix_charge
  real (8) :: rad_z  ! fixed Z for impurity - Z_Carbon --> 6, Z^2 --> 36

end module rad_module


MODULE Ecuyer_random
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

  ! These are unsigned integers in the C version
  INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890

CONTAINS

  SUBROUTINE init_seeds(i1, i2, i3)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i1, i2, i3

    s1 = i1
    s2 = i2
    s3 = i3
    IF (IAND(s1,-2) == 0) s1 = i1 - 1023
    IF (IAND(s2,-8) == 0) s2 = i2 - 1023
    IF (IAND(s3,-16) == 0) s3 = i3 - 1023

    RETURN
  END SUBROUTINE init_seeds

  FUNCTION taus88() RESULT(random_numb)
    ! Generates a random number between 0 and 1.  Translated from C function in:
    ! Reference:
    ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
    ! generators', Math. of Comput., 65, 203-213.

    ! The cycle length is claimed to be about 2^(88) or about 3E+26.
    ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

    IMPLICIT NONE
    REAL (dp) :: random_numb

    INTEGER   :: b

    ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
    !      to the left if j > 0, otherwise to the right.

    b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
    s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
    b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
    s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
    b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
    s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
    random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dp + 0.5_dp

    RETURN
  END FUNCTION taus88

  SUBROUTINE init_seeds_ext(sv)
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: sv(3)

    IF (IAND(sv(1),-2) == 0)  sv(1) = sv(1) - 1023
    IF (IAND(sv(2),-8) == 0)  sv(2) = sv(2) - 1023
    IF (IAND(sv(3),-16) == 0) sv(3) = sv(3) - 1023

    RETURN
  END SUBROUTINE init_seeds_ext

  FUNCTION taus88_ext(sv) RESULT(random_numb)
    ! Generates a random number between 0 and 1.  Translated from C function in:
    ! Reference:
    ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
    ! generators', Math. of Comput., 65, 203-213.

    ! The cycle length is claimed to be about 2^(88) or about 3E+26.
    ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: sv(3)
    REAL (dp) :: random_numb

    INTEGER   :: b, i1, i2, i3

    ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
    !      to the left if j > 0, otherwise to the right.

    i1 = sv(1)
    i2 = sv(2)
    i3 = sv(3)
    b  = ISHFT( IEOR( ISHFT(i1,13), i1), -19)
    i1 = IEOR( ISHFT( IAND(i1,-2), 12), b)
    b  = ISHFT( IEOR( ISHFT(i2,2), i2), -25)
    i2 = IEOR( ISHFT( IAND(i2,-8), 4), b)
    b  = ISHFT( IEOR( ISHFT(i3,3), i3), -11)
    i3 = IEOR( ISHFT( IAND(i3,-16), 17), b)
    random_numb = IEOR( IEOR(i1,i2), i3) * 2.3283064365E-10_dp + 0.5_dp
    sv(1) = i1
    sv(2) = i2
    sv(3) = i3

    RETURN
  END FUNCTION taus88_ext

END MODULE Ecuyer_random

!! random number generator module
module random_xgc
  use Ecuyer_random, only: init_seeds, init_seeds_ext, taus88, taus88_ext
  use sml_module
  implicit none
  type seeds_type
     integer, pointer :: s(:)
  end type seeds_type
  type(seeds_type), allocatable :: sv(:)

contains
  !random number - seed initialize
  subroutine init_ranx()
    implicit none
#ifdef _OPENMP
    integer :: i, gi
    allocate( sv(0:sml_nthreads-1) )
    !$OMP PARALLEL DO PRIVATE( I, GI )
    do i=0,sml_nthreads-1
       allocate(sv(i)%s(3))
       gi = sml_mype*sml_nthreads + i
       sv(i)%s(1) = 1234*gi
       sv(i)%s(2) = 2345*gi + 6789
       sv(i)%s(3) = 4321*gi + 10
       call init_seeds_ext(sv(i)%s)
    end do
#else
    call init_seeds(1234*sml_mype,2345*sml_mype+6789,4321*sml_mype+10)
#endif
  end subroutine init_ranx

  ! random number generator
  function ranx()
    implicit none
    real (8) :: ranx
#ifdef _OPENMP
    integer :: thread_id
    integer :: omp_get_thread_num
    thread_id = omp_get_thread_num()
    ranx=taus88_ext( sv(thread_id)%s )
#else
    ranx=taus88()
#endif
  end function ranx

end module random_xgc


module lim_module
!  character (LEN=64) :: lim_filename  !> limiter file name , need clean up lim_module
  real (kind=8), allocatable :: lim_r(:,:), lim_z(:,:), lim_weight(:,:), lim_en(:,:), &
      lim_org_r(:), lim_org_z(:), lim_ds(:,:), lim_psi(:,:)
  real (kind=8) :: lim_dz, lim_zmin, lim_zmax, lim_r0_down, lim_r0_up,  lim_psi_min
  integer :: lim_zindex(2), lim_store_mz, lim_mdata, lim_store_mz_rl(2)
end module


module neu_module
#ifdef DEGAS2
  use grid_class
  use lim_module
#endif
  character (LEN=128) :: neu_sepfile, neu_limfile
  integer :: neu_col_mode, neu_adjust_n0, neu_cx_period, &
             neu_ion_period, neu_varying_mfp
  integer :: neu_col_period, neu_elastic_col_on

  logical :: neu_update_elec

  real (kind=8) :: neu_t0

  integer, parameter :: neu_sep_mtheta=200
  integer, parameter :: diag_flow_npsi=80 ! this name will be changed.
  integer :: neu_sep_mtheta_file
  real (kind=8) :: neu_sep_r(neu_sep_mtheta), neu_sep_z(neu_sep_mtheta)
  real (kind=8), allocatable :: neu_sep_r_file(:), neu_sep_z_file(:)
  real (kind=8) :: neu_n0, neu_delta_n, neu_delta_theta, neu_mfp, &
                   neu_theta_x, neu_theta_x2, neu_recycle_rate, neu_mfp0, neu_base_den, neu_sep_dtheta
  real (8), allocatable :: neu_weight_sum_lost(:)
  real (kind=8) :: neu_old_inside_weight,   &
                   neu_mom_sum_lost(3), neu_energy_sum_lost,      &
                   neu_actual_ionized, neu_actual_accum_ionized,  &
                   neu_ionize2_psi, neu_weight_accum_lost, neu_weight_accum_ionized, &
                   neu_peak_theta, neu_delta_theta_lowerx1, neu_delta_theta_lowerx2, &
                   neu_delta_theta_upperx1, neu_delta_theta_upperx2
  real (kind=8) :: neu_theta_edge(2), neu_rad_edge(2), neu_theta_edge2(2), neu_rad_edge2(2)
  integer :: neu_ionize_mode, neu_start_time
  real (kind=8) :: neu_temp_factor

  !---------- neutral2 ----------------
  integer, parameter :: neu_grid_mpsi=50, neu_grid_mtheta=40
  integer, parameter :: neu_mtheta=200
  real (kind=8) :: neu_r0(neu_mtheta), neu_z0(neu_mtheta), &
       neu_jacob(neu_mtheta), neu_dtheta, neu_psi_edge, neu_temp0
  real (kind=8) :: neu_dt
  integer :: neu_istep_max, neu_num, neu_mpol, neu_monte_num, neu_flow_on, neu_grid_mode
  real (kind=8) :: neu_grid_den(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_temp(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_wsum(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_esum(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_vol(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_flow(neu_grid_mpsi,neu_grid_mtheta), &
       neu_grid_vsum(neu_grid_mpsi,neu_grid_mtheta)
  real (kind=8) :: neu_grid_dtheta, neu_grid_dpsi, neu_grid_dpsi2, &
                   neu_grid_min_psi, neu_grid_max_psi
  integer :: neu_mode2_period
  real (kind=8) :: neu_ionize_elec_cooling_ev
  real (kind=8) :: neu_nr ! average neutral density at psi=neu_grid_max_psi normalized to neu_n0

  integer (8) :: neu_cx_events(neu_grid_mpsi), neu_ion_events(neu_grid_mpsi)
  real (kind=8) :: neu_ion_events_weights(neu_grid_mpsi), neu_cx_events_weights(neu_grid_mpsi)
  real (kind=8) :: neu_lim_r(neu_sep_mtheta), neu_lim_z(neu_sep_mtheta)

  real (8) :: neu_weight_lost_max ! # of particles
  real (8) :: neu_lost_rate_max ! # of particle per sec


  logical :: neu_enforce_no_neutral ! for radiation cooling
end module neu_module


module input_module
#ifndef USE_OLD_READ_INPUT
!-----------------------------------------------------------------------
!
! Purpose: This module is responsible for reading in input namelist
!          file as a string, broadcasting to all processes, and then
!          finding start of a particular namelist group.
!
! Author:  P. Worley and E. D'Azevedo, January 2013
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
  implicit none
  private                   ! Make the default access private
  save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
  public read_input_file
  public read_namelist_index

!-----------------------------------------------------------------------
! Public data ---------------------------------------------------------
!-----------------------------------------------------------------------
  character (len=:), allocatable, public :: input_string

contains

  subroutine read_input_file
     use sml_module, only : sml_comm, sml_mype
     implicit none
     include 'mpif.h'

     logical :: comment, first
     integer :: input_lth, ierror, i
     character (len=256) :: line

     if (sml_mype == 0) then
        open(unit=14,file='input',access="sequential",form="formatted",action='read')

        ! read in line-by-line and determine required size of input string,
        ! then allocate string
        rewind(14)
        input_lth = 0
        ierror = 0
        do while (ierror .eq. 0)
           line = ' '
           comment = .false.
           read(14,'(A)',iostat=ierror) line
           if (ierror == 0) then
              do i = 1, len_trim(line)
                 if (line(i:i) .eq. '!') comment = .true.
                 if (comment) line(i:i) = ' '
              end do
              input_lth = input_lth + len_trim(line) + 3
           endif
        enddo

        ! allocate string
        allocate(character(len=input_lth) :: input_string)

        ! read in line-by-line and concatenate it to input string
        rewind(14)
        input_string = ""
        first = .true.
        ierror = 0
        do while (ierror .eq. 0)
           line = ' '
           comment = .false.
           read(14,'(A)',iostat=ierror) line
           if (ierror == 0) then
              do i = 1, len_trim(line)
                 if (line(i:i) .eq. '!') comment = .true.
                 if (comment) line(i:i) = ' '
              end do
              if (first) then
                 input_string = trim(line)
                 first = .false.
              else
                 input_string = trim(input_string) // "  " // trim(line)
              endif
           endif
        enddo

        close (14)

        ! set length to actual length and broadcast to other processes
        input_lth = len(input_string)
        call mpi_bcast( input_lth, 1, MPI_INTEGER, 0, sml_comm, ierror)

     else

        ! get required size of input character array,then allocate it
        call mpi_bcast( input_lth, 1, MPI_INTEGER, 0, sml_comm, ierror)
        allocate(character(len=input_lth) :: input_string)

     endif

     ! broadcast input character array
     call mpi_bcast( input_string, input_lth, MPI_CHARACTER, 0, sml_comm, ierror)

   end subroutine read_input_file

   integer function read_namelist_index(input_string,namelist_string)
     use sml_module, only : sml_comm, sml_mype
     implicit none
     include 'mpif.h'

     character (len=*), intent(in) :: input_string
     character (len=*), intent(in) :: namelist_string

     integer :: input_len, nml_len_m1, i, ierror

     input_len = len_trim(input_string)
     nml_len_m1 = len_trim(namelist_string) - 1

     i=1
     do while ((input_string(i:i+nml_len_m1) .ne. namelist_string) &
        .and. (i+nml_len_m1 < len_trim(input_string)))

        do while ((input_string(i:i) .ne. "&") &
           .and. (i+nml_len_m1 < len_trim(input_string)))
           i = i + 1
        enddo

        i = i + 1

     enddo

     if (input_string(i:i+nml_len_m1) .ne. namelist_string) then
        if (sml_mype == 0) then
           write(6,*) &
"namelist read/internal file/end of file reached without finding group:", namelist_string
!pw        call flush(6)
        endif
        call mpi_abort(sml_comm, 1, ierror)
     endif

     read_namelist_index = i-1

   end function read_namelist_index
#endif
end module input_module

#include "assert_mod.F90"
#include "comm_mod.F90"

