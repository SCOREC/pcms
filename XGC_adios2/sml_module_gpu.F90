module sml_module_gpu
use precision_mod_gpu
use sml_module, only :  &
   sml_e_charge_host => sml_e_charge, &
   sml_epsilon0_host => sml_epsilon0, &
   sml_prot_mass_host => sml_prot_mass, &
   sml_elec_mass_host => sml_elec_mass, &
   sml_n_vf_diag_host => sml_n_vf_diag, &
   sml_nlarmor_host => sml_nlarmor, &
   sml_nrk_host => sml_nrk, &
   sml_boundary_diagonal_host => sml_boundary_diagonal, &
   sml_j2ev_host => sml_j2ev,           &
   sml_ev2j_host => sml_ev2j

implicit none
  !! delta-f weight evolution switch. 
  !! false for whole (grad_b, curl_b,and ExB), true for exb only
  logical :: sml_dwdt_exb_only 

  !! delta-f weight evolution switch. 
  !! false for (1-w) factor true for (1) factor. default is .false.
  logical :: sml_dwdt_fix_bg  

  logical :: sml_turb_efield    ! set zero turbulece efield
  logical :: sml_00_efield      ! set zero 00 efield

  !! delta-f switch. 0 for off, 1 for on. delta-f simulation 
  !! is not verified yet, especially for output file.
  logical :: sml_deltaf
  logical :: sml_electron_on 
  logical :: sml_extra_dwdt  ! extra dwdt 
  integer :: sml_deltaf_f0_mode

  !! Bounce routine switch  0 for off, 
  !!                        1 for inner boundary, 
  !!                        2 for both boundaries
  integer :: sml_bounce         
  integer :: sml_bounce_zero_weight

  !! Neutral routine switch  .f. for off
  logical :: sml_neutral

  !!  simulation boundary in R-Z space
  real (kind=work_p) :: sml_2pi, sml_2pi_wedge_n
  real (kind=work_p) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z 



  !! for sml_deltaf_f0_mode==-1  1/Ln , 1/Lt - artificial f0
  real (kind=work_p) :: sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e  



  real(kind=work_p), parameter :: sml_e_charge = sml_e_charge_host
  real(kind=work_p), parameter :: sml_epsilon0 = sml_epsilon0_host
  real(kind=work_p), parameter :: sml_prot_mass = sml_prot_mass_host
  real(kind=work_p), parameter :: sml_elec_mass = sml_elec_mass_host

  real(kind=work_p), parameter :: sml_j2ev = sml_j2ev_host
  real(kind=work_p), parameter :: sml_ev2j = sml_ev2j_host

  real(kind=work_p), parameter :: sml_boundary_diagonal = sml_boundary_diagonal_host
  integer, parameter :: sml_n_vf_diag = sml_n_vf_diag_host
  integer, parameter :: sml_nlarmor = sml_nlarmor_host
  integer, parameter :: sml_nrk = sml_nrk_host

  !! Inner and outer boundary for initial loading
  real(kind=work_p) :: sml_inpsi, sml_outpsi 

  integer :: sml_comm !! communicator
  integer :: sml_comm_null !! MPI_COMM_NULL for adios posix
  integer :: sml_totalpe !! total number of processors
  integer :: sml_mype !! process index
  integer :: sml_nthreads
  integer :: sml_gstep
  integer :: sml_ipc, ipc_gpu
  integer :: sml_sheath_mode
 
  real (kind=work_p) :: sml_time !! simulation time 
  !! -1 for inverted B   1 for normal B <-- given by sml_invert_B
  real (kind=work_p) :: sml_bp_sign       
  !! -bp_sign for co-current bt, 
  !!  bp_sign for couter-current bt <- given by sml_co_curr_bt
  real (kind=work_p) :: sml_bt_sign       

  logical :: sml_ignore_drift_near_wall
  real (kind=work_p):: sml_ignore_drift_r0
  real (kind=work_p):: sml_ignore_drift_z0
  real (kind=work_p):: sml_ignore_drift_slope1
  real (kind=work_p):: sml_ignore_drift_slope2


  attributes(constant) :: sml_bd_min_r, sml_bd_max_r,sml_bd_min_z,sml_bd_max_z 
  attributes(constant) :: sml_bounce, sml_bounce_zero_weight
  attributes(constant) :: sml_2pi, sml_2pi_wedge_n, ipc_gpu, epc_gpu
  attributes(constant) :: sml_turb_efield, sml_00_efield
  attributes(constant) :: sml_f0_1_Ln, sml_f0_1_Lt, sml_f0_1_Lt_e
  attributes(constant) :: sml_deltaf_f0_mode, sml_dwdt_fix_bg
  attributes(constant) :: sml_dwdt_exb_only, sml_deltaf, sml_electron_on, sml_extra_dwdt, sml_neutral
  attributes(constant) :: sml_inpsi, sml_outpsi
  attributes(constant) :: sml_comm, sml_comm_null 
  attributes(constant) :: sml_totalpe,sml_mype, sml_nthreads, sml_gstep, sml_sheath_mode, sml_ipc
  attributes(constant) :: sml_time, sml_bp_sign, sml_bt_sign
  attributes(constant) :: sml_ignore_drift_near_wall, sml_ignore_drift_r0, sml_ignore_drift_z0, &
                        sml_ignore_drift_slope1, sml_ignore_drift_slope2
  contains

  attributes(host) &
  subroutine update_device_sml()
  use sml_module, only :  &
     sml_bd_min_r_host => sml_bd_min_r, &
     sml_bd_max_r_host => sml_bd_max_r, &
     sml_bd_min_z_host => sml_bd_min_z, &
     sml_bd_max_z_host => sml_bd_max_z, &
     sml_bounce_host => sml_bounce, &
     sml_bounce_zero_weight_host => sml_bounce_zero_weight, &
     sml_2pi_host => sml_2pi, &
     sml_2pi_wedge_n_host => sml_2pi_wedge_n, &
     sml_00_efield_host => sml_00_efield, &
     sml_turb_efield_host => sml_turb_efield, &
     sml_f0_1_Ln_host => sml_f0_1_Ln, &
     sml_f0_1_Lt_host => sml_f0_1_Lt, &
     sml_f0_1_Lt_e_host => sml_f0_1_Lt_e, &
     sml_deltaf_f0_mode_host => sml_deltaf_f0_mode, &
     sml_deltaf_host => sml_deltaf, &
     sml_electron_on_host => sml_electron_on, &
     sml_extra_dwdt_host => sml_extra_dwdt, &
     sml_neutral_host => sml_neutral, &
     sml_dwdt_exb_only_host => sml_dwdt_exb_only, &
     sml_dwdt_fix_bg_host => sml_dwdt_fix_bg, &
     sml_comm_host => sml_comm, &
     sml_comm_null_host => sml_comm_null, &
     sml_totalpe_host => sml_totalpe, &
     sml_mype_host => sml_mype, &
     sml_nthreads_host => sml_nthreads, &
     sml_gstep_host => sml_gstep, &
     sml_ipc_host => sml_ipc, &
     sml_sheath_mode_host => sml_sheath_mode, &
     sml_time_host => sml_time, &
     sml_inpsi_host => sml_inpsi, &
     sml_outpsi_host => sml_outpsi, &
     sml_bp_sign_host => sml_bp_sign, &
     sml_bt_sign_host => sml_bt_sign, &
     sml_ignore_drift_near_wall_host => sml_ignore_drift_near_wall, &
     sml_ignore_drift_r0_host => sml_ignore_drift_r0, &
     sml_ignore_drift_z0_host => sml_ignore_drift_z0, &
     sml_ignore_drift_slope1_host => sml_ignore_drift_slope1, &
     sml_ignore_drift_slope2_host => sml_ignore_drift_slope2

  implicit none 

  sml_bd_min_r = sml_bd_min_r_host
  sml_bd_max_r = sml_bd_max_r_host
  sml_bd_min_z = sml_bd_min_z_host
  sml_bd_max_z = sml_bd_max_z_host

  sml_bounce = sml_bounce_host
  sml_bounce_zero_weight = sml_bounce_zero_weight_host

  sml_2pi = sml_2pi_host
  sml_2pi_wedge_n = sml_2pi_wedge_n_host

  sml_00_efield = sml_00_efield_host
  sml_turb_efield = sml_turb_efield_host

  sml_f0_1_Ln = sml_f0_1_Ln_host
  sml_f0_1_Lt = sml_f0_1_Lt_host
  sml_f0_1_Lt_e =  sml_f0_1_Lt_e_host

  sml_deltaf_f0_mode = sml_deltaf_f0_mode_host
  sml_deltaf = sml_deltaf_host
  sml_electron_on = sml_electron_on_host
  sml_extra_dwdt = sml_extra_dwdt_host
  sml_neutral = sml_neutral_host
  sml_dwdt_exb_only = sml_dwdt_exb_only_host
  sml_dwdt_fix_bg = sml_dwdt_fix_bg_host

  sml_comm = sml_comm_host
  sml_comm_null = sml_comm_null_host
  sml_totalpe = sml_totalpe_host
  sml_mype = sml_mype_host
  sml_nthreads = sml_nthreads_host
  sml_gstep = sml_gstep_host
  sml_ipc = sml_ipc_host
  ipc_gpu=sml_ipc_host
  sml_sheath_mode = sml_sheath_mode_host
  sml_inpsi = sml_inpsi_host
  sml_outpsi = sml_outpsi_host

  sml_time = sml_time_host
  sml_bp_sign = sml_bp_sign_host
  sml_bt_sign = sml_bt_sign_host

  sml_ignore_drift_near_wall = sml_ignore_drift_near_wall_host
  sml_ignore_drift_r0 = sml_ignore_drift_r0_host
  sml_ignore_drift_z0 = sml_ignore_drift_z0_host
  sml_ignore_drift_slope1 = sml_ignore_drift_slope1_host
  sml_ignore_drift_slope2 = sml_ignore_drift_slope2_host
  if(sml_extra_dwdt_host) print *, "Warning: GPU implementation is not correct for this case"
  return
  end subroutine update_device_sml



  attributes(host) &
  subroutine update_host_sml()
  use sml_module, only :  &
     sml_bd_min_r_host => sml_bd_min_r, &
     sml_bd_max_r_host => sml_bd_max_r, &
     sml_bd_min_z_host => sml_bd_min_z, &
     sml_bd_max_z_host => sml_bd_max_z, &
     sml_bounce_host => sml_bounce, &
     sml_bounce_zero_weight_host => sml_bounce_zero_weight, &
     sml_2pi_host => sml_2pi, &
     sml_2pi_wedge_n_host => sml_2pi_wedge_n, &
     sml_00_efield_host => sml_00_efield, &
     sml_turb_efield_host => sml_turb_efield, &
     sml_f0_1_Ln_host => sml_f0_1_Ln, &
     sml_f0_1_Lt_host => sml_f0_1_Lt, &
     sml_f0_1_Lt_e_host => sml_f0_1_Lt_e, &
     sml_deltaf_f0_mode_host => sml_deltaf_f0_mode, &
     sml_deltaf_host => sml_deltaf, &
     sml_electron_on_host => sml_electron_on, &
     sml_dwdt_exb_only_host => sml_dwdt_exb_only, &
     sml_dwdt_fix_bg_host => sml_dwdt_fix_bg, &
     sml_comm_host => sml_comm, &
     sml_comm_null_host => sml_comm_null, &
     sml_totalpe_host => sml_totalpe, &
     sml_mype_host => sml_mype, &
     sml_time_host => sml_time, &
     sml_inpsi_host => sml_inpsi, &
     sml_outpsi_host => sml_outpsi, &
     sml_bp_sign_host => sml_bp_sign, &
     sml_bt_sign_host => sml_bt_sign, &
     sml_ipc_host => sml_ipc,         &
     sml_ignore_drift_near_wall_host => sml_ignore_drift_near_wall, &
     sml_ignore_drift_r0_host => sml_ignore_drift_r0, &
     sml_ignore_drift_z0_host => sml_ignore_drift_z0, &
     sml_ignore_drift_slope1_host => sml_ignore_drift_slope1, &
     sml_ignore_drift_slope2_host => sml_ignore_drift_slope2


  implicit none

  sml_bd_min_r_host = sml_bd_min_r
  sml_bd_max_r_host = sml_bd_max_r
  sml_bd_min_z_host = sml_bd_min_z
  sml_bd_max_z_host = sml_bd_max_z

  sml_bounce_host = sml_bounce
  sml_bounce_zero_weight_host = sml_bounce_zero_weight

  sml_2pi_host = sml_2pi
  sml_2pi_wedge_n_host = sml_2pi_wedge_n

  sml_00_efield_host = sml_00_efield
  sml_turb_efield_host = sml_turb_efield

  sml_f0_1_Ln_host = sml_f0_1_Ln
  sml_f0_1_Lt_host = sml_f0_1_Lt
  sml_f0_1_Lt_e_host =  sml_f0_1_Lt_e

  sml_deltaf_f0_mode_host = sml_deltaf_f0_mode
  sml_deltaf_host = sml_deltaf
  sml_electron_on_host = sml_electron_on
  sml_dwdt_exb_only_host = sml_dwdt_exb_only
  sml_dwdt_fix_bg_host = sml_dwdt_fix_bg

  sml_comm_host = sml_comm
  sml_comm_null_host = sml_comm_null
  sml_totalpe_host = sml_totalpe
  sml_mype_host = sml_mype

  sml_inpsi_host = sml_inpsi
  sml_outpsi_host = sml_outpsi

  sml_time_host = sml_time
  sml_bp_sign_host = sml_bp_sign
  sml_bt_sign_host = sml_bt_sign

  sml_ignore_drift_near_wall_host = sml_ignore_drift_near_wall
  sml_ignore_drift_r0_host = sml_ignore_drift_r0
  sml_ignore_drift_z0_host = sml_ignore_drift_z0
  sml_ignore_drift_slope1_host = sml_ignore_drift_slope1
  sml_ignore_drift_slope2_host = sml_ignore_drift_slope2

  return
  end subroutine update_host_sml

end module sml_module_gpu
