module bnc_module_gpu
  use bnc_module, only : bnc_nr_host => bnc_nr
  use precision_mod_gpu
  integer, parameter :: bnc_nr=bnc_nr_host
  real (kind=work_p) :: bnc_min_r, bnc_max_r 
  real (kind=work_p) :: bnc_dr
  real (kind=work_p) :: bnc_z_psi_min(bnc_nr)
  ! real (8) :: bnc_arg_pass_r

  attributes(constant) :: bnc_min_r, bnc_max_r
  attributes(constant) :: bnc_dr, bnc_z_psi_min

  contains

  attributes(host) &
  subroutine update_device_bnc()
  use bnc_module, only : &
    bnc_z_psi_min_host => bnc_z_psi_min, &
    bnc_dr_host => bnc_dr, &
    bnc_min_r_host => bnc_min_r, &
    bnc_max_r_host => bnc_max_r

    bnc_z_psi_min = bnc_z_psi_min_host

    bnc_dr = bnc_dr_host
    bnc_min_r = bnc_min_r_host
    bnc_max_r = bnc_max_r_host

  end subroutine update_device_bnc

end module bnc_module_gpu
