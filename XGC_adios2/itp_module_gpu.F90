module itp_module_gpu
use precision_mod_gpu
implicit none

real(kind=work_p) :: itp_min_psi,itp_max_psi
attributes(constant) :: itp_min_psi, itp_max_psi

contains
  attributes(host) &
  subroutine  update_device_itp()
  use itp_module, only : &
    itp_min_psi_host => itp_min_psi, &
    itp_max_psi_host => itp_max_psi

  itp_min_psi = itp_min_psi_host
  itp_max_psi = itp_max_psi_host

  return
  end subroutine update_device_itp

end module itp_module_gpu
