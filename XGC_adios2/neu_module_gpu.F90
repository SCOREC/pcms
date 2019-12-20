module neu_module_gpu
use neu_module, only:  neu_weight_sum_lost
use precision_mod_gpu
use cudafor
implicit none

  integer, parameter :: neu_weight_sum_lost_gpu_dim = 256
  real (kind=work_p) :: neu_weight_sum_lost_gpu(neu_weight_sum_lost_gpu_dim)
  attributes(device) :: neu_weight_sum_lost_gpu

contains

  attributes(host) &
  subroutine update_device_neu()
  implicit none

  integer :: ierr,icount


  neu_weight_sum_lost_gpu = 0

  return
  end subroutine update_device_neu  

  attributes(host) &
  subroutine update_host_neu()
  implicit none

  real(kind=work_p) :: tmp_neu_weight_sum_lost(neu_weight_sum_lost_gpu_dim)
  integer :: lb


  if (allocated(neu_weight_sum_lost)) then
   if (size(neu_weight_sum_lost) >= 1) then
    tmp_neu_weight_sum_lost = neu_weight_sum_lost_gpu
    lb = lbound(neu_weight_sum_lost,1)
    neu_weight_sum_lost(lb) = neu_weight_sum_lost(lb) +  &
                            sum(tmp_neu_weight_sum_lost)
   endif
  endif

  return
  end subroutine update_host_neu
end module neu_module_gpu
