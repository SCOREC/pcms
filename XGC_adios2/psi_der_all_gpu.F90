attributes(device) &
subroutine psi_der_all_gpu(r,z,ret)
use bicub_mod_gpu, only : bicub_interpol2
use precision_mod_gpu
implicit none
  real (kind=work_p) :: r,z,ret(6)
integer, parameter :: i00 = 1
integer, parameter :: i10 = 2
integer, parameter :: i01 = 3
integer, parameter :: i20 = 4
integer, parameter :: i11 = 5
integer, parameter :: i02 = 6

call bicub_interpol2(r,z,                                           &
  ret(i00),ret(i10),ret(i01),                                          &
  ret(i11),ret(i20),ret(i02))

end subroutine psi_der_all_gpu
