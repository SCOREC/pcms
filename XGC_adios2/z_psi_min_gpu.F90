attributes(device) &
real (kind=8) function z_psi_min_gpu(r)
  use sml_module_gpu, only : sml_mype
  use bnc_module_gpu, only : bnc_min_r, bnc_nr, bnc_dr, bnc_z_psi_min
  use precision_mod_gpu
  implicit none
  real (kind=work_p),intent(in) :: r
  integer :: i
  real (kind=work_p) :: aa,bb
  
  i=int((r-bnc_min_r)/bnc_dr) +1
  i=min(bnc_nr-1,max(1,i))
  bb=(r-bnc_min_r)/bnc_dr + 1D0 -i
  aa=1D0-bb
!  print *, 'zaa', r, i, aa,bb,bnc_z_psi_min(i)
  z_psi_min_gpu=bnc_z_psi_min(i)*aa + bnc_z_psi_min(i+1)*bb

end function z_psi_min_gpu
