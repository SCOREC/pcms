attributes(device) &
subroutine efield_gpu(i,fld,time,itr,p)
  use sml_module_gpu,  only :  &
     sml_00_efield, sml_turb_efield
  use fld_module
  !use eq_module
  use precision_mod_gpu
  use psn_class_gpu, only : efield_gk_elec_gpu
  implicit none
  type(fld_type), intent(inout) :: fld
  integer, intent(in) :: i,itr
  real (kind=work_p), intent(in) :: time, p(3)

  
  if(sml_00_efield .or. sml_turb_efield) then
     !Gyrokinetic E
     call efield_gk_elec_gpu(i,fld,itr,p)
  else
       fld%Er=0D0
       fld%Ez=0D0
       fld%Ephi=0D0
       fld%Epot=0D0
       fld%Er00=0D0
       fld%Ez00=0D0
       fld%ddpotdt=0D0
  endif

end subroutine efield_gpu
