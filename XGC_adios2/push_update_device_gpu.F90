attributes(host) &
subroutine push_update_device_gpu( )
 use sml_module, only : sml_mype
! use bicub_mod, only : psi_bicub

 use sml_module_gpu, only : update_device_sml 
 use eq_module_gpu, only : update_device_eq 
 use itp_module_gpu, only : update_device_itp 
 use one_d_cub_mod_gpu, only : update_device_one_d 
 use bicub_mod_gpu, only : update_device_bicub 
 use ptl_module_gpu, only : update_device_ptl 
 use diag_module_gpu, only : update_device_diag 
 use bnc_module_gpu, only : update_device_bnc 

 integer, parameter :: idebug = 0

 if (idebug >= 1) then
!   if (sml_mype == 0) print*,'before update_device_sml()'
 endif
! call update_device_sml()

 if (idebug >= 1)  then
!   if (sml_mype == 0) print*,'before update_device_eq()'
 endif
 call update_device_eq()

 if (idebug >= 1) then
!   if (sml_mype == 0) print*,'before update_device_itp()'
 endif
 call update_device_itp()

 if (idebug >= 1) then
!   if (sml_mype == 0) print*,'before update_device_one_d()'
 endif
 call update_device_one_d()

 if (idebug >= 1) then
!   if (sml_mype == 0) print*,'before update_device_bicub()'
 endif
 call update_device_bicub()

 if (idebug >= 1) then
!   if (sml_mype == 0) print*,'before update_device_ptl()'
 endif
! call update_device_ptl( )

 if (idebug >= 1) then
!   if (sml_mype) print*,'before update_device_diag()'
 endif
 call update_device_diag()

 if (idebug >= 1) then
!   if (sml_mype) print*,'before update_device_bnc()'
 endif
 call update_device_bnc()

 return
end subroutine push_update_device_gpu
