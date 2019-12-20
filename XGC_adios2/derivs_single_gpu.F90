!obtain derivatives of phase variable : actual calculation is done in derivs_sp. 
! prepare E-field and B-field
attributes(device) &
subroutine derivs_single_gpu(ptli,dy,i,time,fld,ith,diag_on,itr,p,tb)

  use sml_module_gpu
  use fld_module, only : fld_type
  use ptl_module_gpu, only : type_gpu
  use precision_mod_gpu

  implicit none
  integer,intent(in) :: ptli, itr
  real (kind=work_p), intent(inout) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (kind=work_p), intent(in) :: time, p(3)
  type(fld_type), intent(inout) :: fld
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  type(tbuf), intent(in) :: tb
  !
  logical rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)   ! variables for diagnosis 

  integer, parameter :: idebug = 0

  ! Save space information
  fld%r=tb%ph(1)
  fld%z=tb%ph(2)
  fld%phi=tb%ph(3)

  ! obtain B-field information 
  call field_gpu(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field for each particle : use position information from 'charge'
     call efield_gpu(i,fld,time,itr,p)
!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp_elec_gpu(fld,time,ptli,dy,diag_on,vf_diag,tb)
     
  else
     call remove_particle_gpu(i,-1,tb)
!     call remove_particle_gpu(sp%ptl(i)%gid, &
!                              sp%ptl(i)%ph, sp%ptl(i)%ct, &
!                              sp%type,  i, -1)

!     if (idebug >= 1) then
!       print *, 'particle eliminated due to rz_outside',  &
!                i, sml_mype, sp%type, sp%ptl(i)%gid
!     endif
  endif
     

  if(diag_on) call diag_1d_port1_gpu(ptli,dy,type_gpu,vf_diag,ith,tb)

end subroutine derivs_single_gpu
