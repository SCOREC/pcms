! single particle push -- get new_phase using rk4 with initial E-field
attributes(device) &
subroutine push_single_gpu(i,y,new_phase,dt,ith,diag_on,itr,p,tb)
  use sml_module_gpu
  use ptl_module_gpu, only : ptl_gid_gpu, piw1, piw2
  use fld_module
  use precision_mod_gpu
  !gpu use perf_monitor
  implicit none

  integer, intent(in) :: i, itr
  real (kind=work_p), intent(in) :: y(ptl_nphase), p(3)
  real (kind=work_p), intent(inout) :: new_phase(ptl_nphase)
  real (kind=work_p), intent(in) :: dt
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  type(tbuf) :: tb
  !
  integer :: ptli
  type(fld_type) :: fld
  real (kind=work_p) :: dy(ptl_nphase),dyt(ptl_nphase),dym(ptl_nphase)
  integer :: rtn,j,i1
  logical :: rz_outside
  real (kind=work_p) :: vf_diag(sml_n_vf_diag)
  real (kind=work_p) , external :: psi_interpol
  !gpu character (len=5) :: err_str(0:1)
  real (kind=work_p) :: hdt, time, th  ! -- time should be input 
  real (kind=work_p) :: bp,B, E_mag(3)

  !gpu err_str(1)='ion'
  !gpu err_str(0)='elec'


  !pushing particle "i"


  time=sml_time ! need to update -- rk4 time ###

  hdt=dt*0.5D0
  th=time + hdt

  if(ptl_gid_gpu(i)>0) then
     
     !set ptli -- ptli%ct does not change
!     ptli%ct=sp%ptl(i)%ct
!     ptli%ph=sp%ptl(i)%ph
     ptli = i    
     !get derivs with updating E-field information - assigned on fld
     !diag-ports are called when diag_on is .true.
    call derivs_single_gpu(ptli,dy,i,time,fld,ith,diag_on,itr,p,tb)


#ifdef PURE_RK2
     ! Only simple RK2
     
     new_phase = y + dt * dy
     call restrict_weight_gpu(new_phase(piw1:piw2))
#else

     ! RK2 - RK4 hybrid -- time advance with RK4 with time-constant E-field
     
     ! get E-field in magnetic field
     bp=sqrt(fld%br**2 + fld%bz**2)     
     B=sqrt(bp**2 + fld%Bphi**2 )
     E_mag(2)=(fld%Er*fld%Br + fld%Ez*fld%Bz)/bp
     E_mag(3)=(E_mag(2)*bp   + fld%Ephi*fld%Bphi)/B    ! parallel field
     E_mag(1)=(fld%Er*fld%dpsidr + fld%Ez*fld%dpsidz )/(fld%dpsidr**2 + fld%dpsidz**2)
     
     ! get derivs with existing E-field
     tb%ph(:)=y
     call derivs_single_with_e_gpu(ptli,dy ,i,time,E_mag,tb)
     
     tb%ph(:) = y + hdt * dy     
     call derivs_single_with_e_gpu(ptli,dyt,i,th,E_mag,tb)
     
     tb%ph(:) = y + hdt * dyt
     call derivs_single_with_e_gpu(ptli,dym,i,th,E_mag,tb)
     
     tb%ph(:) = y + dt * dym
     dym = dyt + dym
     call derivs_single_with_e_gpu(ptli,dyt,i,time+dt,E_mag,tb)
     
     ! Obtain new_phase
     new_phase = y + dt/6D0 * ( dy + dyt + 2D0*dym )
!     call restrict_weight_gpu(new_phase(piw1:piw2))
#endif     
  endif
    
end subroutine push_single_gpu
