attributes(device) &
subroutine remove_particle_gpu( i,flag, tb)
  use ptl_module_gpu, only : ptl_gid_gpu,  &
                             ptl_nphase,type_gpu,ptl_nconst,pim,piw0
  use sml_module_gpu, only : sml_mype
  use neu_module_gpu, only : neu_weight_sum_lost_gpu
  use precision_mod_gpu
  implicit none

  integer, intent(in) :: i,flag
  type (tbuf), intent(in) :: tb

!  integer :: sp%type
!  integer(8) :: sp%ptl(i)%gid
!  real(8), dimension(ptl_nphase) :: sp%ptl(i)%ph
!  real(8), dimension(ptl_nconst) :: sp%ptl(i)%ct
  !error diagnosis
  real (kind=work_p) :: b_interpol_gpu, b,ekin,pitch
  integer :: ip
  real (kind=work_p) :: dummy
  integer, parameter :: idebug = 0


!  sp%ptl(i)%gid = sp%ptl(i)%gid
!  sp%ptl(i)%ph = sp%ptl(i)%ph
!  sp%ptl(i)%ct = sp%ptl(i)%ct
!  sp%type = sp%type
  
  if(ptl_gid_gpu(i) <= 0 ) then
!    if (idebug >= 1) print *, 'minus gid particle in remove_particle_gpu'
    return
  endif
  ptl_gid_gpu(i)=-ptl_gid_gpu(i)

  !### debug
   if(ptl_gid_gpu(i) > 0 ) then
     if (idebug >= 1) print *, 'something wrong in remove_particle_gpu'
  endif

!if (idebug >= 1) print *, 'sml_neutral', sml_neutral

  if(flag==-2 .and. sml_mype==0) then
     b=b_interpol_gpu(tb%ph(1),tb%ph(2),0.0_work_p)
     call rho_mu_to_ev_pitch2_gpu(tb%ph(4), &
                 tb%ct(pim), b,ekin,pitch,type_gpu)
!     write(400+sp%type,*) sp%ptl(i)%ph(1:2),&
!          ekin, pitch,&          
!          sp%phase0(1,i),sp%phase0(2,i) 
  endif

  if(sml_neutral .and. flag==-1.and.(ptl_gid_gpu(i)<0).and.(type_gpu==1))then
!     neu_weight_sum_lost_gpu(ith)=neu_weight_sum_lost_gpu(ith) + ptl_ct_gpu(i,piw0)
     ip = mod( i, size(neu_weight_sum_lost_gpu)) + 1
     dummy = atomicAdd( neu_weight_sum_lost_gpu(ip), tb%ct(piw0))
  endif

  if(flag==-1 .and. sml_mype==0) then
!      write(450+sp%type,*) sp%phase(1,i), sp%phase(2,i),&
!           ekin, pitch,&          
!           sp%phase0(1,i),sp%phase0(2,i) 
  endif

!   if(flag==-3) then
!       print *,  'unexpected search fail in search_pid', &
!                 sp%phase(1,i), sp%phase(2,i)
!   endif
end subroutine remove_particle_gpu
