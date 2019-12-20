attributes(global) &
subroutine push_kernel_gpu(istep,ipc,phase0, &
                   diag_on,dt_now,gpu_ibegin,gpu_iend)
use cudafor
use dimensions_mod_gpu, only : nthreads_dim
use sml_module_gpu, only : sml_deltaf, sml_nrk, sml_bounce
use ptl_module_gpu, only : ptl_nphase, ptl_gid_gpu, ptl_ph_gpu
use eq_module_gpu, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z
use precision_mod_gpu

implicit none

real(kind=work_p) :: phase0(ptl_nphase,*)
real(kind=work_p) :: phase0_i(ptl_nphase)
integer, value :: istep, ipc !! RK4 index
logical, value :: diag_on
real(kind=work_p), value :: dt_now
integer, value :: gpu_ibegin,gpu_iend
real (kind=work_p) :: phi, p(3)
real(kind=work_p), dimension(ptl_nphase) :: old_phase, new_phase
integer :: i, ith, rtn, itr 
integer, parameter :: idebug = 0

  ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) +  &
      ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
      (blockDim%x * blockDim%y)


!  if (ith.eq.1)  write(*, *) 'push_i threadIdx%x=',threadIdx%x,dt_now


  if ((ith < 1).or.(ith > nthreads_dim)) then
      return
  endif

!  if (diag_on) then
!    diag_1d_f_pv1(:,:,:,ith) = 0
!    if (sml_deltaf) then
!       diag_1d_df_pv1(:,:,:,ith) = 0
!    endif
!  endif




   do i=gpu_ibegin+(ith-1),gpu_iend,nthreads_dim
  
        ! for alive particles only
        if(ptl_gid_gpu(i)>0) then
           
           !******************************************************
           ! actual particle push           
           !******************************************************
           phase0_i = phase0(:,i)
           call push_single_gpu(i,phase0_i, &
                                new_phase,dt_now,ith,diag_on,itr,p)
           phase0(:,i) = phase0_i 
           
           ! check r-z boundary validity and update psi variables
           if(new_phase(1)<eq_min_r .or. new_phase(1)>eq_max_r .or.  &
              new_phase(2)<eq_min_z .or. new_phase(2)>eq_max_z)then
!              call remove_particle_gpu(sp%ptl(i)%gid, &
!                    sp%ptl(i)%ph,sp%ptl(i)%ct,sp%type,i,-1)
               call remove_particle_gpu(i,-1,tb)
!                if (idebug >= 1) then
!                 print *, 'particle eliminated due to rz_outside :',  &
!                 i, sml_mype, sp%type, sp%ptl(i)%gid, &
!                 new_phase(1),new_phase(2)
!                endif
           else                            
              ! bounce 
              if(ipc==sml_nrk .and. sml_bounce/=0) then
                 old_phase(:)=phase0(:,i)
                 call bounce_gpu(new_phase,old_phase,rtn)
                 if(rtn<0)  then
                     call remove_particle_gpu( i,-2,tb)
                 endif
              endif
              
              !******************************************************
              ! time advance one step
              !******************************************************
              ptl_ph_gpu(:,i)= new_phase(:)
              
           endif
           phi=ptl_ph_gpu(3,i)
           if( phi >= sml_2pi_wedge_n .or. phi < 0D0 ) then
           ptl_ph_gpu(3,i)=modulo(phi,sml_2pi_wedge_n)
           endif

        endif
    enddo
  
end subroutine push_kernel_gpu
