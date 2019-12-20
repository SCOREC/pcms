attributes(global) &
subroutine pushe_kernel_gpu(istep,ipc,phase0, &
                   diag_on,dt_now,gpu_ibegin,gpu_iend,use_perm,iperm)
use cudafor
!use dimensions_mod_gpu !, only : nthreads_dim
use sml_module_gpu, only : sml_deltaf, sml_nrk, sml_bounce, sml_2pi_wedge_n
use ptl_module_gpu, only : ptl_nphase, ptl_gid_gpu, ptl_ph_gpu, tr_save_gpu
use grid_class_gpu, only : grid_delta_phi, search_tr2_gpu, search_tr_check_guess_gpu
use eq_module_gpu, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z
use precision_mod_gpu

implicit none

real(kind=work_p), intent(inout) :: phase0(gpu_ibegin:gpu_iend,ptl_nphase)
integer, intent(in) :: iperm(gpu_ibegin:gpu_iend)
real(kind=work_p) :: phase0_i(ptl_nphase)
integer, value :: istep, ipc !! RK4 index
logical, value, intent(in) :: diag_on
logical, value, intent(in) :: use_perm
real(kind=work_p), value, intent(in) :: dt_now
integer, value, intent(in) :: gpu_ibegin,gpu_iend
real(kind=work_p), dimension(ptl_nphase) :: old_phase, new_phase
integer :: i, k, iv, ith, rtn, j 
integer, parameter :: idebug = 0 

real (kind=work_p) :: phi_tmp
real (kind=work_p) :: phi_mid, x(2), phi, mu, rho, xff(2)
real (kind=work_p) :: p(3)
integer :: itr, ip, nthreads_dim
type(tbuf) :: tb

  ith = 1+ ((threadIdx%x-1) + (threadIdx%y-1)*blockDim%x) +  &
      ((blockIdx%x-1) + (blockIdx%y-1)*gridDim%x )* &
      (blockDim%x * blockDim%y)
  nthreads_dim = blockDim%x*gridDim%x

!  if (ith.eq.1) print*, 'threadIdx%x=', ptl_ph_gpu(10,2)


  if ((ith < 1).or.(ith > nthreads_dim)) then
      return
  endif

   do j=gpu_ibegin+(ith-1),gpu_iend,nthreads_dim
        ! optionally use an on-the-fly gather from iperm
        if (use_perm) then
          i = iperm(j)
        else
          i = j
        endif
        tb%ph=ptl_ph_gpu(i,:)
        tb%ct=ptl_ct_gpu(i,:)
        ! for alive particles only
        if(ptl_gid_gpu(i)>0) then

           ! get proper toroidal angle index and weight
           x=tb%ph(1:2)
           phi=tb%ph(3)
           phi_mid=(floor(phi/grid_delta_phi) + 0.5_work_p) * grid_delta_phi
!           mu=sp%ptl(i)%ct(pim)
!           rho=gyro_radius(x,mu)  !gyro radius
!           sp%rhoi(i)=rho

             ! get field following posision at 1/2 angle
 
!   if(idebug.eq.1)print *,'before field_following_pos2_gpu'

           call field_following_pos2_gpu(x,phi,phi_mid,xff)

#ifdef USE_TR_CHECK
           call search_tr_check_guess_gpu(xff,tr_save_gpu(i),itr,p)
           if (itr < 0) then
             call search_tr2_gpu(xff,itr,p)
           endif
#else
           call search_tr2_gpu(xff,itr,p)
#endif


           !remove particle or sheath calculation
           if(itr<0) then
              if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                 call remove_particle_gpu(i,-1,tb)
              else
                 call sheath_calculation_gpu(i,ipc,0,itr,p,ith,tb)
              endif
           endif

#ifdef USE_TR_CHECK
            tb%tr = itr
#endif
!           tr_save_gpu(i)=itr
!           p_save_gpu(:,i)=p
 
          
           !******************************************************
           ! actual particle push           
           !******************************************************

       select case(ipc)
       case(1)
             phase0(i,:)=tb%ph(:)
       end select


           phase0_i = phase0(i,:)
           call push_single_gpu(i,phase0_i, &
                                new_phase,dt_now,ith,diag_on,itr,p,tb)
!           phase0(:,i) = phase0_i 
           
           ! check r-z boundary validity and update psi variables
           if(new_phase(1)<eq_min_r .or. new_phase(1)>eq_max_r .or.  &
              new_phase(2)<eq_min_z .or. new_phase(2)>eq_max_z)then

              call remove_particle_gpu(i,-1,tb)
                if (idebug >= 1) then
!                 print *, 'particle eliminated due to rz_outside :' ,  &
!                 i, sml_mype, sp%type, sp%ptl(i)%gid, &
!                 new_phase(1),new_phase(2)
                endif
           else                            
              ! bounce 
              if(ipc==sml_nrk .and. sml_bounce/=0) then
                 old_phase(:)=phase0(i,:)
                 call bounce_gpu(new_phase,old_phase,rtn)
                 if(rtn<0)  then
                     call remove_particle_gpu( i,-2,tb)
                 endif
              endif
              
              !******************************************************
              ! time advance one step
              !******************************************************
              tb%ph(:) = new_phase(:)
              
           endif

           if(tb%ph(3)>= sml_2pi_wedge_n .or. tb%ph(3)< 0D0 ) then
              tb%ph(3)=modulo(tb%ph(3),sml_2pi_wedge_n)
           endif
           ptl_ph_gpu(i,:) = tb%ph(:)
#ifdef USE_TR_CHECK
           tr_save_gpu(i) = tb%tr
#endif

        endif
    enddo
  
end subroutine pushe_kernel_gpu
