attributes(host) &
subroutine pushe_gpu(istep,ihybrid,ncycle,grid,psn,sp,diag_on_input)
  use grid_class
  use psn_class
  use ptl_module
  use perf_monitor
  use sml_module
  use omp_module , only : split_indices
  use eq_module
  use cudafor
  use precision_mod_gpu

  use ptl_module_gpu, only :  &
     phase0_gpu, ptl_ph_gpu,        &
#ifdef USE_TR_CHECK
     tr_save_gpu, &
#endif
     ptl_ct_gpu,        &
!     iperm_gpu, xstart_gpu,   &
     update_device_species_type, &
     update_host_species_type

  use psn_class_gpu, only : &
    update_device_psn_type

!  use bicub_mod_gpu, only : &
!    update_device_bicub

  use gen_perm_gpu_mod, only : &
    gen_perm_gpu

  use reorder_gpu_mod, only : &
    reorder1d_gpu,       &
    reorder2d_gpu

  use util_mod_gpu, only : get_gpu_streamid

  implicit none
  integer, intent(in) :: istep, ihybrid
  integer, intent(in) :: ncycle
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type),target :: sp
  logical :: diag_on,need_copyin,need_copyout,diag_on_input
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  integer :: istatus
  real(kind=work_p) :: gpu_ratio
  integer :: perm_gpu_freq, sort_gpu_freq
  logical :: perm_gpu, sort_gpu, access_thru_perm
  character(len=255) :: sml_gpu_ratio_str
  integer :: gpu_ibegin, gpu_iend
  integer :: cpu_ibegin, cpu_iend

  integer :: ierr,stream1
  integer :: nblocks
  type(dim3) :: tgrid, tblock
  integer :: lb1,ub1,lb2,ub2
    logical  :: rz_outside
!  type(ptl_type), dimension(:), pointer ::  ptl
 
  real (kind=work_p) :: dt_now,dt,time_now,new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer, parameter :: idebug =0 
  real(kind=work_p), parameter :: dzero = 0.0_work_p

  integer :: icycle
  integer :: epc, i, rtn

  real (kind=work_p) :: phi_mid, x(2), phi, mu, rho, xff(2)
  real (kind=work_p) :: p(3)
  integer :: itr, ip
!  real(kind=work_p) :: phase_tmp(ptl_nphase)

  attributes(device) :: d_ibegin,d_iend !, phase_tmp


! -------------------------------------------
! local variables related to particle sorting
! -------------------------------------------
  logical, parameter :: use_sort_particles = .true.
  logical, parameter :: use_reorder_array = .false.
  ! Sorting by triangle seems to be somewhat faster
  logical, parameter :: use_sort_by_triangle = .true.
  integer :: ilo,ihi,jlo,jhi
  integer :: nx, ny, n, xydim
  real(kind=work_p) :: xmin,ymin, inv_dx,inv_dy
  integer*8, allocatable, dimension(:) :: gid
  real(kind=work_p), allocatable, dimension(:) :: xcoord,ycoord 
  real(8), allocatable, dimension(:,:) :: tphase0
  type(ptl_type), allocatable, dimension(:) :: tptl
  integer, allocatable, dimension(:) :: iperm
  integer, allocatable, dimension(:) :: iperm_gpu, xstart_gpu
  integer :: itotal
  integer :: mm, nn, lld
  integer :: streamid
  attributes(device) :: d_ibegin,d_iend, iperm_gpu, xstart_gpu

  if(sp%num==0) return  ! nothing to push

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

  if (use_sort_particles) then
    call t_startf("pushe_sort_particles")
! --------------
! sort particles
! --------------
    ilo = lbound( grid%guess_table,1)
    ihi = ubound( grid%guess_table,1)
    jlo = lbound( grid%guess_table,2)
    jhi = ubound( grid%guess_table,2)

    allocate( iperm(sp%num))

    if (use_sort_by_triangle) then
      call gen_perm_tri(sp%ptl(:), grid, iperm, sp%num)
    else 
      allocate( gid(sp%num), xcoord(sp%num), ycoord(sp%num) )

!$omp parallel do  private(ith,i)
      do ith=1,sml_nthreads
        do i=i_beg(ith),i_end(ith)
           gid(i) = sp%ptl(i)%gid
           xcoord(i) = sp%ptl(i)%ph(1)
           ycoord(i) = sp%ptl(i)%ph(2)

           iperm(i) = i
        enddo
      enddo

      call gen_perm( ilo,ihi,jlo,jhi,   &
               grid%guess_min, grid%inv_guess_d, &
               sp%num, gid, xcoord, ycoord, iperm )

!    write(20+sml_mype,*) 'iperm(1:5) ', iperm(1:5)
!    write(20+sml_mype,*) 'iperm(n-5:n) ', iperm(sp%num-1:sp%num)


      deallocate( gid, xcoord, ycoord )
    endif
!   ------------------------
!   rearrange particle data in sp
!    
!   which arrays should be permutated
!   ptl
!
!   should other arrays such as phase0(:) or rhoi(:) 
!   be rearranged????
!   ------------------------


!   --------------
!   permute sp%ptl
!   --------------
    allocate(tptl(sp%num))
!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
        tptl(i) = sp%ptl( iperm(i) )
      enddo
    enddo

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
      sp%ptl(i) = tptl(i)
     enddo
    enddo
    deallocate( tptl )



!  -----------------
!  permute sp%phase0
!  -----------------
    allocate(tphase0(lbound(sp%phase0,1):ubound(sp%phase0,1),sp%num))

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
      do i=i_beg(ith),i_end(ith)
        tphase0(:,i) = sp%phase0(:, iperm(i) )
      enddo
    enddo

!$omp parallel do private(ith,i)
    do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
      sp%phase0(:,i) = tphase0(:,i)
     enddo
    enddo
    deallocate( tphase0 )



    deallocate( iperm )
    call t_stopf("pushe_sort_particles")
   endif

!  ptl => sp%ptl
  ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)

  diag_on=diag_on_input

  if ((sml_gpu_ratio < 0.0) .or. (sml_gpu_ratio > 1.0)) then
     sml_gpu_ratio = 0.7
     call getenv("SML_GPU_RATIO",sml_gpu_ratio_str)
     if (len(trim(sml_gpu_ratio_str)) >= 1) then
        read(sml_gpu_ratio_str,*,iostat=istatus) gpu_ratio
        if (istatus.eq.0) then
           gpu_ratio = max( 0.0_work_p, min(1.0_work_p, dble(gpu_ratio)))
           sml_gpu_ratio =  gpu_ratio
           if (idebug >= 1) then
              if (sml_mype.eq.0) write(*,*) 'sml_gpu_ratio ', sml_gpu_ratio
           endif
        endif
     endif
  endif

! -----------------------------------------------
! suggest always choose gpu_ibegin = 1
! so that transfer to GPU is in contiguous block
! -----------------------------------------------
  gpu_ibegin = 1
  gpu_iend = int( sml_gpu_ratio*sp%num )
!  gpu_iend = min( gpu_iend, ubound(ptl_ph_gpu,1))
  cpu_ibegin = gpu_iend + 1
  cpu_iend = sp%num

!   print *, 'gpu_ibegin', cpu_ibegin,cpu_iend
  dt=sml_dt*0.5_work_p/real(sml_ncycle_half)

if (gpu_ibegin <= gpu_iend .or. cpu_ibegin <= cpu_iend) then
 if(gpu_ibegin <= gpu_iend) then

  call t_startf("push_upd_dev_gpu")
!  call push_update_device_gpu()
!if(sml_mype==0) print *, 'before update_device_ptl'
  call update_device_ptl( )
!if(sml_mype==0) print *, 'before update_device_diag'
  call update_device_diag()
!if(sml_mype==0) print *, 'after update_device_diag'
  call t_stopf("push_upd_dev_gpu")

  call t_startf("upd_dev_psn")
!print *, 'before update_device_psn_type'
  call update_device_psn_type(psn )
!print *, 'after update_device_psn_type'
  call t_stopf("upd_dev_psn")

  call t_startf("upd_dev_species_type")
! print *, 'before update_device_species_type'
  call update_device_species_type(sp, gpu_ibegin,gpu_iend)
  call update_device_sml()
  call update_device_neu()
  call t_stopf("upd_dev_species_type")
!   print *, 'gpu_ibegin', cpu_ibegin,cpu_iend
  xmin = grid_guess_min(1)
  ymin = grid_guess_min(2)
  
  inv_dx = grid_inv_guess_d(1)
  inv_dy = grid_inv_guess_d(2)

  nx = size(grid_guess_table,1)
  ny = size(grid_guess_table,2)

  xydim = size(ptl_ph_gpu,1)
  mm = gpu_iend - gpu_ibegin + 1

  ! -------------------------
  ! launch kernel computation
  ! -------------------------
  tblock%x = block_dim
  tblock%y = 1
  tblock%z = 1
  tgrid%x = (mm + block_dim - 1)/block_dim
  tgrid%y = 1
  tgrid%z = 1

  if (diag_on) then
    call setval_gpu( size(diag_1d_f_pv1),diag_1d_f_pv1,dzero)
    if(diag_eflux_on) call setval_gpu( size(diag_1d_eflux_pv),diag_1d_eflux_pv,dzero)

    if  (sml_deltaf) then
      call setval_gpu( size(diag_1d_df_pv1),diag_1d_df_pv1,dzero)
      if(diag_eflux_on) call setval_gpu( size(diag_2d_dflux_pv),diag_2d_dflux_pv,dzero)
    endif
  endif

  if(diag_heat_on) call setval_gpu(size(diag_heat_pv),diag_heat_pv,dzero)
  if(diag_heat_on) call setval_gpu(size(diag_heat_pv_psi),diag_heat_pv_psi,dzero)

  if (idebug >= 1) then
     if (sml_mype.eq.0) write(*,*) 'size diag_1d_f_pv1', diag_1d_npv1, nx, ny, inv_dx, inv_dy
  endif
 endif

  if (sml_perm_gpu_freq < 1) then
     perm_gpu = .false.
     perm_gpu_freq = ncycle
  else
     perm_gpu = .true.
     perm_gpu_freq = sml_perm_gpu_freq
     allocate(iperm_gpu(gpu_ibegin:gpu_iend), stat=ierr)
     call assert(ierr.eq.0,'alloc(iperm_gpu) ',ierr)
  endif

  if (sml_sort_gpu_freq < 1) then
     sort_gpu = .false.
     sort_gpu_freq = ncycle
  else
     sort_gpu = .true.
     sort_gpu_freq = sml_sort_gpu_freq
  endif

  access_thru_perm = .false.
 do icycle=1, ncycle
!
     if(gpu_ibegin <= gpu_iend) then
      !! Run the iperm generation and/or particle sorting
      if (((perm_gpu) .and. (mod(icycle,perm_gpu_freq).eq.0)) &
        & .or. ((sort_gpu) .and. (mod(icycle,sort_gpu_freq).eq.0))) then

        call t_startf("gen_perm_gpu")

        ! If we're accessing through iperm, then it hang around for the whole run
        if (.not. perm_gpu) then
          allocate(iperm_gpu(gpu_ibegin:gpu_iend), stat=ierr)
          call assert(ierr.eq.0,'alloc(iperm_gpu) ',ierr)
        endif

        if (use_sort_by_triangle) then
          allocate(xstart_gpu(grid%ntriangle+1), stat=ierr)
          call assert(ierr.eq.0,'alloc(xstart_gpu) ',ierr)
          call gen_perm_tri_gpu( gpu_iend-gpu_ibegin+1, grid%ntriangle, &
                              gpu_ibegin, gpu_iend,iperm_gpu, xstart_gpu)
        else
          allocate(xstart_gpu(nx*ny+1), stat=ierr)
          call assert(ierr.eq.0,'alloc(xstart_gpu) ',ierr)

          call gen_perm_gpu( nx, ny, xmin, ymin,   &
                        inv_dx, inv_dy, &
                        gpu_iend-gpu_ibegin+1, ptl_gid_gpu, ptl_ph_gpu(:,1:2), xydim, iperm_gpu, xstart_gpu )
        endif

        access_thru_perm = .true.
        call t_stopf("gen_perm_gpu")

        if ((sort_gpu) .and. (mod(icycle,sort_gpu_freq).eq.0)) then
          call t_startf("reordering_gpu")

          nn  = size( ptl_ph_gpu,2)
          lld = size( ptl_ph_gpu,1)

          if (idebug >= 1) then
           if (sml_mype.eq.0) write(*,*) 'size ptl_ph_gpu',size( ptl_ph_gpu,2), size( ptl_ph_gpu,1)
          endif

          call reorder2d_gpu( mm,nn, ptl_ph_gpu, lld, iperm_gpu )

#ifdef USE_TR_CHECK
          call reorder1d_gpu( mm, tr_save_gpu, iperm_gpu )
#endif
          call reorder1d_gpu( mm, ptl_gid_gpu, iperm_gpu )

          nn  = size( ptl_ct_gpu,2)
          lld = size( ptl_ct_gpu,1)

          if (idebug >= 1) then
           if (sml_mype.eq.0) write(*,*) 'size ptl_ct_gpu',size( ptl_ct_gpu,2), size( ptl_ct_gpu,1)
          endif

          call reorder2d_gpu( mm,nn, ptl_ct_gpu, lld, iperm_gpu )

          nn  = size( phase0_gpu,2)
          lld = size( phase0_gpu,1)
          call reorder2d_gpu( mm,nn, phase0_gpu, lld, iperm_gpu )

          if (idebug >= 1) then
           if (sml_mype.eq.0) write(*,*) 'size phase0_gpu',size( phase0_gpu,2), size( phase0_gpu,1)
          endif

          ! if we just sorted, we don't want to access through iperm
          access_thru_perm = .false.
        call t_stopf("reordering_gpu")
        endif

        deallocate(xstart_gpu)

        if (.not. perm_gpu) then
          deallocate(iperm_gpu)
        endif
      endif
     endif
     !                        
     do epc=1,sml_nrk
        call t_startf("electron_loop")

        sml_epc=epc
    
        select case(epc)
        case(1)
           dt_now=0.5_work_p*dt
        case(2)
           dt_now=dt
        end select

        if(ihybrid>1 .or. icycle >1 .or. epc > 1) diag_on=.false.

!        if(gpu_ibegin <= gpu_iend) then
        call t_startf("gpu_processing")

        call t_startf("cuda_thrd_sync1")
        istatus =  cudaThreadSynchronize()
        call t_stopf("cuda_thrd_sync1")

        if (idebug >= 1) then
!         if (sml_mype.eq.0) write(*,*) 'before kernel', sp%ptl(5)%ph, sp%ptl(5)%ct
        endif

        if (idebug >= 1)  then
!          if (sml_mype == 0) write(*,*) 'before launch kernel',icycle, sizeof(phase0_gpu),istatus
        endif

        streamid = get_gpu_streamid()

        call t_startf("pushe_kernel_gpu")
        call  pushe_kernel_gpu<<<tgrid,tblock,0,streamid>>>(istep, &
              epc,phase0_gpu,diag_on,dt_now, &
              gpu_ibegin,gpu_iend, access_thru_perm, iperm_gpu)
        call t_stopf("pushe_kernel_gpu")

        if (idebug >= 1)  then
!          if (sml_mype == 0) write(*,*) 'after launch kernel',icycle, epc, istep
        endif

        call t_stopf("gpu_processing")

        call t_startf("cpu_processing")

!        endif
! --------------------------------------
! concurrently push particles on the cpu
! --------------------------------------
        if (cpu_ibegin <= cpu_iend) then
!          if (idebug >= 1)  then
!             write(*,*) 'cpu_ibegin',cpu_ibegin, cpu_iend, sml_mype
!          endif

           call split_indices(cpu_iend-cpu_ibegin+1, sml_nthreads, i_beg, i_end)
           i_beg = i_beg + (cpu_ibegin-1)
           i_end = i_end + (cpu_ibegin-1)

           if (idebug >= 1) then
              if (sml_mype.eq.0) write(*,*) 'before kernel',sp%ptl(5)%ph, sp%ptl(5)%ct
           endif

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, X, PHI, PHI_MID, &
!$OMP          XFF, ITR, P, &
!$OMP          NEW_PHASE, OLD_PHASE, RTN  )
           do ith=1,sml_nthreads

              call t_startf("PUSHE_CPU_LOOP")
              do i=i_beg(ith),i_end(ith)
                 if(sp%ptl(i)%gid>0)  then

                    ! get proper toroidal angle index and weight
                    x=sp%ptl(i)%ph(1:2)
                    phi=sp%ptl(i)%ph(3)
                    phi_mid=(floor(phi/grid%delta_phi) + 0.5_work_p) * grid%delta_phi

                    ! get field following posision at 1/2 angle
                    call field_following_pos2(x,phi,phi_mid,xff)

                    call search_tr2(grid,xff,itr,p)

                    !remove particle or sheath calculation
                    if(itr<0) then
                       if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                          call remove_particle(sp,i,-1,ith)
                       else

                         call sheath_calculation(grid,psn,sp,i,sp%type,itr,p,ith)
                       endif
                    endif

                    sp%tr_save(i)=itr
                    sp%p_save(:,i)=p

                    select case(epc)
                    case(1)
                       sp%phase0(:,i)=sp%ptl(i)%ph
                    end select

                    call pushe_single(grid,psn,sp,i,sp%phase0(:,i),new_phase,dt_now,ith,diag_on)

                    ! check r-z boundary validity and update psi variables
                    if((new_phase(1)<eq_min_r) .or. &
                       (new_phase(1)>eq_max_r) .or. &
                       (new_phase(2)<eq_min_z) .or. &
                       (new_phase(2)>eq_max_z))       then
                       call remove_particle(sp,i,-1,ith)
!                      if (idebug >= 1) then
!                         print *, 'particle eliminated due to rz_outside :',  &
!                            i, sml_mype, sp%type, sp%ptl(i)%gid,  &
!                            new_phase(1),new_phase(2)
!                      endif
                    else
                       ! bounce 
                       if(epc==sml_nrk .and. sml_bounce/=0) then
                          old_phase(:)=sp%phase0(:,i)
                          call bounce(new_phase,old_phase,rtn)
                          if(rtn<0)  then
                             call remove_particle(sp,i,-2,ith)
                          endif
                       endif

                       !******************************************************
                       ! time advance one step
                       !******************************************************
                       sp%ptl(i)%ph= new_phase(:)

                    endif

                    if(sp%ptl(i)%ph(3)>= sml_2pi_wedge_n .or. sp%ptl(i)%ph(3)< 0D0 ) then
                       sp%ptl(i)%ph(3)=modulo(sp%ptl(i)%ph(3),sml_2pi_wedge_n)
                    endif

                 endif
              enddo !i - particle
              call t_stopf("PUSHE_CPU_LOOP")

           enddo ! ith

        endif
        call t_stopf("cpu_processing")

        call t_startf("cuda_thrd_sync2")
#ifndef NO_CUDA_SYNC2
        istatus =  cudaStreamSynchronize(streamid)
#endif
        call t_stopf("cuda_thrd_sync2")

        call t_stopf("electron_loop")
     enddo
  enddo

  ! FIXME : An allocated check would be better, but for not
  ! this is a double-check on our logic.
  if (perm_gpu) then
    deallocate(iperm_gpu)
  endif
  ! ------------------------------------------
  ! copy data from GPU device back to CPU host
  ! ------------------------------------------
  if (idebug >= 1)  then
!    if (sml_mype == 0) print*,icycle, phase0_gpu(1,1),phase0_gpu(2,1)
  endif


 if(gpu_ibegin <= gpu_iend) then
  call t_startf("push_upd_host_gpu")
  call push_update_host_gpu(  sp, &
                             psn, diag_on_input,  &
                             gpu_ibegin, gpu_iend)
  call t_stopf("push_upd_host_gpu")
  diag_on_input = diag_on
 endif
endif

call t_startf("PUSHE_SEARCH_INDEX")
call chargee_search_index(grid,psn,sp)
call t_stopf("PUSHE_SEARCH_INDEX")

  if (idebug >= 1) then
    if (sml_mype == 0) then
    write(*,9090) sml_gpu_ratio, sum(neu_weight_sum_lost)
 9090 format(' sml_gpu_ratio, sum(neu_weight_sum_lost) ',2(1x,1pe14.5))
    endif
  endif

  sml_epc=1
end subroutine pushe_gpu


