!********************************************************
!! XGC1
!!
!! @Version 3.0  2/20/2013
!! @Authors:
!!  S. Ku
!!  P. Worley, 
!!  E. D'Azevedo, 
!!  M. Adams,
!!  G. Park,
!!  E. Yoon, 
!!  S.H. Koh
!!  J.H. Seo
!!  S. Klasky
!!  J. Lang
!!  EPSI Team
!********************************************************
program xgc1_3
  use sml_module
  use ptl_module
  use grid_class
  use psn_class
  use diag_module
  use smooth_module
  use pol_decomp_module
  use src_module
  use random_xgc
  use perf_monitor
  use neu_module
  use f0_module
  use col_module
  use lbal_module
#ifdef USE_GPU
  use push_mod_gpu, only : init_push_mod_gpu, pushe_gpu
#endif
#ifdef XGC1_EM
#error "Must NOT use -DXGC1_EM when building xgc-es"
#endif
#ifdef ADIOS
  use adios_comm_module
#endif
#ifdef ADIOS2
  use adios2_comm_module
#endif
#ifdef XGC_COUPLING_CORE_EDGE
  use coupling_core_edge
!  use new_coupling_xgc
#endif
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  !
  integer :: istep,ierr,ipc,i
  integer :: ihybrid, epc, icycle, ncycle
  logical :: diag_on
  integer :: final_istep
  real (8) :: local_cost(0:3), max_cost(0:3), total_cost(0:3)
  real (8) :: cumul_cost, cumul_max_cost, cumul_total_cost
  real (8) :: f0_cost, f0_max_cost, f0_total_cost
  real (8) :: ion_max_cnt, ion_total_cnt
  real (8) :: elec_max_cnt, elec_total_cnt
  real (8) :: elec_rebal_trigger, ion_rebal_trigger
  character (len=10) :: ic(0:15)
  integer :: inode1,inode2
#ifdef XGC_COUPLING_CORE_EDGE
  integer :: nn
#endif
  interface
     subroutine mon_stop(i,est_final_step)
       integer, intent(in) :: i
       integer, optional, intent(inout) :: est_final_step
     end subroutine mon_stop
  end interface

  ! MPI initialize 
  call my_mpi_init 

  call init_perf_monitor()

  if(sml_mype==0) print *, 'call petsc_init'
  call mon_start(PETSC_INIT_)
  call petsc_init( ierr )
  call petsc_perf_init( ierr )
  call mon_stop(PETSC_INIT_)
  
  call mon_start(TOTAL_)
  call mon_start(INIT_)
  
  !ADIOS INITIALIZE
#ifdef ADIOS
  if(sml_mype==0) print *, 'call adios_init'
  call t_startf("ADIOS_INIT")
  call adios_comm_init()
  call t_stopf("ADIOS_INIT")
#endif

#ifdef ADIOS2
  call adios2_comm_init('adios2cfg.xml')
#endif
  
  call mon_start(SETUP_)
  call setup(grid,psn,spall)
  call mon_stop(SETUP_)

  !GPU initialize
#ifdef USE_GPU
  call check_point('Init GPU')
  call mon_start(INIT_PUSHMOD_GPU_)
  call init_push_mod_gpu( grid )
  call mon_stop(INIT_PUSHMOD_GPU_)
#endif

  call mon_stop(INIT_)
  
  if (mon_flush_freq > 0) call flush_perf_monitor(-1)

  ! main loop --------------------------
  call check_point('main loop started')

! load balancing cost tracking initialization
  cumul_cost    =  0.D0
  f0_ratio      = -1.D0
  f0_node_cost_updated = .false.

! main loop has been restructured -- FIRST may not required
  call mon_start(FIRST_)
  if(sml_restart .and. sml_f0_grid) then
     call set_gvid0_pid_from_f0(grid%nnode)
  else
     ! f0_ratio == -1.D0, so only particles are load balanced
     ! here even if sml_f0_grid_load_balance is .true.
     call mon_start(SET_WEIGHTS_F_)
     call set_weights(grid,spall)
     call mon_stop(SET_WEIGHTS_F_)
  endif

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#warning 'Load balancing is turned off'
  call mon_start(SHIFT_I_F_)
  call shift_sp(grid,psn,spall(1))
  call mon_stop(SHIFT_I_F_)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(sml_electron_on) then
    call mon_start(SHIFT_E_F_)
    call shift_sp(grid,psn,spall(0))
    call mon_stop(SHIFT_E_F_)
  endif

  if(sml_f0_grid) then    
     if(.not. sml_restart) then
        inode1=gvid0_pid(sml_plane_mype)
        inode2=gvid0_pid(sml_plane_mype+1)-1
        !print *, sml_mype,']]]',inode1,inode2
        call f0_initialize(grid,inode1,inode2,0,f0_nmu)
        if (sml_plane_index==0) then
#ifdef ADIOS
          if (sml_mype==0) print *,"dumping f0_grid"
          call dump_f0_grid(grid)
#endif
        endif
     else
        call f0_init_rest(grid)
     endif
  endif

  if(sml_f0_grid) then
    ! set charge density from restart data
    call t_startf("CHARGEI_F0")
    call chargei_f0(grid,psn)
    call t_stopf("CHARGEI_F0")

    if (sml_electron_on) then
      call t_startf("CHARGEE_F0")
      call chargee_f0(grid,psn)
      call t_stopf("CHARGEE_F0")
    endif
  endif

  !call mon_start(DIAGNOSIS_)
  !call diagnosis(0,1,grid,psn,spall)  ! before first time step: initial conditions. 
  !call mon_stop(DIAGNOSIS_)

  call mon_stop(FIRST_)

  if (mon_flush_freq > 0) call flush_perf_monitor(0)

  call check_point('- do loop start')
  final_istep = sml_mstep+1
  do istep=1,sml_mstep
     call mon_start(MAIN_LOOP_)
     sml_istep=istep 
     sml_gstep=sml_gstep+1
     
     ! load balancing part ---------------------------------------------
     ! getting total number of particle and load imbalance
     call mon_start(MAIN_LOOP_RED_)
     if(sml_f0_grid) then
        f0_cost = sum(f0_node_cost)
     else
        f0_cost = 0.0D0
     endif
     local_cost(3)=f0_cost
     local_cost(2)=cumul_cost
     local_cost(1)=real(spall(1)%num,8)
     if(sml_electron_on) then
        local_cost(0)=real(spall(0)%num,8)
        call mpi_allreduce(local_cost(0),total_cost(0),4,mpi_real8,mpi_sum,sml_comm,ierr)
        call mpi_allreduce(local_cost(0),max_cost(0),4,mpi_real8,mpi_max,sml_comm,ierr)
     else
        call mpi_allreduce(local_cost(1),total_cost(1),3,mpi_real8,mpi_sum,sml_comm,ierr)
        call mpi_allreduce(local_cost(1),max_cost(1),3,mpi_real8,mpi_max,sml_comm,ierr)
     endif
     f0_max_cost       = max_cost(3)
     f0_total_cost     = total_cost(3)
     ! cumul cost estimates are meaningless except when f0_node_cost_updated == .true.
     cumul_max_cost    = max_cost(2)
     cumul_total_cost  = total_cost(2)
     ion_max_cnt       = max_cost(1)
     ion_total_cnt     = total_cost(1)
     call mon_stop(MAIN_LOOP_RED_) 

     ion_rebal_trigger = sml_max_imbalance*((real(spall(1)%min_max_num,8)/ion_total_cnt)*sml_totalpe)
     if(sml_mype==0) write(6,1110) istep, ion_rebal_trigger, &
                       (ion_max_cnt/ion_total_cnt)*sml_totalpe,int(ion_total_cnt,8)
1110 format('step,trigger,ratio,# of ion  ',I5,1x,F10.4,1x,F10.4,1X,I14)
     if(sml_electron_on) then
        elec_max_cnt       = max_cost(0)
        elec_total_cnt     = total_cost(0)
        elec_rebal_trigger = sml_max_imbalance*((real(spall(0)%min_max_num,8)/elec_total_cnt)*sml_totalpe)
        if(sml_mype==0) write(6,1111) istep, elec_rebal_trigger, &
                          (elec_max_cnt/elec_total_cnt)*sml_totalpe,int(elec_total_cnt,8)
1111 format('step,trigger,ratio,# of elec ',I5,1x,F10.4,1x,F10.4,1X,I14)
     endif

     if (f0_total_cost > 0.0D0) then
        f0_ratio = (f0_max_cost/f0_total_cost)*sml_totalpe
     else
        f0_ratio = -1.D0
     endif

     if (f0_node_cost_updated) then
        if(sml_mype==0) write(6,1112) istep, f0_ratio, f0_max_cost, cumul_max_cost
1112 format('step,f0(ratio,max),total(max) ',I5,1X,F10.4,1X,F10.4,1X,F10.4)
     else
        if(sml_mype==0) write(6,1113) istep, f0_ratio, f0_max_cost
1113 format('step,f0(ratio,max)',I5,1X,F10.4,1X,F10.4)
     endif

     pol_decomp_new=.false.

#ifndef ADIOS_ONLY
     if ( sml_pol_decomp ) then
        ! a) If sml_f0_grid_lbal_period > 0, then any step after a collision operation for
        !    which mod(istep,sml_f0_grid_lbal_period) == 1 and f0_ratio > 0 will 
        !    have the collision cost rebalanced, subject to a limit on the acceptable
        !    particle load imbalance. (The step after the first collision step also 
        !    invokes load balancing when f0_ratio > 0.)
        ! b) Before the first collision step f0_node_cost and f0_ratio are <= 0 and 
        !    collision cost is not considered during load balancing, only particles.
        ! c) f0_node_cost is updated only after a collision computation 
        !    (f0_node_cost_updated == .true.), but the f0 node load imbalance will also change 
        !    whenever the partition is recalculated, e.g. because of particle load imbalance.
        !    When sml_f0_grid_lbal_period > 0 and f0_ratio is > 0, f0 cost will be 
        !    load balanced whenever load balancing is triggered, for whatever reason.

        if (istep > sml_f_source_period) then
           if (((istep == (sml_f_source_period+1)) .or. &
               (mod(istep,sml_f0_grid_lbal_period) == 1)) &
               .and. f0_grid_load_balance) then
              ! update history and update particle load imbalance bound when collision 
              ! load rebalancing is enabled
              call update_ptl_constraint(cumul_max_cost, f0_grid_ptl_imbal)
           endif
        endif

        if (.not. sml_electron_on) then
           if ((ion_max_cnt > sml_max_imbalance*real(spall(1)%min_max_num,8)) .or. &
               (((mod(istep,sml_f0_grid_lbal_period) == 1) .or. (istep == (sml_f_source_period+1))) &
                .and. f0_grid_load_balance .and. (istep > sml_f_source_period))) then

              ! load balance
              call mon_start(SET_WEIGHTS_)
              gvid0_pid_old=gvid0_pid
              call set_weights(grid,spall)   ! enabled 
              call mon_stop(SET_WEIGHTS_)

              if (sml_f0_grid) then
                 call t_startf("F0_REDISTRIBUTE")
                 call f0_redistribute(grid,0,f0_nmu)
                 call t_stopf("F0_REDISTRIBUTE")
              endif

              pol_decomp_new=.true.

           endif
        else
           ! while rebalancing focuses on electrons and on f0 cost, also need to check whether 
           ! ions exceed load imbalance limit so that memory issues do not occur
           if ((elec_max_cnt > sml_max_imbalance*real(spall(0)%min_max_num,8)) .or. &
               (ion_max_cnt > sml_max_imbalance*real(spall(1)%min_max_num,8))  .or. &
               (((mod(istep,sml_f0_grid_lbal_period) == 1) .or. (istep == (sml_f_source_period+1))) &
                .and. f0_grid_load_balance .and. (istep > sml_f_source_period))) then

              ! load balance
              call mon_start(SET_WEIGHTS_)
              gvid0_pid_old=gvid0_pid
              call set_weights(grid,spall)   ! enabled 
              call mon_stop(SET_WEIGHTS_)

              if (sml_f0_grid) then
                 call t_startf("F0_REDISTRIBUTE")
                 call f0_redistribute(grid,0,f0_nmu)
                 call t_stopf("F0_REDISTRIBUTE")
              endif

              pol_decomp_new=.true.

           endif
        endif
     endif

     if (pol_decomp_new) then
        call t_startf("SHIFT_I_R")
        call shift_sp(grid,psn,spall(1))
        call t_stopf("SHIFT_I_R")

        if(sml_electron_on) then
           call t_startf("SHIFT_E_R")
           call shift_sp(grid,psn,spall(0))
           call t_stopf("SHIFT_E_R")
        endif

        pol_decomp_new=.false.
     endif
     
     if (f0_node_cost_updated) then
        ! restart cost tracking after collision step
        cumul_cost = 0.D0
        f0_node_cost_updated = .false.
     endif

     cumul_cost = cumul_cost - mpi_wtime()

     ! particle pushing RK2 loop ---------------------------------------
     do ipc=1, sml_nrk
        call mon_start(IPC_LOOP_)
        sml_ipc=ipc
        
        !* Obtain ion charge density
        call mon_start(CHARGEI_)
        call chargei(grid,psn,spall(1))  ! weight updated here
        call mon_stop(CHARGEI_)        
        

        !* Obtain electron charge density

        if(sml_electron_on) then           
           call mon_start(ELECTRONS_)
           
           ! solve poisson equation with n_i(t+dt) 
           ! and n_e_adiabatic and store it dpot
           call mon_start(POISSON_E1_)
           call solve_poisson(grid,psn,2)  ! ignore electron response (non-adiabatic)
           call get_dpot_ff(grid,psn,ipc)  ! small version of get_potential_grad
           call mon_stop(POISSON_E1_)
           
           
           ! Obtain electron charge density 
           ! using adiabatic response
           call mon_start (CHARGEE_)
           call chargee(grid,psn,spall(0)) ! weight updated here
           call mon_stop (CHARGEE_)  
           
           
           ! call mon_start(HYBRID_LOOP_)
           ! do loop for iteration
           do ihybrid=2, sml_nhybrid
              
              ! Solve Poisson equation for updating ddpotdt 
              call mon_start(POISSON_E2_)
              call solve_poisson(grid,psn,1)
              call get_dpot_ff(grid,psn,0)
              call mon_stop(POISSON_E2_)
              
              
              ! Obtain electron charge density             
              call mon_start(CHARGEE_)
              call chargee(grid,psn,spall(0))
              call mon_stop(CHARGEE_)  
           enddo
           
           ! logical sheath calculation
           if(sml_sheath_mode /= 0 .and. sml_sheath_adjust .and. sml_ipc==2) then
             ! Adjust logical sheath
             call sheath_adjust(grid,psn)
           endif 
           call mon_stop(ELECTRONS_)
        endif
        
        ! Solve Poisson equation - final field for particle push --------
        call mon_start(POISSON_)
        call solve_poisson(grid,psn,1)
        if(sml_electron_on) call save_dpot(grid,psn,ipc)
        call mon_stop(POISSON_)

#ifdef XGC_COUPLING_CORE_EDGE_FIELD
        nn=grid%nnode
        if(sml_plane_mype.eq.0)then
          call cce_send_field(psn%dpot(:,0),psn%dpot(:,1),psn%pot0(:),sml_plane_mype)
          call cce_receive_field(sml_plane_mype)
          call cce_process_field(psn%dpot(:,0),psn%dpot(:,1),psn%pot0(:),sml_plane_mype)
        endif
        !if(cce_side.eq.0)then 
          call mpi_bcast(psn%dpot(:,0),nn,mpi_real8,0,sml_plane_comm,ierr)
          call mpi_bcast(psn%dpot(:,1),nn,mpi_real8,0,sml_plane_comm,ierr) 
          call mpi_bcast(psn%pot0(:)  ,nn,mpi_real8,0,sml_plane_comm,ierr)
        !endif
#endif
        
        call mon_start(GET_POT_GRAD_)
        call get_potential_grad(grid,psn)
        call mon_stop(GET_POT_GRAD_)  

        if(sml_electron_on) then
        !##############################################################
        ! non-blocking gather_field_info can start here
        ! (+ GPU pushe)
        !##############################################################
        endif

        ! Time advance of ion phase variables
        call determine_diag_on(istep,ipc,diag_on)
        if(diag_3d_more .and. diag_on) then
           call t_startf("diag_3d_additional")
           call diag_3d_additional(grid,psn,spall(1))
           call t_stopf("diag_3d_additional")
        endif

        ! Push electron particles ----------------------------------
        ! electron subcycling routine
        if(sml_electron_on) then
           ihybrid=1 ! for legacy
           call mon_start(SOL_ELEC_PHASE_)
           call save_or_load_electron_phase(spall(0),ipc,ihybrid) ! Also f0 is restored.
           call mon_stop(SOL_ELEC_PHASE_)
              
           select case(ipc)
           case(1)
              ncycle=sml_ncycle_half
           case(2)
              ncycle=sml_ncycle_half*2
           end select
              
           !##################################################################
           !electron push with sub-cycling
           call mon_start(GAT_FIELD_INFO_)
           call gather_field_info(grid,psn) ! gather other plane E-field info
           call mon_stop(GAT_FIELD_INFO_)
           ! gather_field_info should be block here
           !#################################################################

           call mon_start(PUSHE_)
#ifdef USE_GPU
           call pushe_gpu(istep,ihybrid,ncycle,grid,psn,spall(0),diag_on)
#else
           call pushe(istep,ihybrid,ncycle,grid,psn,spall(0),diag_on)
#endif
           call mon_stop(PUSHE_)
          
        endif
        

        ! Push ion particles ------------------------------------------
        ! Push-ion is after push-electron for GPU parallelism
        call mon_start(PUSH_I_)
        call determine_diag_on(istep,ipc,diag_on)
        call push(istep,ipc,grid,psn,spall(1),spall(1)%phase0,spall(1)%ptl,diag_on)
        call mon_stop (PUSH_I_)

        ! electron push should be finished before this point ---------
        
        call mon_start(DIAGNOSIS_)
        call diagnosis(istep,ipc,grid,psn,spall)  ! most of diagnosis output happens here. Some are in f_source 
        call mon_stop(DIAGNOSIS_)
        

        ! Move ions to their domain
        call mon_start(SHIFT_I_)
        call mon_start(MEM_CLEAN_I_)
        call memory_cleaning_simple(spall) !## ion and electron sepratately?
        call mon_stop(MEM_CLEAN_I_)
        call shift_sp(grid,psn,spall(1))  
        call mon_stop(SHIFT_I_)



        if(sml_electron_on) then
           call mon_start(SHIFT_E_)
           call shift_sp(grid,psn,spall(0))
           call mon_stop(SHIFT_E_)
        endif

        call mon_stop(IPC_LOOP_)
     enddo  !end ion ipc-loop 
     


     ! update f0 , f  - after shift 
     if(sml_f0_grid) then
        call mon_start(F0_GRID_)
        ! need to search triangle index -- redundant
        call t_startf("F0_CHARGEI_SEARCH_INDEX")
        call chargei_search_index(grid,psn,spall(1))
        call t_stopf("F0_CHARGEI_SEARCH_INDEX")

        if (sml_electron_on) then
           call t_startf("F0_CHARGEE_SEARCH_INDEX")
           call chargee_search_index(grid,psn,spall(0))
           call t_stopf("F0_CHARGEE_SEARCH_INDEX")
        endif
        
        !call update_f0_sp(grid,spall(1))      ! moved to f_source  --> after update_w_ion
        !if(sml_electron_on) call update_f0_sp(grid,spall(0))

        ! coulomb collision, neutral collision, and source - update f0_grid
        !rh Check of time step is done in f_source now!
        call f_source(grid,psn, spall)


        ! get toroidal average of f0g
        ! this routine can be located between update_f0_sp and f_source
        ! --> this will give slight axisymmetry for dt
        if(sml_symmetric_f0g) then
           call t_startf("SYMMETRIC_F0G")
           call symmetric_f0g  ! corresponding w1 will be calculated in charge
           call t_stopf("SYMMETRIC_F0G")
        endif


        ! update charge density 
        call t_startf("CHARGEI_F0")
        call chargei_f0(grid,psn)
        call t_stopf("CHARGEI_F0")

        if (sml_electron_on) then
           call t_startf("CHARGEE_F0")
           call chargee_f0(grid,psn)
           call t_stopf("CHARGEE_F0")
        endif

        if(mod(istep,sml_f_source_period)==0) then 
           call t_startf("RESET_F0")
           call reset_f0_f
           call t_stopf("RESET_F0")
        endif

        call mon_stop(F0_GRID_)
     endif
     
     ! coulomb collision, neutral collision, source for full-f method
     if(.not. sml_f0_grid) then
        !Calculate f0 for ion-ion collision 
        call mon_start(COLLISION_)
        call collision(istep,spall(1))
#ifdef VPIC_COL
        if(mod(istep,col_period) .eq. 0 .and. col_mode .eq. 3) call vpic_col(grid,spall) !i-i,i-e,e-i,e-e
#endif
        call mon_stop(COLLISION_)
        
        ! ion/electron-neutral collision
        if(sml_neutral.and.(.not. sml_deltaf)) then
           call mon_start(NEUTRAL_)
           call neutral_col(istep,grid,psn,spall(1),spall(1))
           call mon_stop(NEUTRAL_)
        endif
        
        if(sml_source .and. mod(sml_gstep, src_period)==0) then
           call mon_start(SOURCE_)
           call source(sml_gstep, spall(1))
           call mon_stop(SOURCE_)
        endif
     endif


     ! update current time
     sml_time=sml_time+sml_dt

     cumul_cost = cumul_cost + mpi_wtime()
#endif

     ! periodic restart dump
#if defined(ADIOS) 
     if ((mod(sml_gstep,sml_restart_write_period)==0) .or. (istep == final_istep)) then
        call mon_start(RESTART_WRITE_)

        call mon_start(MEM_CLEAN_RSTRT_)       
        call shift_sp(grid,psn,spall(1))
        if(sml_electron_on) call shift_sp(grid,psn,spall(0))
        call memory_cleaning_simple(spall)
        call mon_stop(MEM_CLEAN_RSTRT_)

        call restart_write(spall,grid, psn)
        if(sml_mype==0 .and. sml_neutral) call background_edensity0_output(grid,psn) ! for neutral collision
        call mon_stop(RESTART_WRITE_)
     endif
#endif
     call mon_stop(MAIN_LOOP_,final_istep)
     
     ! Timer information
     if (mon_flush_freq > 0) then
        if (mod(istep,mon_flush_freq) .eq. 0) &
             call flush_perf_monitor(istep)
     endif

     if (istep >= final_istep) exit
  enddo !end of main loop
  
  call mon_start(FINALIZE_)
  
  if(sml_mype==0) then
     close(30)
     close(31)
     close(71)
     close(22)
  endif

  if(sml_turb_poisson) call delete_solver( psn%solverH, ierr )  
  if(.NOT. sml_use_simple00 ) call delete_solver( psn%solver00, ierr )

#ifdef ADIOS
  call adios_comm_finalize()
#endif
#ifdef ADIOS2
  call adios2_comm_finalize()
#endif

  ! free memories
  call smooth_pol_delete(smooth0L)
  call smooth_pol_delete(smoothH)
  call smooth_pol_delete(smoothdiag)
  call smooth_r_delete(smooth_r1)
  ! free particle memories
  ! ----??

  ! free solver memories
  ! ----??

  !
  call finalize_pol_decomp

  call mon_stop(FINALIZE_)
  call mon_stop(TOTAL_)

  call finish_perf_monitor(ierr)
  ! MPI_finalize
  call my_mpi_finalize
  
1000 format (8(e10.3))
!1001 format (8(f9.3,'%'))
end program xgc1_3


#include "main_extra.F90"
