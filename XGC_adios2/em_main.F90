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
#ifdef USE_GPU
  use push_mod_gpu, only : init_push_mod_gpu, pushe_gpu
#endif
#ifdef ADIOS
  use adios_read_mod
#endif
  use xgc_interfaces
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  !
  integer :: istep,ipc,i
  PetscErrorCode::ierr
  integer :: ihybrid, epc, icycle, ncycle
  logical :: diag_on
  integer :: final_istep
  integer (8) :: idum1(0:1), idum2(0:1), idum3(0:1)
  character (len=10) :: ic(0:15)
  integer :: inode1,inode2
#ifdef ADIOS
  integer :: adios_read_method = ADIOS_READ_METHOD_BP
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
  call petsc_init( ierr ); CHKERRQ(ierr)
  call petsc_perf_init( ierr ); CHKERRQ(ierr)
  call mon_stop(PETSC_INIT_)  
  call mon_start(TOTAL_)
  call mon_start(INIT_)
  
  !ADIOS INITIALIZE
#ifdef ADIOS
  if(sml_mype==0) print *, 'call adios_init'
  call mon_start(ADIOS_INIT_)
  call adios_init('adioscfg.xml'//char(0), sml_comm, ierr)
  if(sml_mype==0) print *, 'use adios reader v2'
  call adios_read_init_method(adios_read_method, sml_comm, 'verbose=3', ierr)


  call mon_stop(ADIOS_INIT_)
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

! main loop has been restructured -- FIRST may not required
  call mon_start(FIRST_)
  if(sml_restart .and. sml_f0_grid) then
     call set_gvid0_pid_from_f0(grid%nnode)
  else
     call mon_start(SET_WEIGHTS_F_)
     call set_weights(grid,spall)
     call mon_stop(SET_WEIGHTS_F_)
  endif

  call mon_start(SHIFT_I_F_)
  call shift_sp(grid,psn,spall(1))
  call mon_stop(SHIFT_I_F_)

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

  ! first   
!  call check_point('- first run of charge')

!  call mon_start(CHARGEI_F_)
!  call chargei(grid,psn,spall(1))
!  call mon_stop(CHARGEI_F_)

!  if(sml_electron_on) then
!     call mon_start(CHARGEE_F_)
!     call chargee(grid,psn,spall(0))
!     call mon_stop(CHARGEE_F_)
!  endif
  if (.not.sml_use_ts_solver) then
  ! Solve Poisson equation
!  call check_point('- first run of poisson')

!  call mon_start(POISSON_F_)
!  call solve_poisson(grid,psn,1)
!  call mon_stop(POISSON_F_)

!  call mon_start(GET_POT_GRAD_F_)
!  call get_potential_grad(grid,psn)
!  call mon_stop(GET_POT_GRAD_F_)
  end if
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
     if(sml_electron_on) then
        idum3=spall(0:1)%num
        call mpi_allreduce(idum3,idum1,2,mpi_integer8,mpi_sum,sml_comm,ierr)
        call mpi_allreduce(idum3,idum2,2,mpi_integer8,mpi_max,sml_comm,ierr)
     else
        idum3(1)=spall(1)%num
        call mpi_allreduce(idum3(1),idum1(1),1,mpi_integer8,mpi_sum,sml_comm,ierr)
        call mpi_allreduce(idum3(1),idum2(1),1,mpi_integer8,mpi_max,sml_comm,ierr)
     endif
     call mon_stop(MAIN_LOOP_RED_) 

     if(sml_mype==0) print *, 'step,ratio,# of ion  ',istep, real(idum2(1))/(real(idum1(1))/sml_totalpe),idum1(1)
     if(sml_electron_on) then
        if(sml_mype==0) print *, 'step,ratio,# of elec ',istep, real(idum2(0))/(real(idum1(0))/sml_totalpe),idum1(0)
     endif

#ifndef ADIOS_ONLY
!pw     if( sml_pol_decomp .and. real(idum2(1))/(real(idum1(1))/real(sml_totalpe)) > sml_max_imbalance ) then
      if(.not. sml_electron_on .and.  sml_pol_decomp .and. (real(idum2(1),8) > sml_max_imbalance*real(spall(1)%min_max_num,8)) ) then
!     if (.false.) then
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

     else if (sml_electron_on) then
!pw        if( sml_pol_decomp .and. real(idum2(0))/(real(idum1(0))/real(sml_totalpe)) > sml_max_imbalance ) then
         if( sml_pol_decomp .and. (real(idum2(0),8) > sml_max_imbalance*real(spall(0)%min_max_num,8)) ) then
!        if (.false.) then
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
        else
           pol_decomp_new=.false.
        endif
     else
        pol_decomp_new=.false.
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
     
     ! particle pushing RK2 loop ---------------------------------------
     do ipc=1, sml_nrk
        call mon_start(IPC_LOOP_)
        sml_ipc=ipc
        
        !* Obtain ion charge density
        call mon_start(CHARGEI_)
        call chargei(grid,psn,spall(1))  ! weight updated here
        call chargee_hyb_mpisum(grid,psn)
        call mon_stop(CHARGEI_)        

        call push_fluid(grid,psn,ipc)

        call mon_start(GET_POT_GRAD_)
        call get_potential_grad(grid,psn)
        call mon_stop(GET_POT_GRAD_)
        
        call t_startf("GET_POT_GRAD_EPARA")
        call get_potential_grad_epara(grid,psn)
        call t_stopf("GET_POT_GRAD_EPARA")

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
  if (sml_use_ts_solver) then
     call ts_delete(psn%ts,ierr);CHKERRQ(ierr) 
  else
     call KSPDestroy(psn%kspmass,ierr)
     call delete_solver( psn%solver00, ierr )
  end if
#ifdef ADIOS
  call adios_finalize(sml_mype,ierr)
  call adios_read_finalize_method (adios_read_method, ierr)
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
  
!1000 format (8(e10.3))
!1001 format (8(f9.3,'%'))
end program xgc1_3

#include "main_extra.F90"



