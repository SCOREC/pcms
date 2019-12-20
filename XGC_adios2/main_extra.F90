


!routines for debug-----------
!! NaN check
subroutine nan_check(tag,sp)
  use ptl_module
  use sml_module
  implicit none
  type(species_type) :: sp
  integer :: i,j
  real (kind=8) :: var
  character (len=30) :: tag

  do i=1, sp%num
     if(sp%ptl(i)%gid>0) then
        do j=1, ptl_nphase
           var=sp%ptl(i)%ph(j)
           if(var > 1D0 .or. var < 2D0 )  then
           else
              print *,'NaN found in',tag, sml_mype,i,j,' : ',sp%ptl(i)%ph(:) 
              stop
           endif
        enddo
     endif
  enddo
end subroutine nan_check

!! subroutine for mu<0 check, only for debug
subroutine minus_mu_check(tag,sp)
  use ptl_module
  use sml_module
  implicit none
  type(species_type) :: sp
  integer :: i,j
  real (kind=8) :: var
  character (len=30) :: tag
  do i=1, sp%num
     if(sp%ptl(i)%gid>0) then
        var=sp%ptl(i)%ct(pim)
        if(var <0D0 )  then
           print *,'minus mu found ',tag, sml_mype,i,' : ',sp%ptl(i)%ph(:) 
           print *,'-----------------'
           stop
        endif
     endif
  enddo
end subroutine 


!-------------------------------------------------------------------------------------
subroutine mon_start(i)
  use perf_monitor
  use sml_module
  use rem_module
#if defined(CAM_TIMERS)
  use perf_mod, only: t_barrierf, t_startf, t_adj_detailf
#endif
  implicit none
  include 'mpif.h'

  integer, intent(in) :: i

  integer ierr
  logical exist

  namelist /rem_param/ rem_walltime, rem_walltime_src

#if defined(CAM_TIMERS)
  if (i <= mon_N2) then
    if (mon_sync(i)) then
      call t_barrierf(mon_str(i+mon_NX), sml_comm)
    endif

    call t_startf(mon_str(i))
    call t_adj_detailf(mon_adj(i))
  endif
#else
  if ((i <= mon_N) .and. (mon_sync(i))) call mpi_barrier(sml_comm,ierr)
#endif

#if !defined(NO_PETSC)
  if (i <= mon_N) call PetscLogEventBegin( event(i), ierr )
#else
  if (i <= mon_N) call cpu_time( mon_time(i) )
#endif
#if defined(XGC_MSG)
  call mpi_barrier(sml_comm,ierr)
  if(sml_mype==0) print *, mon_str(i)//' START'
#endif

! Logic for estimating end of execution time
if (i == RESTART_WRITE_) then
  if (sml_mype == 0) then
    ! if estimating remaining time, capture restart write timestamp
    if (rem_estimate) then
      rem_restart_write_start = mpi_wtime()
    endif
  endif
endif

if (i == MAIN_LOOP_) then
  ! reset restart write flag prior to starting loop iteration
  rem_restart = .false.
!
  if (sml_mype == 0) then
    ! if first time, check that a start time was recorded and
    ! whether the file Walltime.Remaining exists and contains 
    ! an estimate of how much time remained when the job started
    if (rem_walltime_src == -1) then
      if (rem_starttime >= 0.0) then
        ! look for remaining walltime estimate in input file first
        open(unit=14,file='input',action='read')
        read(14,nml=rem_param,iostat=ierr)
        close(14)
        if (ierr == 0) then
          rem_walltime_src = 1
        else
          rem_walltime = -1
        endif
        ! if a Walltime.Remaining file exists, use this estimate
        ! instead
        inquire(FILE="Walltime.Remaining",EXIST=exist)
        if (exist) then
          open(unit=14,file='Walltime.Remaining',action='read')
          read(14,nml=rem_param,iostat=ierr)
          close(14)
          if (ierr == 0) then
            if (rem_walltime_src == -1) then
              ! if 'provenance' not indicated, which should never 
              ! occur, assume that set in pre_aprun script
              rem_walltime_src = 2
            endif
          endif
        endif
        if (rem_walltime <= 0) then
          rem_walltime_src = 0
          rem_estimate = .false.
        else
          rem_estimate = .true.
        endif
      else
        rem_walltime_src = 0
        rem_estimate = .false.
      endif
      call mpi_bcast(rem_walltime_src, 1, mpi_integer, 0, sml_comm, ierr)
    endif

    ! capture step loop timestamp
    if (rem_estimate) then
      rem_step_loop_start = mpi_wtime()
    endif
  else
    if (rem_walltime_src == -1) then
      call mpi_bcast(rem_walltime_src, 1, mpi_integer, 0, sml_comm, ierr)
      if (rem_walltime_src > 0) then
        rem_estimate = .true.
      endif
    endif
  endif
endif

end subroutine mon_start

!-------------------------------------------------------------------------------------
subroutine mon_stop(i,est_final_step)
  use perf_monitor
  use sml_module
  use rem_module
#if defined(CAM_TIMERS)
  use perf_mod, only: t_stopf, t_adj_detailf
#endif
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i
  integer, optional, intent(inout) :: est_final_step

  integer :: steps_remaining, save_final_step_est
  real (kind=8) :: restart_write_stop
  real (kind=8) :: step_loop_stop
  real (kind=8) :: step_loop_time
  real (kind=8) :: time_remaining
  real (kind=8) :: step_loop_avg

  integer ierr
  logical exist

#if defined(NO_PETSC)
  real (kind=8) :: t2
#endif
  
  namelist /rem_param/ rem_walltime, rem_walltime_src

#if !defined(NO_PETSC)
  if (i <= mon_N) call PetscLogEventEnd( event(i), ierr )
#else
  if (i <= mon_N) then
    call cpu_time(t2)
    mon_sum(i)=mon_sum(i)+t2-mon_time(i)
  endif
#endif

#if defined(CAM_TIMERS)
  if (i <= mon_N2) then
    call t_adj_detailf(-mon_adj(i))
    call t_stopf(mon_str(i))
  endif
#endif

! more logic for estimating how much time left, and to act accordingly
if (i == RESTART_WRITE_) then
  rem_restart = .true.
  if (sml_mype == 0) then
    ! if estimating remaining time, calculate time spent in latest restart write
    ! and update statistics
    if (rem_estimate) then
      restart_write_stop = mpi_wtime()
      rem_restart_write = restart_write_stop - rem_restart_write_start
      if (rem_restart_write_max < rem_restart_write) then
        rem_restart_write_max = rem_restart_write
      endif
    endif
  endif
endif

if (i == MAIN_LOOP_) then
  ! abort remaining time estimation if est_final_step not included in call
  if (.not. present(est_final_step)) then
    rem_estimate = .false.
  endif
  if (sml_mype == 0) then
    ! if estimating remaining time, calculate time spent in latest step loop
    ! (minus a restart write) and update statistics
    if (rem_estimate) then
      step_loop_stop = mpi_wtime()
      time_remaining = dble(rem_walltime) - (step_loop_stop - rem_starttime)
      step_loop_time = step_loop_stop - rem_step_loop_start

      if (rem_step_loop_cnt == -1) then
        ! if first time, looking for estimate of halfway point (in step count)
        steps_remaining = int(time_remaining / (2*step_loop_time)) - 1
        if (steps_remaining < 1) then
          ! since first time, try to go at least one more step
          rem_final_step_est = sml_istep + 1
        else
          rem_final_step_est = sml_istep + steps_remaining
        endif

        ! send target to other processes, then initialize counter, max, and running sum
        call mpi_bcast(rem_final_step_est, 1, mpi_integer, 0, sml_comm, ierr)
        rem_step_loop_cnt = 0
        rem_step_loop_max = 0.0
        rem_step_loop_sum = 0.0
        est_final_step = rem_final_step_est
        rem_final_step_upd = .false.
        ! mon_flush_freq update needs to reflect looking for halfway point;
        ! not updating after this - fixed flush frequency is more important than 
        ! total number of flushes
        if (mon_flush_count > 0) then
          mon_flush_freq  = max(1,(2*est_final_step)/mon_flush_count)
        endif
      else
        ! if not first time, then ...

        ! after restart write or when have reached most recent target
        if ((rem_restart) .or. (est_final_step == sml_istep)) then
          ! check whether now have a better estimate of the start time
          if (rem_walltime_src < 3) then
            inquire(FILE="Walltime.Remaining",EXIST=exist)
            if (exist) then
              open(unit=14,file='Walltime.Remaining',action='read')
              read(14,nml=rem_param,iostat=ierr)
              close(14)
              if (ierr == 0) then
                ! recompute time_remaining
                time_remaining = dble(rem_walltime) - (step_loop_stop - rem_starttime)
                ! stop checking after this time
                rem_walltime_src = 3
              endif
            endif
          endif
        endif

        !save current estimate
        save_final_step_est = rem_final_step_est

        ! and calculate new (conservative) estimate of end time.
        ! note 1: removing restart write cost from loop cost when loop included a restart, 
        !  and using max restart cost when need to estimate number of loop iterations 
        !  until run out of time. 
        ! note 2: every time end-of-run (or half way) estimate is reached a
        !  restart write is forced (if adios enabled)
        ! note 3: every restart involves an update of the estimate, so will not have 
        !  multiple restarts per estimate update, and halfway estimate (probably) 
        !  won't be unrealistic because it will be overridden once do first restart
        ! note 4: after reaching midpoint, target real end time (with 'conservative'
        !  end-of-time condition). 
        ! note 5: if taking profile data checkpoints, will update estimated time
        !  for every checkpoint as well
        ! note 6: monitoring as go, and looking to send an unexpected termination message
        !  if some loops take unexpectedly long (still to do)
        if ((rem_restart) .and. (rem_restart_write >= 0.0)) then
          step_loop_time = step_loop_time - rem_restart_write
        endif
        if (rem_step_loop_max < step_loop_time) then
          rem_step_loop_max = step_loop_time
        endif
        rem_step_loop_sum = rem_step_loop_sum + step_loop_time
        rem_step_loop_cnt = rem_step_loop_cnt + 1
        step_loop_avg = (rem_step_loop_sum / rem_step_loop_cnt)
        steps_remaining = &
          int ((time_remaining - 2*(rem_step_loop_max + rem_restart_write_max)) &
               / max(step_loop_avg,step_loop_time)); 
        if (steps_remaining < 1) then
          rem_final_step_est = sml_istep
        else
          rem_final_step_est = sml_istep + steps_remaining
        endif
!
        ! after restart write or when have reached most recent target, send new 
        ! target to other processes, then reset counter and running sum
        if ((rem_restart) .or. (est_final_step == sml_istep)) then
          call mpi_bcast(rem_final_step_est, 1, mpi_integer, 0, sml_comm, ierr)
          rem_step_loop_cnt = 0
          rem_step_loop_sum = 0.0
          est_final_step = rem_final_step_est
          rem_final_step_upd = .false.
          rem_restart_any = .true.
        ! not updating mon_flush_freq - fixed flush frequency is more important 
!pw       if (mon_flush_count > 0) then
!pw         mon_flush_freq  = max(1,est_final_step/mon_flush_count)
!pw       endif
!pw      else if (steps_remaining < 1) then
!pw       call mpi_send (to everyone), possibly via a tree algorithm, even with irecv/probe
!pw       still force a restart somehow?
        endif

        if (rem_final_step_upd) then
          ! return saved estimate if broadcast during last flush_perf_monitor
          ! and if current estimate not broadcast because of restart write or from reaching 
          ! current target
          est_final_step = save_final_step_est
          rem_final_step_upd = .false. 
        endif

      endif
    endif
  else
    ! if estimating remaining time, look for updates from process 0
    if (rem_estimate) then
      if (rem_step_loop_cnt == -1) then
        ! if first time, looking for estimate of halfway point (in step count)
        call mpi_bcast(rem_final_step_est, 1, mpi_integer, 0, sml_comm, ierr)
        rem_step_loop_cnt = 0
        est_final_step = rem_final_step_est
        ! mon_flush_freq update needs to reflect looking for halfway point
        if (mon_flush_count > 0) then
          mon_flush_freq  = max(1,(2*est_final_step)/mon_flush_count)
        endif
      else if ((rem_restart) .or. (est_final_step == sml_istep)) then
        ! if after write_restart or found targeted step, get update from process 0
        call mpi_bcast(rem_final_step_est, 1, mpi_integer, 0, sml_comm, ierr)
        est_final_step = rem_final_step_est
        rem_restart_any = .true.
        ! not updating mon_flush_freq - fixed flush frequency is more important 
!pw     if (mon_flush_count > 0) then
!pw       mon_flush_freq  = max(1,est_final_step/mon_flush_count)
!pw     endif
      else if (rem_final_step_upd) then
        ! return value updated during most recent call to flush_perf_monitor
        est_final_step = rem_final_step_est
      endif
      ! in all cases, mark that an update is not available
      rem_final_step_upd = .false. 
    endif
  endif
endif

#if defined(XGC_MSG)
  call mpi_barrier(sml_comm,ierr)
  if(sml_mype==0) print *, mon_str(i)//' END'
#endif

end subroutine mon_stop

subroutine range_check(spall)
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i,isp,rz_outside
  real (kind=8) :: r,z

  do isp=ptl_isp, ptl_nsp
     do i=1, spall(isp)%num
        r=spall(isp)%ptl(i)%ph(1)
        z=spall(isp)%ptl(i)%ph(2)
        rz_outside=0
        if(r<eq_min_r) then
           r=eq_min_r
           rz_outside=1
        else if (r>eq_max_r)then
           r=eq_max_r
           rz_outside=1
        endif
        if(z<eq_min_z) then
           z=eq_min_z
           rz_outside=1
        else if (z>eq_max_z)then
           z=eq_max_z
           rz_outside=1
        endif
     enddo
     if(rz_outside==1)  then
        print *, 'Outside',r,z
        call err_count
     endif
  enddo
end subroutine range_check

  
subroutine err_count
  implicit none
  integer :: count=0
  integer :: sml_max_error_msg_allow
  save count

  sml_max_error_msg_allow=3000
! count may be ok with memory conflict. no need to be exact
!!!$OMP CRITICAL (ERR_COUNTX)  
  count=count+1
  if(count>sml_max_error_msg_allow) then
     print *, '###############################################################'
     print *, '# error count is exceeded maximum number allowed',sml_max_error_msg_allow,count
     print *, '###############################################################'
     stop
  endif
!!!$OMP END CRITICAL (ERR_COUNTX)
end subroutine

subroutine write_runinfo
  use sml_module
  use diag_module
  implicit none

  if(.not. sml_restart) then
     open(unit=1,file='runinfo.dat',status='replace')
  else
     open(unit=1,file='runinfo.dat',position='append')
  endif

  write(1,1000) sml_run_count, sml_gstep+1, sml_gstep/diag_1d_period+1
  close(1)
  
1000 format(8I10)
end subroutine write_runinfo

subroutine check_spnum(str,sp)
  use sml_module
  use ptl_module
  implicit none
  character (len=*) :: str
  type(species_type) :: sp

  if(sp%num > sp%maxnum) then
    print *, str, ':num, maxnum, mype', sp%num, sp%maxnum, sml_mype
  endif

end subroutine check_spnum




