!----------------------PERFORMANCE STUFF------------------------------------------
module perf_monitor
use mpi
#ifdef WITH_GPTL
use gptl, only: gptlprocess_namelist, gptlsetoption, gptlverbose, gptlinitialize, &
                gptlstart, gptlstop, gptlfinalize, gptlpr_summary_file
#endif
use par_in, only: diagdir
use file_io, only: get_unit_nr


implicit none

integer :: timer_comm
character(len=100) :: timer_outdir=""

contains

  subroutine timer_setdir(dir)
    character(len=*), intent(in) :: dir
    integer :: rank, ierr, fileunit
#ifdef WITH_GPTL
    call mpi_comm_rank(timer_comm, rank, ierr)
    if (rank.eq. 0) then
       call get_unit_nr(fileunit)
       open(fileunit, file=TRIM(dir)//"/test", iostat=ierr)
       if (ierr.ne.0) then
          call execute_command_line ('mkdir -p ' // adjustl(trim(dir)))
       else
          close(fileunit,status="DELETE")
       end if

       write (timer_outdir, "(A)") trim(dir)
       print *,"Saving performance in ", dir

    endif
#endif
  end subroutine timer_setdir


  subroutine init_perf_monitor(comm)
    integer, intent(in) :: comm
#ifdef WITH_GPTL
    integer :: i, rank, ierr
    logical :: masterproc, exists
    timer_comm = comm
    call timer_setdir(trim(diagdir)//"/timing/")

    inquire(FILE="perf_gptl", EXIST=exists)
    if (exists) then
       call gptlprocess_namelist("perf_gptl", 100, ierr)
    endif
    ierr = gptlsetoption(gptlverbose, 0)
    ierr = gptlinitialize()
#endif
  end subroutine init_perf_monitor


  subroutine start_timer(label)
    character(len=*), intent(in) :: label
#ifdef WITH_GPTL
    integer :: ierr
    ierr = gptlstart(label)
#endif
  end subroutine start_timer


  subroutine stop_timer(label)
    character(len=*), intent(in) :: label
#ifdef WITH_GPTL
    integer :: ierr
    ierr = gptlstop(label)
#endif
  end subroutine stop_timer


  subroutine flush_perf_monitor(istep)
    integer, intent(in) :: istep
#ifdef WITH_GPTL
    integer :: ierr, rank
    logical :: exists
    character (len=100) :: gfilename

    ierr = gptlstart("sync1_t_prf")
    call mpi_barrier(timer_comm,ierr)
    ierr = gptlstop("sync1_t_prf")

    if (istep < 0) then
       write(gfilename, "(A, A)") trim(timer_outdir), "init.summary.txt"
    else
       write(gfilename, "(A, i9.9, A)") trim(timer_outdir), istep, ".summary.txt"
    endif

    ierr = gptlpr_summary_file(timer_comm, trim(gfilename))
    if (ierr /= 0) then
       print *, ierr, 'summary.F90: error from gptlpr_summary_file'
       stop 1
    end if

    ierr = gptlstart("sync2_t_prf")
    call mpi_barrier(timer_comm, ierr)
    ierr = gptlstop("sync2_t_prf")
#endif

  end subroutine flush_perf_monitor


  subroutine finish_perf_monitor()
#ifdef WITH_GPTL
    integer :: ierr
    logical :: exists
    character (len=100) :: gfilename
    write(gfilename, "(A, A)") trim(timer_outdir), "all.txt"
    ierr = gptlpr_summary_file(timer_comm, trim(gfilename))
    ierr = gptlfinalize()
#endif
  end subroutine finish_perf_monitor

end module perf_monitor
