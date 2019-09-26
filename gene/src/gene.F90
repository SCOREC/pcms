#include "redef.h"
#include "intrinsic_sizes.h"
!>Program to start several gene runs simultaneously. 
!!\todo embellish/sort the standard output (all parallel genes write into the same file at the moment!)
Program gene
  use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
  use mpi
  use gene_scan
  use communications
  use par_mod
  use parameters_IO
  use check_parameters
  use gene_subroutine
  use file_io
  use pinning_mod
  use diagnostics, only: NRGFILE, cat_output
  use, intrinsic :: iso_c_binding
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
#ifdef MEMORY_CHECKER
  use memory_checker
#endif
  implicit none

#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ attributes offload:mic :: perfon,perfoff,perfinit,perfout
#endif

  integer :: ierr
#ifdef MEMORY_CHECKER
  integer(kind=c_long) :: max_allocated,still_allocated
  integer(kind=c_int) :: unfreed_blocks
#endif
  double precision:: wtime_start
  logical:: mult_par=.true.
  integer:: i_proc,i_group
  logical:: performance_output=.false.
  integer:: gene_comm, comm_parall
  character(len=FILEEXT_MAX):: f_ext
  character(len=FILENAME_MAX):: ch_in

  !set git_branch/git_master here as par_other.F90 might not be
  !recompiled after code modifications (gene_subroutine, however,
  !almost always is)
#if defined(GIT_BRANCH)
  git_branch = GIT_BRANCH
  git_master = GIT_MASTER
#endif

#ifdef WITHOMP
  call mpi_init_thread(MPI_THREAD_MULTIPLE, omp_level,ierr)
#ifdef __INTEL_COMPILER
  call pinning()
#endif
#else
  call mpi_init(ierr)
  omp_level = MPI_THREAD_SINGLE
#endif
  if (ierr /= 0) stop 'mpi_init failed!'

  call check_for_scan(mult_par)

  call initialize_comm_scan(gene_comm,comm_parall)
  LIKWID_INIT
#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
    !DEC$ offload begin target(mic:WHICHMIC)
    PERFINIT
    PERFON('mmain')
    !DEC$ end offload
#endif

  PERFINIT         ! Initialize performance monitoring
  PERFON('GENE')

  if (mype_gl.eq.0) then
     write(*,*)
     write(*,"(A)") "*****************************************************"
     write(*,"(3A)") "***** This is GENE 11 (release ",&
          trim(release),") *******"
     if (git_branch.ne.'') then !show abbreviated hash
        write(*,"(A,A7,A)") "***** GIT branch hash: ", git_branch, '                *******'
        write(*,"(A,A7,A)") "***** GIT master hash: ", git_master, '                *******'
     endif
     write(*,"(A)") "*****************************************************"
     write(*,*)
  endif
  
  wtime_start=MPI_Wtime() 
  
  call check_for_diagdir

  if(.not.mult_par) then
     !only one gene simulation
     f_ext='.dat'
     ch_in=''
     call rungene(gene_comm, par_in_dir, f_ext, ch_in)
  else
     !list of gene simulations
     call run_scan(gene_comm, comm_parall)
  end if
  call erase_stop_file
  call create_finished_file

  if (mype_gl.eq.0) then
     write(*,"(A,F10.3,A)") "Total wallclock time for GENE: ", MPI_Wtime()-wtime_start, " sec"
     if(mult_par) write(*,"(A,F10.3,A)") "Percentage of idle time at the end of the scan due to load imbalancing: ", &
          idle_global/n_parallel_sims/n_procs_sim/(MPI_Wtime()-wtime_start)*100, " %"
  end if

  PERFOFF !GENE

#if defined(WITH_MIC_NONLIN) && defined(WITH_OFFLOAD) && defined(WITHPERF)
  !DEC$ offload begin target(mic:WHICHMIC)
  PERFOFF
  PERFOUT('mmain')
  !DEC$ end offload
#endif

  ! output performance data
  if (performance_output.and.(n_procs_sim.le.64)) then
     do i_group=0,n_parallel_sims-1
        do i_proc=0,n_procs_sim-1
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           if ((i_proc.eq.mype).and.(i_group.eq.my_sim)) then
              write(*,"(A,I4)") "Process ",mype
              PERFOUT('t_loop')
              !close(TEMPFILE)
             call flush(OUTPUT_UNIT)
           end if
        end do
     end do
  else
     if (mype_gl.eq.0) then
        PERFOUT('GENE') 
     endif
  end if
  LIKWID_CLOSE

  call finalize_comm_scan(gene_comm,comm_parall)
  call mpi_finalize(ierr)

#ifdef MEMORY_CHECKER
  call get_malloc_stat(max_allocated,unfreed_blocks,still_allocated)
  write(*,"(I7,A,I10,A)") mype,": max_allocated = ",max_allocated," bytes"
#endif
End Program gene
