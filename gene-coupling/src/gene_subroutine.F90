!#include "epik_user.inc"
#include "redef.h"
module gene_subroutine
  use par_mod
  use communications
  use perf_opt
  use eigenvalue_comp
  use initial_value_comp
  use compute_dt
  use parameters_IO
  use fullmatrix_aux
  use neo_equilibrium
  !the following modules are mainly required for external
  !data access:
  use diagnostics, only: cat_output
  use time_averages, only: avgfluxes, avgflux_stime
  use diagnostics_df, only: avgprof_stime, avgprof, istep_prof
  use diagnostics_neoclass, only: neoflux
  use checkpoint, only: checkpoint_in, read_checkpoint,&
       &reset_chpt_time,chpt_read_h5
  use geometry, only: magn_geometry, dvdx, sqrtgxx_fs,&
       &minor_r, rhostar, q0, shat, trpeps, Lref, Bref,&
       &read_tracer_namelist
  use profiles_mod, only: mref, nref, Tref, unset_in_prof_defaults
  use collisions, only: coll
  use diagnostics_fmsvd
  use codemods_mod


  ! This is the GPTL timing if built, otherwise it does nothing
  use perf_monitor

#ifdef ADIOS
  use adios_write_mod, only: adios_finalize
!  use discretization, only: mype
  use adios_io, only: adios_initialize
#endif


#ifdef COUPLE_XGC
   use coupling_core_gene, only: cce_initialize, cce_destroy
   use adios_comm_module, only: staging_read_method, staging_read_params
   use adios_read_mod
#else

#ifdef INIT_XGC
    use coupling_core_gene, only: cce_initialize, cce_destroy
    use adios_comm_module, only: staging_read_method, staging_read_params
    use adios_read_mod
#endif

#endif

  Implicit None

! --- allow 'external' manipulation of variables ---
!!  private
!!  public:: rungene, simtime_start, time

contains

  Subroutine rungene(MPI_COMM_SIM, in_dir, in_file_ex, chpt_in)

    integer :: MPI_COMM_SIM
    character(Len=128):: in_dir
    character(len=20):: in_file_ex
    character(len=*):: chpt_in
    double precision :: realtime
    logical :: sav_print_ini_msg
    !EXTERNAL :: mt_trace_start

#ifdef ADIOS
    integer :: rank, mpi_err, adios_err
#endif

    PERFON('gsub')

    par_in_dir=in_dir
    file_extension=in_file_ex

    !POMP$ INST INIT

    sav_print_ini_msg= print_ini_msg !save value which might be set from calling program

    if(in_file_ex.eq.'.dat') then
       print_ini_msg = .true.
    else
       !switch off messages for more than one parameter file
       print_ini_msg = .false.
    end if

    Call initialize_comm_sim(MPI_COMM_SIM)        ! Init constants for processor array

    IF (my_mpi_comm_world_member) THEN !inactive processors skip everything

       call init_perf_monitor(my_mpi_comm_world)
       call start_timer("FULL_SIMLATION")
       call start_timer("GENERAL_INIT")

       PERF_CONTEXT_START('autopar')
       PERFON('gini')
       Call get_systime(realtime_start)

       if (par_in_dir.ne.'skip_parfile') call read_parameters(par_in_dir)

       CALL write_codemods

       !checkpoint
       if(chpt_in.eq.'no') then
          read_checkpoint=.false.
          checkpoint_in=''
       elseif (chpt_in.ne.'') then
          checkpoint_in=chpt_in
       !else take extension being set, e.g., in parameters_IO
       end if
       if ((mype == 0).and.(print_ini_msg)) then
          if (.not.x_local) then
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,"(A)") "*********** nonlocal in radial direction ************"
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,*)
          elseif (.not.y_local) then
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,"(A)") "*********** nonlocal in toroidal direction ************"
             WRITE(*,"(A)") "*****************************************************"
             WRITE(*,*)
          endif
       endif

       if ((mype==0).and.(.not.cat_output)) call write_parameters

       call my_barrier

       call start_timer("PARALLEL_OPTIMIZATION")
       call optimize_parall_perf
       call stop_timer("PARALLEL_OPTIMIZATION")

       print_ini_msg = sav_print_ini_msg

       Call get_systime(realtime)
       initime = realtime-realtime_start

       PERFOFF
       PERF_CONTEXT_END


#ifdef ADIOS
       !call adios_set_application_id("GENE", adios_err)
       call adios_initialize(my_mpi_comm_world)
#endif

#ifdef COUPLE_XGC
      call adios_read_init_method(staging_read_method, my_mpi_comm_world, trim(staging_read_params), adios_err)
      call cce_initialize
#else

#ifdef INIT_XGC
      call adios_read_init_method(staging_read_method, my_mpi_comm_world, trim(staging_read_params), adios_err)
      call cce_initialize
#endif

#endif


       call stop_timer("GENERAL_INIT")

       select case(comp_type)
       case('EV') !eigenvalue computation
          call comp_eigenvalues

       case('IV') !initial value computation
          call set_dt_max
          if (precomp_nc) call initialize_g1_nc

          call start_timer("INITIAL_VALUE")
          call initial_value
          call stop_timer("INITIAL_VALUE")

          if (precomp_nc) call finalize_g1_nc

       case('NC')
          call comp_nc

       case('MO') !matrix output only
          call output_operator

       case('SV') !matrix output only
          if(mype==0) write(*,*) "Starting field/moment SVD."
          call svd_projection_fmom

       end select

       If (mype == 0) Call write_parameters

       Call get_systime(realtime)
       realtime = realtime - realtime_start
       If ((mype == 0).and.print_ini_msg) &
            &Write(*,"(A,F10.3,a)") "Time for GENE simulation: ",realtime, " sec"

       if (par_in_dir.ne.'skip_parfile') then
          deallocate(spec)
          call unset_in_prof_defaults()
       endif


#ifdef COUPLE_XGC
       call cce_destroy()
       call adios_read_finalize_method(staging_read_method, adios_err)
#else

#ifdef INIT_XGC
       call cce_destroy()
       call adios_read_finalize_method(staging_read_method, adios_err)
#endif

#endif

#ifdef ADIOS
       call adios_finalize(mype, adios_err)
#endif


       call stop_timer("FULL_SIMLATION")
       call finish_perf_monitor()
    ENDIF !my_mpi_comm_world_member?

    Call finalize_comm_sim

    PERFOFF  !gene_subroutine

  End Subroutine rungene

end module gene_subroutine
