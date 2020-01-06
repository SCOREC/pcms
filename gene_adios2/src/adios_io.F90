! This module replaces GENE's use of hlst_adios_checkpoint.F90.
! The original HLST (High-level support team) version is not flexible enough to do what is needed in the ECP WDM project.

#include "intrinsic_sizes.h"

module adios_io

    use mpi
    use adios_write_mod
    use adios_read_mod
    use par_in, only: file_extension
    use par_mod, only: &
        adios_xml, diagdir, &
        nx0, nky0, nz0, nv0, nw0, n_spec, n_fields, &
        li1, li2, lj1, lj2, lk1, lk2, ll1, ll2, lm1, lm2, ln1, ln2, &
        itime, time, dt, g_1
    use aux_fields, only: emfields
    use antenna, only: Apar_pre_antenna, antenna_type

    ! This is the GPTL timing if built, otherwise it does nothing
    use perf_monitor

    implicit none


    ! Common ADIOS variables, should only need these in situations like the mom group
    integer :: adios_err
    integer(8) :: adios_handle, adios_groupsize, adios_totalsize


    ! Dimensions for each MPI process's data
    integer, dimension(1:6) :: local_dims


    ! File names
    character(len=FILENAME_MAX) :: adios_field_file
    character(len=FILENAME_MAX) :: adios_field_file_ky0
    character(len=FILENAME_MAX) :: adios_mom_file
    character(len=FILENAME_MAX) :: adios_restart_file


contains

  subroutine adios_initialize(comm)
    integer, intent(in) :: comm
    integer :: adios_err

    local_dims(1) = li2 - li1 + 1
    local_dims(2) = lj2 - lj1 + 1
    local_dims(3) = lk2 - lk1 + 1
    local_dims(4) = ll2 - ll1 + 1
    local_dims(5) = lm2 - lm1 + 1
    local_dims(6) = ln2 - ln1 + 1

    ! GENE appends "/" to diagdir automatically elsewhere (parameters_IO.F90)
    write(adios_field_file, '(A,A)') trim(diagdir), "field"//trim(file_extension)//".bp"
    write(adios_field_file_ky0, '(A,A)') trim(diagdir), "field-ky0"//trim(file_extension)//".bp"
    write(adios_mom_file, '(A,A)')   trim(diagdir), "mom_"//trim(file_extension)//".bp"

    call adios_init(adios_xml, comm, adios_err)

  end subroutine adios_initialize


  subroutine adios_checkpoint_read_init(comm, param)
    integer, intent(in) :: comm
    character(len=*), intent(in), optional :: param
    integer :: adios_err

    if (present(param)) then
       call adios_read_init_method(ADIOS_READ_METHOD_BP, comm, param, adios_err)
    else
       call adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "verbose=2", adios_err)
    end if
  end subroutine adios_checkpoint_read_init


  subroutine adios_checkpoint_write(fname, comm)
    integer, intent(in) :: comm
    character(len=*), intent(in) :: fname
    character(len=FILENAME_MAX) :: outfile
    integer :: adios_err
    write(outfile, "(a, a, i0, a)") fname, "_", itime, ".bp"

    call start_timer("CHECKPOINT_WRITE_ADIOS")
    call adios_open(adios_handle, "restart", TRIM(outfile), "w", comm, adios_err)
#include "gwrite_restart.fh"
    call adios_close(adios_handle, adios_err)
    call stop_timer("CHECKPOINT_WRITE_ADIOS")

  end subroutine adios_checkpoint_write


  subroutine adios_g_read(g, bounds, comm, adios_handle, infile)
    complex, dimension(:, :, :, :, :, :), intent(out) :: g
    integer(8), dimension(6), intent(in) :: bounds
    integer(8), intent(in) :: adios_handle
    integer, intent(in) :: comm
    character(len=*), intent(in), optional :: infile
    integer(8) :: bb_sel
    integer :: step=0, nstep=1, adios_err

    call adios_selection_boundingbox(bb_sel, 6, bounds, shape(g, kind=8))
    call start_timer("READ RESTART (g)")
    call adios_schedule_read(adios_handle, bb_sel, "g", step, nstep, g, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call adios_read_close(adios_handle, adios_err)
    call stop_timer("READ RESTART (g)")
    call adios_selection_delete(bb_sel)
  end subroutine adios_g_read


  subroutine adios_checkpoint_dims(infile, comm, adios_handle, dims, time, dt)
    character(len=*), intent(in) :: infile
    integer, intent(in) :: comm
    integer(8), intent(inout) :: adios_handle
    integer, dimension(6), intent(out) :: dims
    real(8), intent(out) :: time, dt
    integer(8) :: wb_sel
    integer :: method=ADIOS_READ_METHOD_BP, step=0, nstep=1, adios_err

    call adios_selection_writeblock(wb_sel, 0)
    call start_timer("READ RESTART (dims/time)")
    call adios_read_open_file(adios_handle, infile, method, comm, adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_x", step, nstep, dims(1), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_ky", step, nstep, dims(2), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_z", step, nstep, dims(3), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_vz", step, nstep, dims(4), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_vx", step, nstep, dims(5), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "global_num_species", step, nstep, dims(6), adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "time", step, nstep, time, adios_err)
    call adios_schedule_read(adios_handle, wb_sel, "dt", step, nstep, dt, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call stop_timer("READ RESTART (dims/time)")
    call adios_selection_delete(wb_sel)
  end subroutine adios_checkpoint_dims


  subroutine ad_aopen(adios_handle, groupname, filename, comm, step)
    integer(8), intent(inout) :: adios_handle
    character(len=*), intent(in) :: groupname
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(in) :: step
    integer :: adios_err

    if (step.eq.0) then
       call adios_open(adios_handle, groupname, filename, "w", comm, adios_err)
    else
       call adios_open(adios_handle, groupname, filename, "a", comm, adios_err)
    end if
  end subroutine ad_aopen


  ! Write GENE's "field" file. (4.2 in user manual)
  subroutine adios_field_write(field_step, comm)
    integer, intent(in) :: comm, field_step
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2) :: tmp_field
    integer(8) :: adios_handle
    call start_timer("DIAG_FIELD_ADIOS")
    call ad_aopen(adios_handle, "field", adios_field_file, comm, field_step)
#include "gwrite_field.fh"
    call adios_close(adios_handle, adios_err)
    call stop_timer("DIAG_FIELD_ADIOS")
  end subroutine adios_field_write


  subroutine adios_field_write_ky0(field_step, comm)
    integer, intent(in) :: comm, field_step
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2) :: tmp_field
    integer(8) :: adios_handle
    call start_timer("DIAG_FIELD_ADIOS_KY0")
    call ad_aopen(adios_handle, "field_ky0", adios_field_file_ky0, comm, field_step)
#include "gwrite_field_ky0.fh"
    call adios_close(adios_handle, adios_err)
    call stop_timer("DIAG_FIELD_ADIOS_KY0")
  end subroutine adios_field_write_ky0


  subroutine adios_mom_open(mom_step, comm, adios_handle)
    integer, intent(in) :: comm, mom_step
    integer(8) :: adios_handle
    call start_timer("DIAG_MOM_ADIOS")
    call ad_aopen(adios_handle, "mom", adios_mom_file, comm, mom_step)
    call stop_timer("DIAG_MOM_ADIOS")
  end subroutine adios_mom_open


  subroutine adios_mom_close(comm, adios_handle)
    integer, intent(in) :: comm
    integer(8) :: adios_handle
    call start_timer("DIAG_MOM_ADIOS")
    call adios_close(adios_handle, adios_err)
    call stop_timer("DIAG_MOM_ADIOS")
  end subroutine adios_mom_close

end module adios_io
