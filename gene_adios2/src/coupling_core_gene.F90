module coupling_core_gene

  use coupling_base_setup
#ifdef ADIOS2
  use adios2
  use adios2_comm_module, only: engines, adios2obj
#else
  use adios_comm_module, only: staging_read_method
  use adios_read_mod
  use adios_write_mod
#endif
  use mpi

!use adios_util, only: adios_xopen, adios_xclose
  implicit none

  public:: check_coupler
  public:: initialize_coupling_engines, dump_dt
  public:: send_density, receive_field, receive_density, receive_density_3d_nogc
  public:: cce_initialize
  public:: cce_process_field, cce_process_density, write_check_file,receive_Gene_density, cce_varpi_grid

  private

  character(5) :: fld_name_XGC
# if defined(ADIOS2)
  type(adios2_variable),save :: dens_id, field_id
  type(adios2_io),save :: coupling_io, dens_io
  logical,save :: read_engine_av=.false., write_engine_av=.false.
# endif
    real, dimension(:,:),allocatable, save :: dens_out
contains

  subroutine cce_initialize()
     use discretization, only: mype

     call set_default_coupler
     open(unit=20, file='coupling.in', status='old', action='read')
     READ(20, NML=coupling)
     close(unit=20)

     if (cce_side==0) then
        cce_my_side = 'core'
        cce_other_side = 'edge'
     else
        cce_my_side = 'edge'
        cce_other_side = 'core'
     endif

     cce_surface_number = cce_last_surface - cce_first_surface + 1
     cce_node_number = cce_last_node - cce_first_node+1

     if (cce_dpot_index0) then
        fld_name_XGC='dpot0'
        if (mype.eq.0) write(*,'(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        if (mype.eq.0) write(*,'(A)')"This coupling is supposed to have a phase shift in the toridal direction"
        if (mype.eq.0) write(*,'(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     else
        fld_name_XGC='dpot1'
     endif
    
  end subroutine cce_initialize


  subroutine initialize_coupling_engines(node_number, nphi_total,offset_x, block_count)
    use communications, only: my_mpi_comm_world,mpi_comm_x
    use discretization, only: my_pespec, my_pew, my_pez, my_pev
    integer, intent(in) ::node_number, nphi_total,offset_x, block_count

  end subroutine initialize_coupling_engines


  subroutine check_coupler(gene_first, gene_nfs, gene_nnodes)
    integer :: gene_first, gene_nfs, gene_nnodes

    if (gene_first+1.ne.cce_first_surface) then
       print *, "Mismatch on first surface"
    end if

    if (cce_last_surface.gt.gene_first+gene_nfs) then
       print *, "Mismatch on last surface"
    end if

    if (gene_nnodes.ne.cce_node_number) then
       print *, "Mistmatch in number of nodes GENE ", gene_nnodes, " XGC ", cce_node_number
    end if

  end subroutine check_coupler


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Density portion of the coupling
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Send the density to the other side
  subroutine send_density(density, iphi, nphi_total, block_count, block_start, block_end, comm)
    use par_mod, only: time
    use par_other, only: itime
    use discretization, only: mype

    integer, intent(in) :: comm
    integer, intent(in) :: iphi
    integer, intent(in) :: nphi_total
    integer, intent(in) :: block_count, block_start,block_end
    real, dimension(0:nphi_total-1, block_start:block_end), intent(in) :: density

    ! GENE gives the input density to this function with flipped dimension orders
    real, dimension(:, :), allocatable :: tmp
#ifdef ADIOS2
    integer(8), dimension(2) :: gdims, goffset, ldims
!    type(adios2_io) :: io
!    type(adios2_engine) :: engine
!    type(adios2_variable) :: varid
#else
    integer(8) :: adios_handle, adios_groupsize, adios_totalsize
#endif
    integer :: adios_err, i, maxplane
    character(16)::fld_name='gene_density'

    maxplane = nphi_total - 1


    if (cce_comm_density_mode.eq.1.or.cce_comm_density_mode.eq.2) then
#ifdef ADIOS2       
       if (.not.write_engine_av) then
          !engines(1) is the writer
          gdims(1) = cce_node_number
          gdims(2) = nphi_total
          goffset(1) = block_start
          goffset(2) = 0
          ldims(1) = block_count
          ldims(2) = nphi_total

        print *, mype,  ' runs the coupling '


          call adios2_declare_io(dens_io,adios2obj,'density_coupling',adios_err)
          call adios2_define_variable(dens_id, dens_io,fld_name,&
               & adios2_type_dp, 2, gdims, goffset,&
               & ldims,adios2_constant_dims, adios_err)
          call adios2_open(engines(1), dens_io, trim(cce_folder)//'/'&
               & //'gene_density.bp',adios2_mode_write,&
               & comm, adios_err)
          write_engine_av=.true.
          allocate(dens_out(block_start:block_end, 0:maxplane))
       endif
#endif
       if (.not.allocated(dens_out)) allocate(dens_out(block_start:block_end, 0:maxplane))
       do i=0, maxplane
          dens_out(block_start:block_end, i) = density(i, block_start:block_end)
       end do

#ifdef ADIOS2
       call adios2_begin_step(engines(1),adios2_step_mode_append,adios_err)

       do i = 1,10 
       print *,  "First 10 -> dim1 X dim2: ", block_start + i - 1, 0, " is: ", dens_out(block_start + i - 1, 0)
       end do

       do i = 1,10 
       print *,  "Last 10 -> dim1 X dim2: ", block_end - 10 + i, maxplane, " is :", dens_out(block_end - 10 + i, maxplane)
       end do

       call adios2_put(engines(1),dens_id,dens_out,adios_err)
       call adios2_end_step(engines(1), adios_err)
#else 
       call get_cce_filename("density", staging_read_method, cce_my_side)
       call adios_open(adios_handle,'coupling', cce_filename, 'w', comm, adios_err)
#include "gwrite_coupling.fh"
       call adios_close(adios_handle, adios_err)
       deallocate(dens_out)
#endif
    endif

  end subroutine send_density


  subroutine write_check_file(comm)
    use discretization, only: mype
    integer, intent(in) :: comm
    integer :: ierr
#ifdef ADIOS
    call mpi_barrier(comm, ierr)
    if (staging_read_method .eq. ADIOS_READ_METHOD_BP) then
       if (mype.eq.0) then
          call get_cce_filename("density", staging_read_method, cce_my_side)
          open(20, file=cce_lock_filename, status="new", action="write")
          close(20)
       end if
    end if
#endif
  end subroutine write_check_file


  subroutine receive_GENE_density(data_block, block_start, block_end, block_count, iphi, comm)
     integer, intent(in) :: iphi, block_start, block_end, block_count
     integer, intent(in) :: comm
     real, dimension(block_start:block_end, 1) :: data_block
     real, dimension(:,:), allocatable :: tmp
#if 0
     integer(8) :: adios_handle
     integer(8) :: bb_sel
     integer(8), dimension(2) :: bounds, counts
     integer :: adios_err


     ! change second index afterswitching node order in XGC
     bounds(1) = int(block_start, kind=8)
     counts(1) = int(block_count, kind=8)
     bounds(2) = iphi
     counts(2) = 1
     allocate(tmp(block_start:block_end, 1))

     call cce_couple_open(adios_handle, "density", staging_read_method, cce_my_side, comm, adios_err)
     call adios_selection_boundingbox(bb_sel, 2, bounds, counts)
     call adios_schedule_read(adios_handle, bb_sel, "data", 0, 1, tmp, adios_err)
     call adios_perform_reads(adios_handle, adios_err)
     call adios_read_close(adios_handle, adios_err)
     call adios_selection_delete(bb_sel)

     data_block(block_start:block_end, 1) = tmp(block_start:block_end, 1)
     deallocate(tmp)
#endif
  end subroutine receive_GENE_density

  subroutine cce_process_density()
     cce_step = cce_step + 1
  end subroutine cce_process_density


  subroutine receive_field(data_block, block_start, block_end, block_count, nphi_total, comm)
    use discretization, only:mype
     integer, intent(in) :: block_start
     integer, intent(in) :: block_end
     integer, intent(in) :: block_count
     integer, intent(in) :: nphi_total
     integer, intent(in) :: comm
     real, dimension(0:nphi_total-1, block_start:block_end), intent(out) :: data_block


     ! Again tmp is ordered opposite of what we'll hand back to GENE
     real, dimension(:, :), allocatable :: tmp
#ifdef ADIOS2
!     type(adios2_io) :: io
!     type(adios2_engine) :: engine
     type(adios2_variable) :: varid
#else              
     integer(8) :: adios_handle, bb_sel
#endif
     integer(8), dimension(2) :: bounds, counts
     integer :: maxplane, adios_err, i, i_p
     maxplane = nphi_total - 1


     if (cce_side.eq.0.and.cce_comm_field_mode.GT.1) then
        ! change second index afterswitching node order in XGC
        bounds(1) = int(cce_first_node-1+block_start, kind=8)
        counts(1) = int(block_count, kind=8)
        bounds(2) = 0
        counts(2) = nphi_total
        allocate(tmp(block_start:block_end, 0:maxplane))
#ifdef ADIOS2
        if (.not.read_engine_av)then
           call adios2_declare_io(coupling_io, adios2obj, &
                & "field_from_coupling", adios_err)
           call adios2_open(engines(2), coupling_io, & 
                & trim(cce_folder)//'/'//'cpl_field.bp', &
                & adios2_mode_read, comm, adios_err)
           read_engine_av=.true.
           IF (MYPE.EQ.0)print *, 'created engine for field'
        endif

        call adios2_begin_step(engines(2),adios2_step_mode_read,adios_err)
        call adios2_inquire_variable(field_id,coupling_io,'cpl_field',adios_err)
        IF (MYPE.EQ.0) THEN
           if (.not.field_id%valid)then
              print *, fld_name_XGC, ' in field not found'
           else
              print *, fld_name_XGC, ' in field found'
           endif
        endif

        call adios2_set_selection(field_id, 2, bounds, counts, adios_err)
        call adios2_get(engines(2), field_id, tmp, adios2_mode_deferred, adios_err)
        call adios2_end_step(engines(2), adios_err)

        IF (MYPE.EQ.0) THEN
           print *, minval(tmp),maxval(tmp)
        Endif

        do i_p=0,maxplane
           if(any(isnan(tmp(:,i_p))))print *,"nan on plane",i_p,'rank',mype
           if(any(tmp(:,i_p)>huge(tmp(1,1))))print *,"inf on plane",i_p,'rank',mype
        end do
#else
        call cce_couple_open(adios_handle, fld_name_XGC, staging_read_method, cce_other_side, comm, adios_err, lockname="field")
        call adios_selection_boundingbox(bb_sel, 2, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "data", 0, 1, tmp, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_read_close(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)
#endif
        do i=0, maxplane
           data_block(i, block_start:block_end) = tmp(block_start:block_end, i)
        end do

        deallocate(tmp)
     endif

  end subroutine receive_field


  subroutine receive_density(data_block, block_start, block_end, block_count, nphi_total, comm)
     use discretization, only: mype

     integer, intent(in) :: block_start
     integer, intent(in) :: block_end
     integer, intent(in) :: block_count
     integer, intent(in) :: nphi_total
     integer, intent(in) :: comm
     real, dimension(0:nphi_total-1, block_start:block_end), intent(out) :: data_block


     ! Again tmp is ordered opposite of what we'll hand back to GENE
     real, dimension(:, :), allocatable :: tmp


     integer(8) :: adios_handle, bb_sel
     integer(8), dimension(2) :: bounds, counts
     integer :: maxplane, adios_err, i
     maxplane = nphi_total - 1


     if (cce_side.eq.0.and.cce_comm_field_mode.GT.1) then
        ! change second index afterswitching node order in XGC
        bounds(1) = int(block_start, kind=8)
        counts(1) = int(block_count, kind=8)
        bounds(2) = 0
        counts(2) = nphi_total

        allocate(tmp(block_start:block_end, 0:maxplane))
#ifdef ADIOS2
#else
        call cce_couple_open(adios_handle, "init", staging_read_method, &
                           & cce_other_side, comm, adios_err, lockname="init")
        call adios_selection_boundingbox(bb_sel, 2, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "data", 0, 1, tmp, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_read_close(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)
#endif
        do i=0, maxplane
           data_block(i, block_start:block_end) = tmp(block_start:block_end, i)
        end do

        deallocate(tmp)
     endif

  end subroutine receive_density

  subroutine receive_density_3D_noGC(data_block, block_start, block_end, block_count, nphi_total, comm)
     use discretization, only: mype

     integer, intent(in) :: block_start
     integer, intent(in) :: block_end
     integer, intent(in) :: block_count
     integer, intent(in) :: nphi_total
     integer, intent(in) :: comm
     real, dimension(0:nphi_total-1, block_start:block_end), intent(out) :: data_block


     ! Again tmp is ordered opposite of what we'll hand back to GENE
     real, dimension(:, :), allocatable :: tmp


     integer(8) :: adios_handle, bb_sel
     integer(8), dimension(2) :: bounds, counts
     integer :: maxplane, adios_err, i
     maxplane = nphi_total - 1


     if (cce_side.eq.0.and.cce_comm_field_mode.GT.1) then
        ! change second index afterswitching node order in XGC
        bounds(2) = int(cce_first_node-1+block_start, kind=8)
        counts(2) = int(block_count, kind=8)
        bounds(1) = 0
        counts(1) = nphi_total

        !allocate(tmp(block_start:block_end, 0:maxplane))
#ifdef ADIOS2
        data_block=0
#else
        call cce_couple_open(adios_handle, "XGC3D", staging_read_method, &
                           & cce_other_side, comm, adios_err, lockname="XGCinit")
        call adios_selection_boundingbox(bb_sel, 2, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "iden_nogyro", 0, 1, data_block, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_read_close(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)
#endif
        !do i=0, maxplane
        !   data_block(i, block_start:block_end) = tmp(block_start:block_end, i)
        !end do

        !deallocate(tmp)
     endif

  end subroutine receive_density_3D_noGC


  subroutine cce_process_field()
    cce_field_step = cce_field_step + 1
  end subroutine cce_process_field


  function cce_varpi_grid(ipsi) result(varpi)
    real :: varpi
    integer :: ipsi

    !call cce_initialize()
    if (cce_density_model.eq.6)then
       if (cce_npsi>0)then
           print *, 'Not implemented'
           stop
       else
       ! Linear weight
          if (cce_side.EQ.0) then
             if (ipsi.le.cce_first_surface_coupling_axis) then
                varpi=0D0
             else if (ipsi.le.cce_last_surface_coupling_axis) then
                varpi=dble(ipsi-cce_first_surface_coupling_axis)/ &
                    & dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
             else if (ipsi.le.cce_first_surface_coupling) then
                varpi=1D0
             else if (ipsi.gt.cce_last_surface_coupling) then
                varpi=0D0
             else
                varpi=1D0-dble(ipsi-cce_first_surface_coupling)/ &
                    & dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
             endif

          else
             print *, 'GENE is for the core, put cce_side=0'
             stop
          endif

       endif

    endif

  end function cce_varpi_grid


  subroutine dump_dt(dt)
     use discretization, only: mype
     use par_poloidal_planes, only: cref
     use geometry, only: Lref
     use par_mod, only: time, itime
     real, intent(in) :: dt

     integer(8) :: adios_handle, adios_groupsize, adios_totalsize
     integer :: adios_err
#if 0
     if (mype.eq.0) then
        call get_cce_filename("dt", staging_read_method, cce_my_side)
        call adios_open(adios_handle,'dt', cce_filename, 'w', MPI_COMM_SELF, adios_err)
#include "gwrite_dt.fh"
        call adios_close(adios_handle, adios_err)

        ! write check file
        if (staging_read_method .eq. ADIOS_READ_METHOD_BP) then
           call get_cce_filename("dt", staging_read_method, cce_my_side)
           open(20, file=cce_lock_filename, status="new", action="write")
           close(20)
	end if

     end if
#endif
  end subroutine dump_dt


end module coupling_core_gene
