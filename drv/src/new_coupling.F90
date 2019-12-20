module new_coupling
#ifdef ADIOS2
  use adios2
  use adios2_comm_module
#endif
  use adios_comm_module
  use adios_read_mod
  use adios_write_mod
  use hdf5
  use mpi


!  use sim_param                 ! substitute routine to sml_module
  IMPLICIT NONE
  logical :: cce_nzeromode

  ! program metadata
  character(:), allocatable :: cce_input_file
  character(:), allocatable :: cce_folder
  character(:), allocatable :: cce_output_dir
  logical :: cce_initialized

  ! simulation metadata
  integer :: cce_step
  integer :: cce_side ! 0:core (GENE/XGC), 1: edge(XGC)
  integer :: cce_density_model

  character(:), allocatable :: cce_adios_group
  character(8) :: cce_my_side
  character(8) :: cce_other_side
  character(8) :: cce_density_step_string
  character(8) :: cce_field_step_string
  character(8) :: cce_plane_string

  integer, dimension(:), allocatable :: cce_surface_first_node
  integer, dimension(:), allocatable :: cce_surface_last_node

  integer :: cce_first_surface
  integer :: cce_last_surface
  integer :: cce_surface_count
  integer :: cce_all_surface_number

  ! density field global vals
  integer :: cce_density_first_node
  integer :: cce_density_last_node
  integer :: cce_density_node_count
  integer :: cce_node_number                    !added for variable signature consistency
  integer :: cce_first_node                     !added for variable signature consistency
  integer :: cce_last_node                      !added for variable signature consistency
!  integer :: staging_read_method = ADIOS_READ_METHOD_BP

  integer :: cce_density_step
  integer :: cce_field_step
  integer :: cce_field_model

  integer :: cce_comm_density_mode
  integer :: cce_comm_field_mode

  ! potential field global vals
  integer :: cce_field_first_node
  integer :: cce_field_last_node
  integer :: cce_field_node_count
  integer :: cce_field_node_number      !variable signature  consistency

  integer :: cce_field_first_surface
  integer :: cce_field_last_surface
  integer :: cce_field_surface_count

  integer :: cce_first_surface_coupling
  integer :: cce_last_surface_coupling

  integer :: cce_first_surface_coupling_axis
  integer :: cce_last_surface_coupling_axis

  ! physical data (fields)
  real, dimension(:), allocatable :: cce_density      ! density is represented as a 1D array
  real, dimension(:), allocatable :: cce_pot0
  real, dimension(:), allocatable :: cce_dpot0
  real, dimension(:), allocatable :: cce_dpot1

  !
  real(8) :: cce_alpha
  real(8) :: cce_dt ! only in GENE

  !
  integer :: cce_npsi
  real(8), dimension(:), allocatable :: cce_varpi
  real(8), dimension(:), allocatable :: cce_psi

  logical :: cce_dpot_index0
  logical :: cce_bcast_dpot
  logical :: cce_bcast_pot0

  !
  character(256) :: cce_lock_filename
!  character(256) :: cce_folder
  character(256) :: cce_filename

  ! namelists
  namelist /coupling/ cce_side, &
    cce_npsi, &
    cce_alpha, &
    cce_density_model, &
    cce_folder, &
    cce_first_surface, &
    cce_last_surface, &
    cce_density_first_node, &           ! cce_first_node
    cce_density_last_node, &            ! cce_last_node
    cce_comm_density_mode, &
    cce_all_surface_number, &
    cce_field_model, &
    cce_comm_field_mode, &
    cce_first_surface_coupling, &
    cce_last_surface_coupling, &
    cce_field_first_surface, &          ! cce_first_surface_field
    cce_field_last_surface, &           ! cce_last_surface_field
    cce_first_surface_coupling_axis, &
    cce_last_surface_coupling_axis, &
    cce_nzeromode, &
    cce_dpot_index0, &
    cce_dt

  namelist /surfaces/ cce_surface_first_node, cce_surface_last_node
  namelist /varpi/ cce_varpi, cce_psi

#ifdef ADIOS2
  type(adios2_io),save :: send_io, read_io
  type(adios2_engine),save :: send_engine, read_engine
  type(adios2_variable),save:: field_id
#endif
  real*8,dimension(:),allocatable,save::field_out

  interface cce_process_field
    module procedure cce_process_field_GENE, cce_process_field_XGC
  end interface

 interface cce_process_density
    module procedure cce_process_density_GENE, cce_process_density_XGC
  end interface 

  public:: check_coupler
  public:: initialize_coupling_engines, dump_dt
  public:: send_density, receive_field, receive_density, receive_density_3d_nogc
  public:: cce_initialize
  public:: cce_process_field, cce_process_density, write_check_file,receive_Gene_density, cce_varpi_grid

  private

  character(5) :: fld_name_XGC
# if defined(ADIOS2)
  type(adios2_variable),save :: dens_id !, field_id
  type(adios2_io),save :: coupling_io, dens_io
  logical,save :: read_engine_av=.false., write_engine_av=.false.
#endif
  real, dimension(:,:),allocatable, save :: dens_out


contains


subroutine cce_set_defaults
    implicit none
       cce_alpha=0.5D0
       cce_density_model=0
       cce_step=1 !In case of restart, one might want to change this
       cce_density_step = 0
       cce_field_step=1
       cce_first_node=10
       cce_last_node=0
       cce_comm_density_mode=2
       cce_all_surface_number=1
       cce_first_surface_coupling=-1
       cce_last_surface_coupling=-1
       cce_field_model=0
       cce_comm_field_mode=0
       cce_npsi=-1
       cce_nzeromode=.false.
       cce_dpot_index0=.false.
end subroutine cce_set_defaults


! GENE's version of cce_initialize
subroutine cce_initialize()
     use sim_param
    integer :: ipsi

     call cce_set_defaults
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

     cce_surface_count = cce_last_surface - cce_first_surface + 1
     cce_node_number = cce_last_node - cce_first_node+1

     if (cce_dpot_index0) then
        fld_name_XGC='dpot0'
        if (mype.eq.0) write(*,'(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        if (mype.eq.0) write(*,'(A)')"This coupling is supposed to have a phase shift in the toridal direction"
        if (mype.eq.0) write(*,'(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     else
        fld_name_XGC='dpot1'
     endif

    if(.not.allocated(cce_density))then

!      call cce_set_defaults

       open(unit=20,file='coupling.in', status='old',action='read')
       READ(20, NML=coupling)
       close(unit=20)

       if(cce_side==0)then
          cce_my_side='core'
          cce_other_side='edge'
       else
          cce_my_side='edge'
          cce_other_side='core'
       endif

       cce_surface_count=cce_last_surface-cce_first_surface+1

       cce_node_number=cce_last_node-cce_first_node+1

       cce_field_node_number=cce_node_number
       cce_field_first_node=cce_first_node
       cce_field_last_node =cce_last_node

       if(cce_node_number.LE.0)then
          allocate(cce_surface_first_node(cce_all_surface_number))
          allocate(cce_surface_last_node(cce_all_surface_number))

          open(unit=20,file='coupling.in', status='old',action='read')
          READ(20, NML=surfaces)
          close(unit=20)

          cce_first_node=cce_surface_first_node(cce_first_surface)
          cce_last_node=cce_surface_last_node(cce_last_surface)
          cce_node_number=cce_last_node-cce_first_node+1

          cce_field_first_node=cce_surface_first_node(1)
          cce_field_last_node =cce_surface_last_node(cce_all_surface_number)
          cce_field_node_number=cce_field_last_node-cce_field_first_node+1
       endif

       allocate(cce_density(cce_first_node:cce_last_node))

       if(cce_comm_field_mode.NE.0)then
#ifdef XGC_COUPLING_CORE_EDGE_FIELD
          allocate(cce_dpot0(cce_field_first_node:cce_field_last_node))
          allocate(cce_dpot1(cce_field_first_node:cce_field_last_node))
          allocate(cce_pot0(cce_field_first_node:cce_field_last_node))
#else
          print *,'You need to compile with -DXGC_COUPLING_CORE_EDGE_FIELD'
#endif
       endif
!   print *,sml_intpl_mype,'first-last',cce_first_node,cce_last_node

       if(cce_first_surface_coupling.lt.1)cce_first_surface_coupling=cce_first_surface
       if(cce_last_surface_coupling.lt.1)cce_last_surface_coupling=cce_last_surface

 !  if(sml_restart)then
 !   cce_step=2*sml_istep
 !   cce_field_step=cce_step
 !  endif

       if(cce_npsi>0)then
          allocate(cce_varpi(cce_npsi))
          allocate(cce_psi(cce_npsi))
          open(unit=20,file='coupling.in', status='old',action='read')
          READ(20, NML=varpi)
          close(unit=20)
          cce_first_surface_coupling=-1
          cce_last_surface_coupling=-1
          do ipsi=1,cce_npsi
             if(cce_first_surface_coupling.eq.-1.and.cce_varpi(ipsi)>0.0001)then
                cce_first_surface_coupling=ipsi
             endif
             if(cce_last_surface_coupling.eq.-1.and.cce_varpi(ipsi)>0.9999)then
                cce_last_surface_coupling=ipsi
             endif
          enddo
          print *,"psi_min=",cce_first_surface_coupling,"psi_max=",cce_last_surface_coupling
          if(cce_side.EQ.0)cce_varpi(:)=1D0-cce_varpi(:)
       endif

    endif

!    call adios2_declare_io(read_io, adios2obj, "density_coupling", err)
!    call adios2_open(read_engine, read_io, trim(cce_folder) // &
!      &'/' // 'density.bp', adios2_mode_read, sml_intpl_comm, err)

end subroutine cce_initialize


! Failed merged cce_initialize attempt
#if 0

subroutine cce_initialize()     ! This subroutine can be merged btw XGC & GENE
  use sim_param                 
    implicit none
    integer :: ipsi

    ! TODO: don't check for allocation, use a logical instead
    ! GENE calls set_default_couplier while XGC does the below
    if (.not.allocated(cce_density)) then
      call cce_set_defaults

      open(unit=20, file=cce_input_file, status='old', action='read')
      READ(20, NML=coupling)
      close(unit=20)

      if (cce_side == 0) then
        cce_my_side = 'core'
        cce_other_side = 'edge'
      else
        cce_my_side = 'edge'
        cce_other_side = 'core'
      endif

      cce_surface_count = cce_last_surface - cce_first_surface + 1
      cce_density_node_count = cce_density_last_node - cce_density_first_node + 1


! To merge GENE & XGC, uncomment
#if 0
      if(cce_dpot_index0) then
             fld_name_XGC="dpot0"
              if(mype==0) write(*, '(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              if(mype==0) write(*, '(A)')"This coupling is supposed to have a shift in the toroidal direction"
              if(mype==0) write(*, '(A)')"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        else
                fld_name_XGC="dpot1"
        endif
#endif

      cce_field_node_count = cce_density_node_count
      cce_field_first_node = cce_density_first_node
      cce_field_last_node = cce_density_last_node

      if (cce_density_node_count < 0) then
        allocate(cce_surface_first_node(cce_all_surface_number))
        allocate(cce_surface_last_node(cce_all_surface_number))

        open(unit=20,file=cce_input_file, status='old',action='read')
        READ(20, NML=surfaces)
        close(unit=20)

        cce_density_first_node = cce_surface_first_node(cce_first_surface)
        cce_density_last_node = cce_surface_last_node(cce_last_surface)
        cce_density_node_count = cce_density_last_node - cce_density_first_node + 1

        cce_field_first_node = cce_surface_first_node(1)
        cce_field_last_node = cce_surface_last_node(cce_all_surface_number)
        cce_field_node_count = cce_field_last_node - cce_field_first_node + 1
      endif

      ! 1:cce_density_node_count
      allocate(cce_density(cce_density_first_node:cce_density_last_node))
!      allocate(cce_density(cce_density_first_node:cce_density_last_node,cce_first_surface:cce_last_surface))

      if (cce_field_model /= 0 .and. cce_comm_field_mode /= 0) then
!        allocate(cce_dpot0(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
 !       allocate(cce_dpot1(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
  !      allocate(cce_pot0(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
        allocate(cce_dpot0(cce_field_first_node:cce_field_last_node))
        allocate(cce_dpot1(cce_field_first_node:cce_field_last_node))
        allocate(cce_pot0(cce_field_first_node:cce_field_last_node))
      endif

      if (cce_first_surface_coupling < 1) cce_first_surface_coupling = cce_first_surface
      if (cce_last_surface_coupling < 1) cce_last_surface_coupling = cce_last_surface

      if ( cce_npsi > 0) then
        allocate(cce_varpi(cce_npsi))
        allocate(cce_psi(cce_npsi))

        open(unit=20,file=cce_input_file, status='old',action='read')
        READ(20, NML=varpi)
        close(unit=20)

        cce_first_surface_coupling = -1
        cce_last_surface_coupling = -1

        do ipsi = 1, cce_npsi
          if (cce_first_surface_coupling == -1 .and. cce_varpi(ipsi) > 0.0001) then
            cce_first_surface_coupling = ipsi
          endif
          if (cce_last_surface_coupling == -1 .and. cce_varpi(ipsi) > 0.9999) then
            cce_last_surface_coupling = ipsi
          endif
        enddo
        print *,"psi_min=",cce_first_surface_coupling,"psi_max=",cce_last_surface_coupling

        if(cce_side == 0) cce_varpi(:) = 1D0 - cce_varpi(:)
      endif
    endif
end subroutine cce_initialize
#endif


subroutine initialize_coupling_engines(node_number, nphi_total,offset_x, block_count)
 ! use communications, only: my_mpi_comm_world
  use sim_param, only: my_pespec, my_pew, my_pez, my_pev, mpi_comm_x, my_mpi_comm_world
  integer, intent(in) ::node_number, nphi_total,offset_x, block_count

end subroutine initialize_coupling_engines


subroutine check_coupler(gene_first, gene_nfs, gene_nnodes)
        integer :: gene_first, gene_nfs, gene_nnodes
        if(gene_first+1.ne. cce_first_surface) then
                print *, "Mistmatch on first surface"
        end if

        if(cce_last_surface.gt.gene_first+gene_nfs) then
                print *, "Mistmatch on last surface"
        end if

        if(gene_nnodes.ne. cce_node_number) then
                print *, "Mistmatch on number of nodes GENE ", gene_nnodes, " XGC ", cce_node_number
        end if
end subroutine check_coupler

subroutine send_density(density, iphi, nphi_total, block_count, block_start, block_end, comm)
    use sim_param, only: time, mype, itime

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
    character(4)::fld_name='data'

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

          call adios2_declare_io(dens_io,adios2obj,'density_coupling',adios_err)
          call adios2_define_variable(dens_id, dens_io,fld_name,&
               & adios2_type_dp, 2, gdims, goffset,&
               & ldims,adios2_constant_dims, adios_err)
          call adios2_open(engines(1), dens_io, trim(cce_folder)//'/'&
               & //'density.bp',adios2_mode_write,&
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
    use sim_param, only: mype
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


subroutine receive_field(data_block, block_start, block_end, block_count, nphi_total, comm)
    use sim_param, only:mype
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
                & "field_coupling", adios_err)
           call adios2_open(engines(2), coupling_io, &
                & trim(cce_folder)//'/'//'field.bp', &
                & adios2_mode_read, comm, adios_err)
           read_engine_av=.true.
           IF (MYPE.EQ.0)print *, 'created engine for field'
        endif

        call adios2_begin_step(engines(2),adios2_step_mode_read,adios_err)
        call adios2_inquire_variable(field_id,coupling_io,'dadat',adios_err)
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
     use sim_param, only: mype

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
     use sim_param, only: mype

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


subroutine dump_dt(dt)
     use sim_param, only: mype, cref, Lref, time, itime
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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! XGC #ifndef GENE_SIDE starts here !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! The GENE_SIDE macro is defined in XGC
subroutine cce_receive_density()
  use adios_read_mod
  use sim_param                 !  use sml_module

    logical :: ex
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    integer*8 :: buf_id, buf_size, total_size
    integer :: err
    integer*8 :: sel2, start2(2), count2(2)
    real*8, dimension(:),allocatable :: arrtmp
#ifdef ADIOS2
    type(adios2_variable) :: varid
    character(4):: fld_name='data'
#endif
!    call t_startf("CCE_RECEIVE_DENSITY")

    cce_density=0D0

    !if(cce_side.eq.1.and.cce_comm_density_mode.GT.1)then
    if(cce_comm_density_mode.eq.2.or.cce_comm_density_mode.eq.3)then
       allocate(arrtmp(cce_node_number))
       arrtmp=0D0
       start2(1) = 0
       count2(1) = cce_node_number
       start2(2) = sml_intpl_mype
       count2(2) = 1

#ifdef ADIOS2
       if (cce_step .eq. 1) then
          call adios2_declare_io(read_io, adios2obj, 'density_coupling',&
               & err)
          call adios2_open(read_engine,read_io,trim(cce_folder)//'/density.bp',&
               &adios2_mode_read, sml_intpl_comm, err)
       endif
       call adios2_begin_step(read_engine, adios2_step_mode_read, err)
       call adios2_inquire_variable(varid,read_io,fld_name, err)
       if (.not.varid%valid)then
          if (sml_mype.eq.0) print *, fld_name, ' variable not found'
       else
          if (sml_mype.eq.0) print *, fld_name, ' found'
       endif
       call adios2_set_selection(varid, 2, start2, count2, err)
       call adios2_get(read_engine,varid,arrtmp,adios2_mode_deferred,err)
       call adios2_end_step(read_engine, err)

#else
      ex=.false.

       write(cce_stepstr,'(I0.5)') cce_step
       write(planestr,'(I0.5)') sml_intpl_mype

#ifndef CCE_TEST_1

       cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//&
            &'_'//trim(cce_stepstr)//'.unlock'
       !    print *,'Wait unlock filename',cce_filename
       do while((.not.ex).and.(staging_read_method.eq.0))
          inquire(file=cce_filename,EXIST=ex)
       end do

       cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//&
            &'_'//trim(cce_stepstr)//'.bp'
       if (staging_read_method.gt.0) then
          cce_filename='density_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
          if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
          call adios_read_open (buf_id, cce_filename, staging_read_method,&
               &sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
       else
#ifdef CCE_DEBUG
          if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
          call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
       end if
       if(err/=0) then
          print *, 'coupling receive error: could not open file', cce_filename
          !stop
       endif

       call adios_selection_boundingbox(sel2, 2, start2, count2)
       call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
       call adios_perform_reads (buf_id, err)
       call adios_read_close (buf_id, err)

#else

       allocate(arrtmp(cce_node_number))
       arrtmp=0D0
       if(sml_intpl_mype.eq.0) arrtmp(90000)=1D0

       cce_filename=trim(cce_folder)//'/density_test_'&
            &//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
       ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',sml_intpl_comm,err)  !sml_comm_null,err)
       buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer

       ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
       ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
       ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
       ADIOS_WRITE_LBL(buf_id,'first_node',cce_first_node,err)
       ADIOS_WRITE_LBL(buf_id,'last_node',cce_last_node,err)
       ADIOS_WRITE_LBL(buf_id,'node_number',cce_node_number,err)
       ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
       ADIOS_WRITE_LBL(buf_id,'cce_model',cce_density_model,err)
       ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
       ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
       !actual data
       ADIOS_WRITE_LBL(buf_id,'data', arrtmp,err)

       ADIOS_CLOSE(buf_id,err)

#endif

#endif

       cce_density(cce_first_node:cce_last_node)=arrtmp(1:cce_node_number)

       !print *,trim(cce_my_side)//'density_'//trim(planestr)//'_'//trim(cce_stepstr)//'_R',arrtmp(1),arrtmp(cce_node_number)
       deallocate(arrtmp)

    endif
!    call t_stopf("CCE_RECEIVE_DENSITY")
end subroutine cce_receive_density


#ifndef GENE_SIDE
subroutine cce_process_density_GENE()
cce_density_step = cce_density_step + 1 
end subroutine cce_process_density_GENE

subroutine cce_process_density_XGC(density)
  use mpi
  use sim_param        !  use sml_module
  implicit none
  real(8), dimension(:), intent(inout) :: density

!  include 'mpif.h'
  real(8) :: alpha

  integer :: ipsi,ipsi0,ipsi1

  ! ADIOS
  integer(8) :: buf_id, buf_size, total_size
  integer :: err

  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
!  call t_startf("CCE_PROCESS_DENSITY")

  write(cce_stepstr,'(I0.5)') cce_step
  write(planestr,'(I0.5)') sml_intpl_mype

  !cce_step=cce_step+1

  select case (cce_density_model)
    case (-1)
      density(cce_first_node:cce_last_node)=0D0
    case (0)
      !Do nothing
    case(1)
      !Linear coupling
      if((cce_side.EQ.1).AND.(cce_first_surface.LT.cce_first_surface_coupling))then
        ipsi0=cce_surface_first_node(cce_first_surface)
        ipsi1=cce_surface_last_node(cce_first_surface_coupling)
        density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
      endif
      do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
        alpha=dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
        ipsi0=cce_surface_first_node(ipsi)
        ipsi1=cce_surface_last_node(ipsi)
        if(cce_side.EQ.1)alpha=1D0-alpha
        density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)+alpha*cce_density(ipsi0:ipsi1)
      enddo
      if((cce_side.EQ.0).AND.(cce_last_surface.GT.cce_last_surface_coupling))then
        ipsi0=cce_surface_first_node(cce_last_surface_coupling)
        ipsi1=cce_surface_last_node(cce_last_surface)
        density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
      endif
    case(2)
      !Average
      alpha=cce_alpha
      ipsi0=cce_first_node
      ipsi1=cce_last_node
      density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)+alpha*cce_density(ipsi0:ipsi1)
    case(3)
      ipsi0=cce_first_node
      ipsi1=cce_last_node
      density(ipsi0:ipsi1)=1D0*density(ipsi0:ipsi1)
    case(4)
      alpha=cce_alpha
      ipsi0=cce_first_node
      ipsi1=cce_last_node
      density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)
    case(5)
      if((cce_side.EQ.1))then
        ipsi0=cce_surface_first_node(cce_first_surface)
        ipsi1=cce_surface_last_node(cce_last_surface_coupling)
        density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
      endif
      if((cce_side.EQ.0))then
        ipsi0=cce_surface_first_node(cce_first_surface_coupling)
        ipsi1=cce_surface_last_node(cce_last_surface)
        density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
      endif
#ifdef XGC_COUPLING_CORE_EDGE_VARPI2
    case(6)
      ipsi0=cce_surface_first_node(cce_first_surface) !cce_first_node
      ipsi1=cce_surface_last_node(cce_last_surface)   !cce_last_node
      !print *,'case(6) density(ipsi0:ipsi1)',ipsi0,ipsi1
      density(ipsi0:ipsi1)=density(ipsi0:ipsi1)+cce_density(ipsi0:ipsi1)
#endif
    case default
      print *,'Unknown coupling model'
      stop
  end select
if(.false.)then
  cce_filename=trim(cce_folder)//'/'//trim(cce_my_side)//'after_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#ifdef CCE_DEBUG
  if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif

  call adios_open(buf_id,trim('coupling')//char(0),trim(cce_filename)//char(0),'w'//char(0), MPI_COMM_SELF, err)
  buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer
  call adios_group_size(buf_id,buf_size,total_size,err)
  call adios_write(buf_id,trim('nphi')//char(0),sml_nphi_total,err)
  call adios_write(buf_id,trim('iphi')//char(0),sml_intpl_mype,err)
  call adios_write(buf_id,trim('first_node')//char(0),cce_first_node,err)
  call adios_write(buf_id,trim('last_node')//char(0),cce_last_node,err)
  call adios_write(buf_id,trim('node_number')//char(0),cce_node_number,err)
  call adios_write(buf_id,trim('cce_side')//char(0),cce_side,err)
  call adios_write(buf_id,trim('cce_model')//char(0),cce_density_model,err)
  call adios_write(buf_id,trim('time')//char(0),sml_time,err)
  call adios_write(buf_id,trim('step')//char(0),sml_gstep,err)
  call adios_write(buf_id,trim('data')//char(0), density(cce_first_node:cce_last_node),err)   !actual data

  call adios_close(buf_id,err)
endif

  cce_step=cce_step+1

!  call t_stopf("CCE_PROCESS_DENSITY")

end subroutine cce_process_density_XGC

 subroutine cce_send_field(dpot0,dpot1,pot0,flag_pot0)
    use sim_param
    implicit none
!    include 'mpif.h'

    real(8), dimension(:), intent(in) :: dpot0,dpot1,pot0
    !real*8, dimension(:,:), intent(in) :: dpot
    integer, intent(in) :: flag_pot0

    character(5)::fld_name, cce_stepstr, planestr
    integer :: myrank, mysize, ierr
    character(512) :: cce_filename

    ! ADIOS
    integer(8) :: buf_id, buf_size, total_size
    integer :: err,mode
#ifdef ADIOS2
!    type(adios2_variable) :: varid
    integer(8), dimension(2) :: gdims, goffset, ldims
#endif
    real(8), dimension(:),allocatable :: arrtmp
    !call t_startf("CCE_SEND_FIELD")

    if(cce_side.eq.1.and.cce_comm_field_mode.GT.0)then

       write(cce_stepstr,'(I0.5)') cce_field_step
       write(planestr,'(I0.5)') sml_intpl_mype


#ifdef ADIOS2
       if (cce_field_step == 1) then
          allocate(field_out(cce_field_node_number))
          gdims(1) = cce_field_node_number
          gdims(2) = sml_nphi_total
          goffset(1) = 0
          goffset(2) = sml_intpl_mype
          ldims(1) = cce_field_node_number
          ldims(2) = 1
          call adios2_declare_io(send_io, adios2obj, 'field_coupling', err)
          call adios2_define_variable(field_id, send_io, 'dadat',  &
               & adios2_type_dp, 2, gdims, goffset,&
               & ldims, adios2_constant_dims, err)
          call adios2_open(send_engine, send_io, trim(cce_folder) // '/'&
               & //'field.bp', adios2_mode_write,&
               & sml_intpl_comm, err)
          call MPI_Comm_rank(sml_intpl_comm, myrank, ierr);
          call MPI_Comm_size(sml_intpl_comm, mysize, ierr);
          print *, sml_mype, 'opens as rank ', myrank,' of size', mysize
       endif
#endif
       if (.not.allocated(field_out)) allocate(field_out(cce_field_node_number))

       if (cce_dpot_index0) then
          fld_name='dpot0'
          field_out(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
          if(cce_nzeromode) field_out(1:cce_field_node_number)=field_out(1:cce_field_node_number)+&
               &pot0(cce_field_first_node:cce_field_last_node)
       else  !else of "if(cce_dpot_index0)then"
          fld_name='dpot1'
          field_out(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
          if(cce_nzeromode) field_out(1:cce_field_node_number)=field_out(1:cce_field_node_number)+&
               pot0(cce_field_first_node:cce_field_last_node)
       endif


#ifdef ADIOS2
       if (sml_mype.eq.0) print *, 'sending field as ', fld_name, kind(field_out)
       call adios2_begin_step(send_engine, adios2_step_mode_append, err)
       call adios2_put(send_engine, field_id, field_out, err)
       call adios2_end_step(send_engine, err)
#else
       if (staging_read_method.gt.0) then
          cce_filename= trim(fld_name)//'_'//trim(cce_my_side)//'.bp'
       else
          cce_filename=trim(cce_folder)//'/'//trim(fld_name)//&
               &'_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
       endif
#ifdef CCE_DEBUG
       if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif

      call adios_open(buf_id,trim('coupling')//char(0),trim(cce_filename)//char(0),'w'//char(0), sml_intpl_comm, err)
      buf_size= 4*8 + 8 + 8*cce_field_node_number  + 100 !last 100 is buffer
      call adios_group_size(buf_id,buf_size,total_size,err)
      call adios_write(buf_id,trim('nphi')//char(0),sml_nphi_total,err)
      call adios_write(buf_id,trim('iphi')//char(0),sml_intpl_mype,err)
      call adios_write(buf_id,trim('first_node')//char(0),cce_field_first_node,err)
      call adios_write(buf_id,trim('last_node')//char(0),cce_field_last_node,err)
      call adios_write(buf_id,trim('node_number')//char(0),cce_field_node_number,err)
      call adios_write(buf_id,trim('cce_side')//char(0),cce_side,err)
      call adios_write(buf_id,trim('cce_model')//char(0),cce_field_model,err)
      call adios_write(buf_id,trim('time')//char(0),sml_time,err)
      call adios_write(buf_id,trim('step')//char(0),sml_gstep,err)
      call adios_write(buf_id,trim('data')//char(0),arrtmp,err)
      call adios_close(buf_id,err)

       !Create an unlock file: one file per step
       call mpi_barrier(sml_intpl_comm, err)
       if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
          cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//&
               &'_'//trim(cce_stepstr)//'.unlock'
          open(20, file=cce_filename, status="new", action="write")
          close(20)
       endif
#endif
    endif


!    call t_stopf("CCE_SEND_FIELD")
  end subroutine cce_send_field








#if 0
subroutine cce_send_field(pot0, dpot0, dpot1, flag_pot0)
  use sim_param
  implicit none

  real(8), dimension(:), intent(in) :: pot0, dpot0, dpot1
  integer, intent(in) :: flag_pot0

  character(5)::cce_stepstr,planestr, fld_name
  character(512) :: cce_filename

  ! ADIOS
  integer(8) :: buf_id, buf_size, total_size
  integer :: err,mode

  real(8), dimension(:),allocatable :: arrtmp
!  call t_startf("CCE_SEND_FIELD")

  !cce_density=0D0

  if(cce_side.eq.1.and.cce_comm_field_mode.GT.0)then

    write(cce_stepstr,'(I0.5)') cce_field_step
    write(planestr,'(I0.5)') sml_intpl_mype

#ifdef ADIOS2
       if (cce_field_step == 1) then
          allocate(field_out(cce_field_node_number))
          gdims(1) = cce_field_node_number
          gdims(2) = sml_nphi_total
          goffset(1) = 0
          goffset(2) = sml_intpl_mype
          ldims(1) = cce_field_node_number
          ldims(2) = 1
          call adios2_declare_io(send_io, adios2obj, 'field_coupling', err)
          call adios2_define_variable(field_id, send_io, 'dadat',  &
               & adios2_type_dp, 2, gdims, goffset,&
               & ldims, adios2_constant_dims, err)
          call adios2_open(send_engine, send_io, trim(cce_folder) // '/'&
               & //'field.bp', adios2_mode_write,&
               & sml_intpl_comm, err)
          call MPI_Comm_rank(sml_intpl_comm, myrank, ierr);
          call MPI_Comm_size(sml_intpl_comm, mysize, ierr);
          print *, sml_mype, 'opens as rank ', myrank,' of size', mysize

       endif
#endif
       if (.not.allocated(field_out)) allocate(field_out(cce_field_node_number))

       if (cce_dpot_index0) then
          fld_name='dpot0'
          field_out(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
          if(cce_nzeromode) field_out(1:cce_field_node_number)=field_out(1:cce_field_node_number)+&
               pot0(cce_field_first_node:cce_field_last_node)
       endif


#ifdef ADIOS2
       if (sml_mype.eq.0) print *, 'sending field as ', fld_name, kind(field_out)
       call adios2_begin_step(send_engine, adios2_step_mode_append, err)
       call adios2_put(send_engine, field_id, field_out, err)
       call adios2_end_step(send_engine, err)
#else
       if (staging_read_method.gt.0) then
          cce_filename= trim(fld_name)//'_'//trim(cce_my_side)//'.bp'
       else
          cce_filename=trim(cce_folder)//'/'//trim(fld_name)//&
               &'_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
       endif
#ifdef CCE_DEBUG
       if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif

      call adios_open(buf_id,trim('coupling')//char(0),trim(cce_filename)//char(0),'w'//char(0), sml_intpl_comm, err)
      buf_size= 4*8 + 8 + 8*cce_field_node_number  + 100 !last 100 is buffer
      call adios_group_size(buf_id,buf_size,total_size,err)
      call adios_write(buf_id,trim('nphi')//char(0),sml_nphi_total,err)
      call adios_write(buf_id,trim('iphi')//char(0),sml_intpl_mype,err)
      call adios_write(buf_id,trim('first_node')//char(0),cce_field_first_node,err)
      call adios_write(buf_id,trim('last_node')//char(0),cce_field_last_node,err)
      call adios_write(buf_id,trim('node_number')//char(0),cce_field_node_number,err)
      call adios_write(buf_id,trim('cce_side')//char(0),cce_side,err)
      call adios_write(buf_id,trim('cce_model')//char(0),cce_field_model,err)
      call adios_write(buf_id,trim('time')//char(0),sml_time,err)
      call adios_write(buf_id,trim('step')//char(0),sml_gstep,err)
      call adios_write(buf_id,trim('data')//char(0),arrtmp,err)
      call adios_close(buf_id,err)

    !Create an unlock file: one file per step
    call mpi_barrier(sml_intpl_comm, err)
       if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
          cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//&
               &'_'//trim(cce_stepstr)//'.unlock'
          open(20, file=cce_filename, status="new", action="write")
          close(20)
       endif
#endif
    endif


!    call t_stopf("CCE_SEND_FIELD")
  end subroutine cce_send_field

#endif





subroutine cce_process_field_XGC(dpot0, dpot1, pot0, flag_pot0)
    implicit none

    integer, intent(in) :: flag_pot0
    real(8), dimension(:), intent(inout) :: pot0
    real(8), dimension(:), intent(inout) :: dpot0
    real(8), dimension(:), intent(inout) :: dpot1

  !real*8 :: alpha
 ! call t_startf("CCE_PROCESS_FIELD")
  !cce_bcast_dpot=.false.
  !cce_bcast_pot0=.false.

  select case (cce_field_model)
    case (0)
      !Do nothing
    case(1)
      cce_bcast_dpot=.true.
      dpot0(:)=cce_dpot0(:)
      dpot1(:)=cce_dpot1(:)
      if(flag_pot0.eq.0)then
        cce_bcast_pot0=.true.
        pot0(:)=cce_pot0(:)
      endif
    case default
      print *,'Unknown coupling model'
      stop
  end select

  cce_field_step=cce_field_step+1
!  call t_stopf("CCE_PROCESS_FIELD")

end subroutine cce_process_field_XGC
           



#endif
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! XGC #ifndef GENE_SIDE ends here !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cce_process_field_GENE()
    cce_field_step = cce_field_step + 1
  end subroutine cce_process_field_GENE



#ifndef GENE_SIDE                      
subroutine cce_varpi_grid(rho_ff)
      real (8), dimension(:), intent(inout) :: rho_ff
      integer :: ipsi,ipsi0,ipsi1
      real (8) :: varpi

      call cce_initialize()
      if(cce_density_model.eq.6)then

        if(cce_npsi>0)then
          if((cce_side.EQ.1).AND.(cce_first_surface.LT.cce_first_surface_coupling))then
            ipsi0=cce_surface_first_node(cce_first_surface)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling-1)
            rho_ff(ipsi0:ipsi1)=0D0!rho_ff(ipsi0:ipsi1)
          endif
          do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
            ipsi0=cce_surface_first_node(ipsi)
            ipsi1=cce_surface_last_node(ipsi)
            rho_ff(ipsi0:ipsi1)=cce_varpi(ipsi)*rho_ff(ipsi0:ipsi1)
          enddo
          if((cce_side.EQ.0).AND.(cce_last_surface.GT.cce_last_surface_coupling))then
            ipsi0=cce_surface_first_node(cce_last_surface_coupling+1)
            ipsi1=cce_surface_last_node(cce_last_surface)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
          endif

        else
           !if(ipsi.le.cce_first_surface_coupling_axis)then
           !   varpi=0D0
           !elseif(ipsi.le.cce_last_surface_coupling_axis)then
           !   varpi=1D0-dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
           !elseif (ipsi.le.cce_first_surface_coupling) then
           !   varpi=1D0
           !elseif (ipsi.gt.cce_last_surface_coupling)then
           !   varpi=0D0
           !else
           !   varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
           !endif

          !Linear weight
          if(cce_side.EQ.0)then !core
            !1. Zero near axis
            ipsi0=cce_surface_first_node(cce_first_surface)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling_axis)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
            !2. Linear increasing weight 
            do ipsi=cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
              varpi=dble(ipsi-cce_first_surface_coupling_axis)/&
                  &dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            !3. Unity in the middle

            !4. Linear decreasinging weight
            do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
              varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              !varpi=1D0-varpi
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            !5. Zero near edge
            ipsi0=cce_surface_first_node(cce_last_surface_coupling)
            ipsi1=cce_surface_last_node(cce_last_surface)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
          elseif(cce_side.EQ.1)then !edge
            ipsi0=cce_surface_first_node(cce_last_surface_coupling_axis)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling)
            rho_ff(ipsi0:ipsi1)=0D0!rho_ff(ipsi0:ipsi1)

            do ipsi=cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
              varpi=dble(ipsi-cce_first_surface_coupling_axis)/&
                   &dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              varpi=1D0-varpi
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
              varpi=dble(ipsi-cce_first_surface_coupling)/&
                   &dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
          endif
        endif
      endif
    end subroutine
#else
  !GENE SIDE
  function cce_varpi_grid(ipsi) result(varpi)
    real :: varpi
    integer :: ipsi

!      call cce_initialize()
    if(cce_density_model.eq.6)then
      if(cce_npsi>0)then
          stop('Not implemented');
      else
        !Linear weight
          if (cce_side.EQ.0)then
             if(ipsi.le.cce_first_surface_coupling_axis)then
                varpi=0D0
             elseif(ipsi.le.cce_last_surface_coupling_axis)then
                varpi=dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
             elseif (ipsi.le.cce_first_surface_coupling) then
                varpi=1D0
             elseif (ipsi.gt.cce_last_surface_coupling)then
                varpi=0D0
             else
                varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
             endif
           else
             stop('GENE is for the core, put cce_side=0');
          endif
        endif
      endif
    end function
#endif



!!!!!!!!!!!!!!!!
!! USELESS - The following arent used
!!!!!!!!!!!!!!!!

#if 0
subroutine cce_unique_id(namestr, method, side)        ! This should replace cce_get_filename and should use hdf5
 character(len=*), intent(in) :: namestr
 character(len=*), intent(in) :: side
 integer, intent(in) :: method

 !character(len=8), parameter :: filename=trim(trim(cce_folder), trim(trim(namestr), trim(side)), trim(cce_stepstr))      ! file name
 character(len=4), parameter :: dsetname = "gene_dset"          ! dataset name

 write(cce_lock_filename, "(A,'_', A)")
 integer(hid_t) :: file_id      ! file identifier
 integer(hid_t) :: dset_id      ! dataset identifier
 integer(hid_t) :: dspace_id    ! dataspace identifier

 integer :: error               ! error flag

 call h5open_f(error)                                           ! initialize hdf5 interface
 call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)     ! create new default file

 call h5fclose_f(dset_id, error)                                ! close file
 call h5close_f(error)                                          ! close interface

end subroutine cce_unique_id
#endif

subroutine cce_destroy()
  if(allocated(cce_surface_first_node)) deallocate(cce_surface_first_node)
  if(allocated(cce_surface_last_node)) deallocate(cce_surface_last_node)
  if(allocated(cce_density)) deallocate(cce_density)
  if(allocated(cce_dpot0)) deallocate(cce_dpot0)
  if(allocated(cce_dpot1)) deallocate(cce_dpot1)
  if(allocated(cce_pot0)) deallocate(cce_pot0)
end subroutine cce_destroy


subroutine cce_advance_density()

end subroutine cce_advance_density


subroutine cce_advance_field()

end subroutine cce_advance_field
        

subroutine cce_preprocess_density(density, plane_id)
    real(8), dimension(:), intent(inout) :: density
    integer, intent(in) :: plane_id
    real(8) :: alpha

    integer :: ipsi
    integer :: ipsi0
    integer :: ipsi1

    character(512) :: cce_filename
end subroutine cce_preprocess_density


subroutine cce_preprocess_field(dpot0, dpot1, pot0, flag_pot0)
    implicit none
    real(8), dimension(:), intent(inout) :: pot0
    real(8), dimension(:), intent(inout) :: dpot0
    real(8), dimension(:), intent(inout) :: dpot1
    integer, intent(in) :: flag_pot0
end subroutine cce_preprocess_field
  
#ifndef GENE_SIDE       

! 1D XGC version of the send density function
subroutine cce_send_density(density)
  use sim_param 
  implicit none
    real(8), dimension(:), intent(in) :: density
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename

    ! ADIOS
    integer(8) :: buf_id, buf_size, total_size
    integer :: err

    real(8), dimension(:),allocatable :: arrtmp
!    call t_startf("CCE_SEND_DENSITY")

    cce_density=0D0

    if(cce_comm_density_mode.eq.1.or.cce_comm_density_mode.eq.2)then

    write(cce_stepstr,'(I0.5)') cce_step
    write(planestr,'(I0.5)') sml_intpl_mype

    allocate(arrtmp(cce_node_number))
    arrtmp(1:cce_node_number)=density(cce_first_node:cce_last_node)
    !cce_density(cce_first_node:cce_last_node)=density(cce_first_node:cce_last_node)
    if (staging_read_method.gt.0) then
    cce_filename='density_'//trim(cce_my_side)//'.bp'
    else
    cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
    endif

    !print *,'Send
    !filename',cce_filename

   !TODO: Open the file and write
    !the density for each node
    !point
#ifdef CCE_DEBUG
    if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif
    
    call adios_open(buf_id,trim('coupling')//char(0),trim(cce_filename)//char(0),'w'//char(0), sml_intpl_comm, err)
    buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer
    call adios_group_size(buf_id,buf_size,total_size,err)
    call adios_write(buf_id,trim('nphi')//char(0),sml_nphi_total,err)
    call adios_write(buf_id,trim('iphi')//char(0),sml_intpl_mype,err)
    call adios_write(buf_id,trim('first_node')//char(0),cce_first_node,err)
    call adios_write(buf_id,trim('last_node')//char(0),cce_last_node,err)
    call adios_write(buf_id,trim('node_number')//char(0),cce_node_number,err)
    call adios_write(buf_id,trim('cce_side')//char(0),cce_side,err)
    call adios_write(buf_id,trim('cce_model')//char(0),cce_density_model,err)
    call adios_write(buf_id,trim('time')//char(0),sml_time,err)
    call adios_write(buf_id,trim('step')//char(0),sml_gstep,err)
    call adios_write(buf_id,trim('data')//char(0),arrtmp,err)   !actual data

    call adios_close(buf_id,err)

!    print *,trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'_W',arrtmp(2),arrtmp(cce_node_number)

    deallocate(arrtmp)

    !Create an unlock file: one file per step
    if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
      cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)
    endif

  endif
!  call t_stopf("CCE_SEND_DENSITY")

end subroutine cce_send_density


#endif


! updated this
subroutine cce_receive_field_XGC(flag_pot0)
  use mpi
  use sim_param
  implicit none
!  include 'mpif.h'

  integer, intent(in) :: flag_pot0

  logical :: ex
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
  integer*8 :: buf_id
  integer :: err
  integer*8 :: sel2, start2(2), count2(2)
  real*8, dimension(:),allocatable :: arrtmp
  !! call t_startf("CCE_RECEIVE_FIELD")
  
  cce_dpot0=0D0
  cce_dpot1=0D0
  cce_pot0=0D0

  if(cce_side.eq.0.and.cce_comm_field_mode.GT.1)then

    ex=.false.

    write(cce_stepstr,'(I0.5)') cce_field_step
    write(planestr,'(I0.5)') sml_intpl_mype

    cce_filename=trim(cce_folder)//'/field_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.unlock'
    do while((.not.ex).and.(staging_read_method.eq.0))
      inquire(file=cce_filename,EXIST=ex)
    end do

    cce_filename=trim(cce_folder)//'/dpot0_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
    if (staging_read_method.gt.0) then
      cce_filename='dpot0_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
    else
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
    end if
    if(err/=0) then
      print *, 'coupling receive error: could not open file', cce_filename
      stop
    endif
    allocate(arrtmp(cce_field_node_number))
    arrtmp=0D0

    start2(1) = 0
    count2(1) = cce_field_node_number
    start2(2) = sml_intpl_mype
    count2(2) = 1
    call adios_selection_boundingbox(sel2, 2, start2, count2)
    call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)
    cce_dpot0(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)

    cce_filename=trim(cce_folder)//'/dpot1_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
    if (staging_read_method.gt.0) then
      cce_filename='dpot1_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
    else
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
    end if
    if(err/=0) then
      print *, 'coupling receive error: could not open file', cce_filename
      stop
    endif
    arrtmp=0D0
    call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)
    cce_dpot1(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)

    if(flag_pot0.eq.0)then
      arrtmp=0D0
      cce_filename=trim(cce_folder)//'/pot0_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
      if (staging_read_method.gt.0) then
        cce_filename='pot0_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
        if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
        call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
      else 
#ifdef CCE_DEBUG
        if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
        call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
      end if
      if(err/=0) then
        print *, 'coupling receive error: could not open file', cce_filename
        stop
      endif
      call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
      cce_pot0(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)
    endif

    deallocate(arrtmp)
  endif

 ! call t_stopf("CCE_RECEIVE_FIELD")

end subroutine cce_receive_field_XGC



#ifdef ADIOS
  subroutine get_cce_filename(namestr, method, side, lockname)
    character(len=*), intent(in) :: namestr
    character(len=*), intent(in) :: side
    integer, intent(in) :: method
    character(len=*), intent(in), optional :: lockname

    character(len=256) :: shortbase
    character(len=256) :: longbase
    character(len=5) :: cce_stepstr

    write(cce_stepstr,'(I0.5)') cce_step

    if (present(lockname)) then
        write(shortbase, "(A, '_', A)") trim(lockname), trim(side)
        write(longbase, "(A, '/', A, '_', A)") trim(cce_folder), trim(shortbase), trim(cce_stepstr)
	write(cce_lock_filename, "(A, '.unlock')")  trim(longbase)
    end if

    write(shortbase, "(A, '_', A)") trim(namestr), trim(side)
    write(longbase, "(A, '/', A, '_', A)") trim(cce_folder), trim(shortbase), trim(cce_stepstr)
    if (.not.present(lockname)) then
       write(cce_lock_filename, "(A, '.unlock')") trim(longbase)
    end if

    if (method.eq.ADIOS_READ_METHOD_BP) then
       write(cce_filename, "(A, '.bp')") trim(longbase)
    else
      write(cce_filename, "(A, '.bp')") trim(shortbase)
    end if

  end subroutine get_cce_filename

  subroutine cce_couple_open(adios_handle, namestr, method, side, comm, adios_err, lockname)
    character(len=*), intent(in) :: namestr
    character(len=*), intent(in) :: side
    integer, intent(in) :: method
    integer(8), intent(inout) :: adios_handle
    integer, intent(in) :: comm
    integer, intent(inout) :: adios_err
    character(len=*), intent(in), optional :: lockname
    logical :: ex
    real(4) :: timeout = -1.0

    if (present(lockname)) then
       call get_cce_filename(namestr, method, side, lockname=lockname)
    else
       call get_cce_filename(namestr, method, side)
    end if


    if (method .eq. ADIOS_READ_METHOD_BP) then
       !print *,'Wait unlock filename',cce_lock_filename
       !print *,'Wait filename',cce_filename
       ex = .false.
       do while(.NOT.ex)
          inquire(file=trim(cce_lock_filename), EXIST=ex)
       end do
       call adios_read_open_file(adios_handle, cce_filename, method, comm, adios_err)

    else
      call adios_read_open(adios_handle, cce_filename, method, comm, ADIOS_LOCKMODE_ALL, timeout, adios_err)
    end if

    if (adios_err.ne.0) then
       print *, 'coupling receive error: could not open file', cce_filename
       stop
    endif
  end subroutine cce_couple_open
#endif




end module new_coupling
