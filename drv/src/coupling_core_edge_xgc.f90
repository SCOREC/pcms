module coupling_core_edge
  implicit none
  private

  ! program metadata
  character(:), allocatable :: cce_input_file
  character(:), allocatable :: cce_output_dir
  logical :: cce_initialized

  ! simulation metadata
  integer :: cce_side ! 0:core (GENE/XGC), 1: edge(XGC)
  integer :: cce_density_model

  character(:), allocatable :: cce_adios_group
  character(5) :: cce_my_side
  character(5) :: cce_other_side

  logical :: cce_bcast_dpot
  logical :: cce_bcast_pot0

  integer, dimension(:), allocatable :: cce_surface_first_node
  integer, dimension(:), allocatable :: cce_surface_last_node

  integer :: cce_first_surface
  integer :: cce_last_surface
  integer :: cce_surface_number
  integer :: cce_all_surface_number

  integer :: cce_first_node
  integer :: cce_last_node
  integer :: cce_node_count

  integer :: cce_step
  integer :: cce_field_step
  integer :: cce_field_model

  integer :: cce_comm_density_mode
  integer :: cce_comm_field_mode

  integer :: cce_field_first_node
  integer :: cce_field_last_node
  integer :: cce_field_node_number

  integer :: cce_first_surface_field
  integer :: cce_last_surface_field

  integer :: cce_first_surface_coupling
  integer :: cce_last_surface_coupling

  integer :: cce_first_surface_coupling_axis
  integer :: cce_last_surface_coupling_axis

  ! physical data (fields)
  real, dimension(:), allocatable :: cce_density
  real, dimension(:), allocatable :: cce_pot0
  real, dimension(:), allocatable :: cce_dpot0
  real, dimension(:), allocatable :: cce_dpot1

  !
  real(8) :: cce_alpha
  real(8) :: cce_dt ! only in gene version

  !
  integer :: cce_npsi
  real(8), dimension(:), allocatable :: cce_varpi
  real(8), dimension(:), allocatable :: cce_psi

  ! namelists
  namelist /coupling/ cce_side, &
    cce_density_model, &
    cce_output_dir, &
    cce_first_surface, &
    cce_last_surface, &
    cce_first_node, &
    cce_last_node, &
    cce_comm_density_mode, &
    cce_all_surface_number, &
    cce_alpha, &
    cce_field_model, &
    cce_comm_field_mode, &
    cce_first_surface_coupling, &
    cce_last_surface_coupling, &
    cce_first_surface_field, &
    cce_last_surface_field, &
    cce_npsi, &
    cce_first_surface_coupling_axis, &
    cce_last_surface_coupling_axis

  namelist /surfaces/ cce_surface_first_node, cce_surface_last_node
  namelist /varpi/ cce_varpi, cce_psi

contains

  subroutine cce_set_defaults
    implicit none
    cce_filename = 'coupling.in'
    ccs_adios_group = 'coupling'
    cce_alpha = 0.5D0
    cce_density_model = 0
    cce_step = 1 !In case of restart, one might want to change this
    cce_field_step = 1
    cce_first_node = 10
    cce_last_node = 0
    cce_comm_density_mode = 2
    cce_all_surface_number = 1
    cce_first_surface_coupling = -1
    cce_last_surface_coupling = -1
    cce_field_model = 0
    cce_comm_field_mode = 0
    cce_npsi = -1
  end subroutine cce_set_defaults

  !TODO clean this and make the minimum common
  subroutine cce_initialize
    implicit none
    integer :: ipsi

    ! TODO: don't check for allocation, use a logical instead
    ! if (.not. cce_initialized) then
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

      cce_surface_number = cce_last_surface - cce_first_surface + 1
      cce_node_count = cce_last_node - cce_first_node + 1

      cce_field_node_number = cce_node_count
      cce_field_first_node = cce_first_node
      cce_field_last_node = cce_last_node

      if (cce_node_count < 0) then
        allocate(cce_surface_first_node(cce_all_surface_number))
        allocate(cce_surface_last_node(cce_all_surface_number))

        open(unit=20,file=cce_input_file, status='old',action='read')
        READ(20, NML=surfaces)
        close(unit=20)

        cce_first_node = cce_surface_first_node(cce_first_surface)
        cce_last_node = cce_surface_last_node(cce_last_surface)
        cce_node_count = cce_last_node - cce_first_node + 1

        cce_field_first_node = cce_surface_first_node(1)
        cce_field_last_node = cce_surface_last_node(cce_all_surface_number)
        cce_field_node_number = cce_field_last_node - cce_field_first_node + 1
      endif

      ! 1:cce_node_count
      allocate(cce_density(cce_first_node:cce_last_node))

      if (cce_field_model /= 0 .and. cce_comm_field_mode /= 0) then
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

  subroutine cce_destroy
    if(allocated(cce_surface_first_node)) deallocate(cce_surface_first_node)
    if(allocated(cce_surface_last_node)) deallocate(cce_surface_last_node)
    if(allocated(cce_density)) deallocate(cce_density)
    if(allocated(cce_dpot0)) deallocate(cce_dpot0)
    if(allocated(cce_dpot1)) deallocate(cce_dpot1)
    if(allocated(cce_pot0)) deallocate(cce_pot0)
  end subroutine cce_destroy

  ! xgc version
  subroutine cce_send_density(density)
    use adios_mod
    use sml_module, only: sml_intpl_mype, sml_nphi_total, sml_time, sml_gstep

    ! program data
    character(:), allocatable :: cce_filename

    ! simulation data
    character(5) :: cce_stepstr
    character(5) :: planestr

    ! simulation data (adios)
    integer(8) :: buf_id,
    integer(8) :: buf_size
    integer(8) :: total_size
    integer :: err

    ! physical data
    real(8), dimension(:), intent(in) :: density

    if (cce_comm_density_mode == 1 .or. cce_comm_density_mode == 2) then

      write(planestr,'(I0.5)') sml_intpl_mype
      write(cce_stepstr,'(I0.5)') cce_step

      cce_density(:) = density(cce_first_node:cce_last_node)

      cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'

      buf_size = (40 + (8 * cce_node_count))

      ! open buffer and set size
      call adios_open(buf_id,cce_adios_group,cce_filename,'w',MPI_COMM_SELF,err)
      call adios_group_size(buf_id,buf_size,total_size,err)

      ! write coupling header
      call adios_write(buf_id,'nphi',sml_nphi_total,err)
      call adios_write(buf_id,'iphi',sml_intpl_mype,err)
      call adios_write(buf_id,'first_node',cce_first_node,err)
      call adios_write(buf_id,'last_node',cce_last_node,err)
      call adios_write(buf_id,'node_number',cce_node_count,err)
      call adios_write(buf_id,'cce_side',cce_side,err)
      call adios_write(buf_id,'cce_model',cce_density_model,err)
      call adios_write(buf_id,'time',sml_time,err)
      call adios_write(buf_id,'step',sml_gstep,err)

      ! write actual data
      call adios_write(buf_id,'data',cce_density,err)

      ! close the buffer
      call adios_close(buf_id,err)

      !Create an unlock file

      cce_filename = trim(cce_output_dir)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)
    endif
  end subroutine cce_send_density

  ! gene version, need to unify
  subroutine cce_send_density(density,iphi,nphi_total,myli0,myli1,myli2,comm)
    use mpi
    use par_mod, only: time
    use par_other, only: itime

    integer, intent(in) :: comm
    integer, intent(in) :: iphi
    integer, intent(in) :: nphi_total
    integer, intent(in) :: myli1
    integer, intent(in) :: myli2
    integer, intent(in) :: myli0
    real, dimension(myli1:myli2), intent(in) :: density

    character(5) :: cce_stepstr
    character(5) :: planestr
    character(:), allocatable :: cce_filename

    ! adios
    integer(8) :: adios_handle
    integer(8) :: adios_groupsize
    integer(8) :: adios_totalsize
    integer :: err

     if (cce_comm_density_mode == 1 .or. cce_comm_density_mode == 2) then

      write(cce_stepstr,'(I0.5)') cce_step
      write(planestr,'(I0.5)') iphi

      cce_filename = trim(cce_output_dir)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//&
        & '_'//trim(cce_stepstr)//'.bp'

      call adios_open(adios_handle,cce_adios_group,cce_filename,'w',comm,err)
      adios_groupsize = 0
      adios_groupsize = adios_groupsize + 4_8 &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 8_8 * (block_count) * (nphi_total) &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 4_8 &
        + 8_8 &
        + 4_8
      call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, err)
      call adios_write (adios_handle, "block", block_count, err)
      call adios_write (adios_handle, "offset_x", block_start, err)
      call adios_write (adios_handle, "nnodes", cce_node_number, err)
      call adios_write (adios_handle, "nphi_total", nphi_total, err)
      call adios_write (adios_handle, "iphi", iphi, err)
      call adios_write (adios_handle, "data", tmp, err)
      call adios_write (adios_handle, "start_fs", cce_first_surface, err)
      call adios_write (adios_handle, "first_node", cce_first_node, err)
      call adios_write (adios_handle, "last_node", cce_last_node, err)
      call adios_write (adios_handle, "node_number", cce_node_number, err)
      call adios_write (adios_handle, "cce_side", cce_side, err)
      call adios_write (adios_handle, "cce_density_model", cce_density_model, err)
      call adios_write (adios_handle, "time", time, err)
      call adios_write (adios_handle, "step", itime, err)

      call adios_close(adios_handle,err)

      !Create an unlock file
      cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//&
        & '_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)

    endif

  end subroutine cce_send_density

  subroutine cce_process_density()
    ! TODO: this definitely should not be here
    cce_step = cce_step + 1
  end subroutine cce_process_density


#endif

#ifndef GENE_SIDE
  subroutine cce_receive_density()
    use adios_read_mod
    use sml_module

    include 'mpif.h'

    logical :: ex
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    integer*8 :: buf_id
    integer :: adios_read_method = ADIOS_READ_METHOD_BP, err
    integer*8 :: sel1=0
    real*8, dimension(:),allocatable :: arrtmp
    integer :: stat
    call t_startf("CCE_RECEIVE_DENSITY")
    cce_density=0D0

    if(cce_comm_density_mode.eq.2.or.cce_comm_density_mode.eq.3)then
      ex=.false.

      write(cce_stepstr,'(I0.5)') cce_step
      write(planestr,'(I0.5)') sml_intpl_mype

      cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
      err=-1
      do while(err/=0)
        do while(.NOT.ex)
          inquire(file=cce_filename,EXIST=ex)
        end do
#ifndef CCE_DEBUG
        open(unit=1234, iostat=stat, file=cce_filename, status='old')
        if (stat == 0) close(1234, status='delete')
#endif
#ifdef CCE_DEBUG
        cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'

        print *,sml_intpl_mype,'Read filename',cce_filename
#else
        cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
#endif
        call adios_read_open_file (buf_id, cce_filename, adios_read_method, MPI_COMM_SELF, err)
        if(err/=0) then
          print *, 'coupling receive error: could not open file', cce_filename
          !call sleep(0.01)
          !stop
        endif
      enddo
      allocate(arrtmp(cce_node_count))
      arrtmp=0D0
      call adios_schedule_read (buf_id, sel1, 'data', 0, 1, arrtmp, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)

      cce_density(cce_first_node:cce_last_node)=arrtmp(1:cce_node_count)

      deallocate(arrtmp)

    endif

    call t_stopf("CCE_RECEIVE_DENSITY")
  end subroutine cce_receive_density

  subroutine cce_process_density(density)
    use sml_module

    real(8), dimension(:), intent(inout) :: density
    real(8) :: alpha

    integer :: ipsi
    integer :: ipsi0
    integer :: ipsi1

    character(5) :: cce_stepstr
    character(5) :: planestr
    character(512) :: cce_filename

    write(cce_stepstr,'(I0.5)') cce_step
    write(planestr,'(I0.5)') sml_intpl_mype

    select case (cce_density_model)
    case (0)
      ! nothing
    case(6)
      ipsi0 = cce_surface_first_node(cce_first_surface) !cce_first_node
      ipsi1 = cce_surface_last_node(cce_last_surface)   !cce_last_node
      density(ipsi0:ipsi1) = density(ipsi0:ipsi1) + cce_density(ipsi0:ipsi1)
    case default
      print *,'Unknown coupling model'
      stop
    end select
    cce_step = cce_step + 1
  end subroutine cce_process_density
#endif

#ifndef GENE_SIDE
  !#ifdef XGC_COUPLING_CORE_EDGE_FIELD
  !Send the density to the other side
  subroutine cce_send_field(dpot0,dpot1,pot0,flag_pot0)
    use sml_module

    include 'mpif.h'

    real*8, dimension(:), intent(in) :: dpot0,dpot1,pot0
    !real*8, dimension(:,:), intent(in) :: dpot
    integer, intent(in) :: flag_pot0

    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename

    ! ADIOS
    integer*8 :: buf_id, buf_size, total_size
    integer :: err,mode

    real*8, dimension(:),allocatable :: arrtmp1,arrtmp2,arrtmp3

    !cce_density=0D0

    if(cce_side.eq.1.and.cce_comm_field_mode.GT.0)then

      write(cce_stepstr,'(I0.5)') cce_field_step
      write(planestr,'(I0.5)') sml_intpl_mype

      allocate(arrtmp1(cce_field_node_number),arrtmp2(cce_field_node_number),arrtmp3(cce_field_node_number))

      arrtmp1(1:cce_field_node_number)=pot0(cce_field_first_node:cce_field_last_node)
      arrtmp2(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
      arrtmp3(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
      !#ifdef SPECIFIC_GENEXGC
      !TODO: write dpot0 + pot0 in arrtmp ============== TODO !!!!
      !arrtmp(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)+pot0(cce_field_first_node:cce_field_last_node)
      !#else
      !arrtmp(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
      !#endif
#ifdef CCE_DEBUG
      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'.bp'
#endif
      ADIOS_OPEN(buf_id,'coupling_fields',cce_filename,'w',MPI_COMM_SELF,err)  !sml_comm_null,err)
      buf_size= 4*8 + 8 + 24*cce_field_node_number  + 100 !last 100 is buffer
      ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
      ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
      ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
      ADIOS_WRITE_LBL(buf_id,'first_node',cce_field_first_node,err)
      ADIOS_WRITE_LBL(buf_id,'last_node',cce_field_last_node,err)
      ADIOS_WRITE_LBL(buf_id,'node_number',cce_field_node_number,err)
      ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
      ADIOS_WRITE_LBL(buf_id,'cce_model',cce_field_model,err)
      ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
      ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
      ADIOS_WRITE_LBL(buf_id,'pot0', arrtmp1,err)
      ADIOS_WRITE_LBL(buf_id,'dpot0', arrtmp2,err)
      ADIOS_WRITE_LBL(buf_id,'dpot1', arrtmp3,err)
      ADIOS_CLOSE(buf_id,err)

      deallocate(arrtmp1)
      deallocate(arrtmp2)
      deallocate(arrtmp3)

#ifdef OLDIMPLEMENTATIOBOBSOLETE
#ifndef SPECIFIC_GENEXGC
      arrtmp(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
#ifdef CCE_DEBUG
      cce_filename=trim(cce_output_dir)//'/dpot1_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
      cce_filename=trim(cce_output_dir)//'/dpot1_'//trim(cce_my_side)//'_'//trim(planestr)//'.bp'
#endif
      ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',MPI_COMM_SELF,err)  !sml_comm_null,err)
      buf_size= 4*8 + 8 + 8*cce_field_node_number  + 100 !last 100 is buffer
      ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
      ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
      ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
      ADIOS_WRITE_LBL(buf_id,'first_node',cce_field_first_node,err)
      ADIOS_WRITE_LBL(buf_id,'last_node',cce_field_last_node,err)
      ADIOS_WRITE_LBL(buf_id,'node_number',cce_field_node_number,err)
      ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
      ADIOS_WRITE_LBL(buf_id,'cce_model',cce_field_model,err)
      ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
      ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
      ADIOS_WRITE_LBL(buf_id,'data', arrtmp,err)
      ADIOS_CLOSE(buf_id,err)

      if(flag_pot0.eq.0)then
        arrtmp(1:cce_field_node_number)=pot0(cce_field_first_node:cce_field_last_node)
#ifdef CCE_DEBUG
        cce_filename=trim(cce_output_dir)//'/pot0_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
        cce_filename=trim(cce_output_dir)//'/pot0_'//trim(cce_my_side)//'_'//trim(planestr)//'.bp'
#endif
        ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',MPI_COMM_SELF,err)  !sml_comm_null,err)
        buf_size= 4*8 + 8 + 8*cce_field_node_number  + 100 !last 100 is buffer
        ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
        ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
        ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
        ADIOS_WRITE_LBL(buf_id,'first_node',cce_field_first_node,err)
        ADIOS_WRITE_LBL(buf_id,'last_node',cce_field_last_node,err)
        ADIOS_WRITE_LBL(buf_id,'node_number',cce_field_node_number,err)
        ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
        ADIOS_WRITE_LBL(buf_id,'cce_model',cce_field_model,err)
        ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
        ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
        ADIOS_WRITE_LBL(buf_id,'data', arrtmp,err)
        ADIOS_CLOSE(buf_id,err)
      endif
      deallocate(arrtmp)
      !Create an unlock file
#endif
#endif

      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)
    endif

  end subroutine cce_send_field

  subroutine cce_receive_field(flag_pot0)
    use adios_read_mod
    use sml_module

    include 'mpif.h'

    integer, intent(in) :: flag_pot0

    logical :: ex
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    integer*8 :: buf_id
    integer :: adios_read_method = ADIOS_READ_METHOD_BP, err
    integer*8 :: sel1=0
    real*8, dimension(:),allocatable :: arrtmp1,arrtmp2,arrtmp3
    integer :: stat

    cce_dpot0=0D0
    cce_dpot1=0D0
    cce_pot0=0D0

    if(cce_side.eq.0.and.cce_comm_field_mode.GT.1)then

      ex=.false.

      write(cce_stepstr,'(I0.5)') cce_field_step
      write(planestr,'(I0.5)') sml_intpl_mype

      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
      do while(.NOT.ex)
        inquire(file=cce_filename,EXIST=ex)
      end do
#ifndef CCE_DEBUG
      open(unit=1234, iostat=stat, file=cce_filename, status='old')
      if (stat == 0) close(1234, status='delete')
#endif
#ifdef CCE_DEBUG
      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
#endif
      call adios_read_open_file (buf_id, cce_filename, adios_read_method, MPI_COMM_SELF, err)
      if(err/=0) then
        print *, 'coupling receive error: could not open file', cce_filename
        stop
      endif
      allocate(arrtmp1(cce_field_node_number),arrtmp2(cce_field_node_number),arrtmp3(cce_field_node_number))
      arrtmp1=0D0
      arrtmp2=0D0
      arrtmp3=0D0
      call adios_schedule_read (buf_id, sel1, 'pot0', 0, 1, arrtmp1, err)
      call adios_schedule_read (buf_id, sel1, 'dpot0', 0, 1, arrtmp2, err)
      call adios_schedule_read (buf_id, sel1, 'dpot1', 0, 1, arrtmp3, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
      cce_pot0(cce_field_first_node:cce_field_last_node)=arrtmp1(1:cce_field_node_number)
      cce_dpot0(cce_field_first_node:cce_field_last_node)=arrtmp2(1:cce_field_node_number)
      cce_dpot1(cce_field_first_node:cce_field_last_node)=arrtmp3(1:cce_field_node_number)
      deallocate(arrtmp1)
      deallocate(arrtmp2)
      deallocate(arrtmp3)

#ifdef OLDIMPLEMENTATIOBOBSOLETE
#ifdef CCE_DEBUG
      cce_filename=trim(cce_output_dir)//'/dpot1_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
      cce_filename=trim(cce_output_dir)//'/dpot1_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
#endif
      call adios_read_open_file (buf_id, cce_filename, adios_read_method, MPI_COMM_SELF, err)
      if(err/=0) then
        print *, 'coupling receive error: could not open file', cce_filename
        stop
      endif
      arrtmp=0D0
      call adios_schedule_read (buf_id, sel1, 'data', 0, 1, arrtmp, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
      cce_dpot1(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)

      if(flag_pot0.eq.0)then
        arrtmp=0D0
#ifdef CCE_DEBUG
        cce_filename=trim(cce_output_dir)//'/pot0_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#else
        cce_filename=trim(cce_output_dir)//'/pot0_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
#endif
        call adios_read_open_file (buf_id, cce_filename, adios_read_method, MPI_COMM_SELF, err)
        if(err/=0) then
          print *, 'coupling receive error: could not open file', cce_filename
          stop
        endif
        call adios_schedule_read (buf_id, sel1, 'data', 0, 1, arrtmp, err)
        call adios_perform_reads (buf_id, err)
        call adios_read_close (buf_id, err)
        cce_pot0(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)
      endif

      deallocate(arrtmp)
#endif
    endif
  end subroutine cce_receive_field


  subroutine cce_process_field(dpot0,dpot1,pot0, flag_pot0)
    use sml_module

    include 'mpif.h'

    real*8, dimension(:)  , intent(inout) :: pot0,dpot0,dpot1
    integer, intent(in) :: flag_pot0

    select case (cce_field_model)
    case (0)
      ! nothing
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

    cce_field_step = cce_field_step + 1

  end subroutine cce_process_field
#endif

#ifdef GENE_SIDE
  subroutine cce_receive_field(iphi,myli0,myli1,myli2,data_block,comm)
    use adios_read_mod
    use mpi

    integer, intent(in) :: iphi,myli0,myli1,myli2
    integer, intent(in) :: comm
    real, dimension(myli1:myli2), intent(out) :: data_block

    logical :: ex
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    integer(8) :: adios_handle
    integer :: adios_read_method = ADIOS_READ_METHOD_BP, adios_err
    integer(8) :: bb_sel
    real, dimension(:),allocatable :: arrtmp

    cce_dpot0=0D0

    if(cce_side.eq.0.and.cce_comm_field_mode.GT.1)then

      ex=.false.

      write(cce_stepstr,'(I0.5)') cce_field_step
      write(planestr,'(I0.5)') iphi
      cce_filename=trim(cce_output_dir)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
      write(*,*)'lookfor',cce_filename
      do while(.NOT.ex)
        inquire(file=cce_filename,EXIST=ex)
      end do
      write(*,*)'found',cce_filename
      cce_filename=trim(cce_output_dir)//'/dpot0_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
      call adios_read_open_file (adios_handle, cce_filename, adios_read_method, comm, adios_err)
      if(adios_err.ne.0) then
        print *, 'coupling receive error: could not open file', cce_filename
        stop
      endif

      call adios_selection_boundingbox(bb_sel, 1, (/int(myli1, kind=8), int(1,kind=8)/), shape(data_block, kind=8))
      call adios_schedule_read(adios_handle, bb_sel, "data", 0 , 1, data_block, adios_err)
      call adios_perform_reads(adios_handle, adios_err)
      call adios_read_close(adios_handle, adios_err)
      call adios_selection_delete(bb_sel)

    endif

  end subroutine cce_receive_field

  subroutine cce_process_field()
    cce_field_step=cce_field_step+1
  end subroutine cce_process_field

#endif

#ifndef GENE_SIDE
  subroutine cce_varpi_grid(rho_ff)
    real (8), dimension(:), intent(inout) :: rho_ff
    integer :: ipsi,ipsi0,ipsi1
    real (8) :: varpi

    call cce_initialize()
    if(cce_density_model.eq.6)then
      !Linear weight
      if(cce_side.EQ.0)then !core
        !1. Zero near axis
        ipsi0=cce_surface_first_node(cce_first_surface)
        ipsi1=cce_surface_last_node(cce_first_surface_coupling_axis)
        rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
        !2. Linear increasing weight
        do ipsi=cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
          varpi=dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
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
          varpi=dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
          ipsi0=cce_surface_first_node(ipsi)
          ipsi1=cce_surface_last_node(ipsi)
          varpi=1D0-varpi
          rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
        enddo
        do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
          varpi=dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
          ipsi0=cce_surface_first_node(ipsi)
          ipsi1=cce_surface_last_node(ipsi)
          rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
        enddo
      endif
      !        endif
    endif
  end subroutine cce_varpi_grid
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
  end function cce_varpi_grid
#endif

end module coupling_core_edge



!<adios-group name="coupling">
!  <var name="nphi" type="integer"/>
!  <var name="iphi" type="integer"/>
!  <var name="first_node" type="integer"/>
!  <var name="last_node" type="integer"/>
!  <var name="nodes_number" type="integer"/>
!  <var name="cce_side" type="integer"/>
!  <var name="cce_density_model" type="integer"/>
!  <var name="time" type="real*8"/>
!  <var name="step" type="integer"/>
!
!!  <global-bounds dimensions="nphi,nodes_number" offsets="iphi,0">
!!     <var name="density" type="real*8" dimensions="1,nodes_number"/>
!!  </global-bounds>
!  <global-bounds dimensions="node_number" offsets="0">
!     <var name="density" type="real*8" dimensions="node_number"/>
!  </global-bounds>
!</adios-group>
!
!<method priority="3" method="MPI" iterations="100" group="coupling"/>
