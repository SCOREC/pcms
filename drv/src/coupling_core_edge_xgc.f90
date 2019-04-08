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
  character(8) :: cce_my_side
  character(8) :: cce_other_side
  character(8) :: cce_density_step_string
  character(8) :: cce_field_step_string
  character(8) :: cce_plane_string

  logical :: cce_bcast_dpot
  logical :: cce_bcast_pot0

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

  integer :: cce_density_step
  integer :: cce_field_step
  integer :: cce_field_model

  integer :: cce_comm_density_mode
  integer :: cce_comm_field_mode

  ! potential field global vals
  integer :: cce_field_first_node
  integer :: cce_field_last_node
  integer :: cce_field_node_count

  integer :: cce_field_first_surface
  integer :: cce_field_last_surface
  integer :: cce_field_surface_count

  integer :: cce_first_surface_coupling
  integer :: cce_last_surface_coupling

  integer :: cce_first_surface_coupling_axis
  integer :: cce_last_surface_coupling_axis

  ! physical data (fields)
  real, dimension(:,:), allocatable :: cce_density
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
    cce_density_first_node, &
    cce_density_last_node, &
    cce_comm_density_mode, &
    cce_all_surface_number, &
    cce_alpha, &
    cce_field_model, &
    cce_comm_field_mode, &
    cce_first_surface_coupling, &
    cce_last_surface_coupling, &
    cce_field_first_surface, &
    cce_field_last_surface, &
    cce_npsi, &
    cce_first_surface_coupling_axis, &
    cce_last_surface_coupling_axis

  namelist /surfaces/ cce_surface_first_node, cce_surface_last_node
  namelist /varpi/ cce_varpi, cce_psi

contains

  subroutine cce_set_defaults
    implicit none
    cce_input_file = 'coupling.in'
    ccs_adios_group = 'coupling'
    cce_alpha = 0.5D0
    cce_density_model = 0
    cce_density_step = 1 !In case of restart, one might want to change this
    cce_field_step = 1
    cce_density_first_node = 10
    cce_density_last_node = 0
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

      cce_surface_count = cce_last_surface - cce_first_surface + 1
      cce_density_node_count = cce_density_last_node - cce_density_first_node + 1

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
      allocate(cce_density(cce_density_first_node:cce_density_last_node,cce_first_surface:cce_last_surface))

      if (cce_field_model /= 0 .and. cce_comm_field_mode /= 0) then
        allocate(cce_dpot0(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
        allocate(cce_dpot1(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
        allocate(cce_pot0(cce_field_first_node:cce_field_last_node,cce_field_first_surface:cce_field_last_surface))
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

  subroutine cce_advance_density
    cce_density_step = cce_density_step + 1
    write(cce_density_step_string,'(I0.5)') cce_density_step
  end subroutine cce_advance_density

  subroutine cce_advance_field
    cce_field_step = cce_field_step + 1
    write(cce_field_step_string,'(I0.5)') cce_field_step
  end subroutine cce_advance_field

  subroutine cce_preprocess_density(density, plane_id)
    real(8), dimension(:), intent(inout) :: density
    integer, intent(in) :: plane_id
    real(8) :: alpha

    integer :: ipsi
    integer :: ipsi0
    integer :: ipsi1

    character(512) :: cce_filename

    select case (cce_density_model)
    case (-1)
      density(cce_density_first_node:cce_density_last_node) = 0D0
    case (0)
      ! nothing
    case(1) !linear coupling
      if ((cce_side == 1) .AND. (cce_first_surface < cce_first_surface_coupling)) then
        ipsi0 = cce_surface_first_node(cce_first_surface)
        ipsi1 = cce_surface_last_node(cce_first_surface_coupling)
        density(ipsi0:ipsi1) = cce_density(ipsi0:ipsi1)
      endif
      do ipsi = cce_first_surface_coupling, cce_last_surface_coupling
        alpha = dble(ipsi-cce_first_surface_coupling) / dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
        ipsi0 = cce_surface_first_node(ipsi)
        ipsi1 = cce_surface_last_node(ipsi)
        if(cce_side == 1) alpha = 1D0 - alpha
        density(ipsi0:ipsi1) = (1D0 - alpha) * density(ipsi0:ipsi1) + alpha * cce_density(ipsi0:ipsi1)
      enddo
      if ((cce_side == 0).AND.(cce_last_surface > cce_last_surface_coupling)) then
        ipsi0 = cce_surface_first_node(cce_last_surface_coupling)
        ipsi1 = cce_surface_last_node(cce_last_surface)
        density(ipsi0:ipsi1) = cce_density(ipsi0:ipsi1)
      endif
    case(2) ! average
      alpha = cce_alpha
      ipsi0 = cce_field_first_node
      ipsi1 = cce_density_last_node
      density(ipsi0:ipsi1) = (1D0 - alpha) * density(ipsi0:ipsi1) + alpha * cce_density(ipsi0:ipsi1)
    case(3) ! identity (doesn't do anything...)
      ipsi0 = cce_density_first_node
      ipsi1 = cce_density_last_node
      density(ipsi0:ipsi1) = density(ipsi0:ipsi1)
    case(4)
      alpha = cce_alpha
      ipsi0 = cce_density_first_node
      ipsi1 = cce_density_last_node
      density(ipsi0:ipsi1) = (1D0 - alpha) * cce_density(ipsi0:ipsi1)
    case(5)
      if ((cce_side == 1)) then
        ipsi0 = cce_surface_first_node(cce_first_surface)
        ipsi1 = cce_surface_last_node(cce_last_surface_coupling)
        density(ipsi0:ipsi1) = cce_density(ipsi0:ipsi1)
      else if ((cce_side == 0)) then
        ipsi0 = cce_surface_first_node(cce_first_surface_coupling)
        ipsi1 = cce_surface_last_node(cce_last_surface_coupling)
        density(ipsi0:ipsi1) = cce_density(ipsi0:ipsi1)
      endif
    case(6)
      ipsi0 = cce_surface_first_node(cce_first_surface) !cce_density_first_node
      ipsi1 = cce_surface_last_node(cce_last_surface)   !cce_density_last_node
      density(ipsi0:ipsi1) = density(ipsi0:ipsi1) + cce_density(ipsi0:ipsi1)
    case default
      print *,'Unknown coupling model'
      stop
    end select
  end subroutine cce_preprocess_density

  subroutine cce_preprocess_field(dpot0, dpot1, pot0, flag_pot0)
    implicit none
    real(8), dimension(:), intent(inout) :: pot0
    real(8), dimension(:), intent(inout) :: dpot0
    real(8), dimension(:), intent(inout) :: dpot1
    integer, intent(in) :: flag_pot0

    select case (cce_field_model)
    case (0)
      ! nothing
    case(1)
      cce_bcast_dpot = .true.
      dpot0(:) = cce_dpot0(:)
      dpot1(:) = cce_dpot1(:)
      if (flag_pot0 == 0) then
        cce_bcast_pot0 = .true.
        pot0(:) = cce_pot0(:)
      endif
    case default
      print *,'Unknown coupling model'
      stop
    end select
    cce_field_step = cce_field_step + 1
  end subroutine cce_preprocess_field

  ! all values local, though this assumes that the
  !  indices used in density(:,:) are global indices instead of local indices...?
  subroutine cce_send_density(density, nd_bgn, nd_cnt, pln_bgn, pln_cnt, comm)
    implicit none

    real(8), dimension(:,:), intent(in) :: density
    integer, intent(in) :: nd_bg
    integer, intent(in) :: nd_cnt
    integer, intent(in) :: pln_bgn
    integer, intent(in) :: pln_cnt
    integer, intent(in) :: comm

    integer :: nd_end
    integer :: pln_end
    character(:), allocatable :: cce_filename

    ! adios
    integer(8) :: buf_id
    integer(8) :: buf_grp_sz
    integer(8) :: buf_sz
    integer :: err

    nd_end = nd_bgn + nd_cnt
    pln_end = pln_bgn + pln_cnt

    if (cce_comm_density_mode == 1 .or. cce_comm_density_mode == 2) then
      write(cce_plane_string,'(I0.5)') pln_bgn
      cce_filename = trim(cce_output_dir)//'/density_'//trim(cce_my_side)//'_'//trim(cce_plane_string)//'_'//trim(cce_density_step_string)//'.bp'
      cce_density(:,:) = density(nd_bgn:nd_end,pln_bgn:pln_end)
      buf_grp_sz = 24 + 8 * (nd_cnt * pln_cnt)
      call adios_open(buf_id, cce_adios_group, cce_filename, 'w', comm, err)
      call adios_group_size (buf_id, buf_grp_sz, buf_sz, err)
      call adios_write (buf_id, "gbl_nd_cnt", cce_density_node_count, err)
      call adios_write (buf_id, "lcl_nd_cnt", nd_cnt, err)
      call adios_write (buf_id, "lcl_nd_bgn", nd_bgn, err)
      call adios_write (buf_id, "gbl_pln_cnt", cce_surface_count, err)
      call adios_write (buf_id, "lcl_pln_cnt", pln_cnt, err)
      call adios_write (buf_id, "lcl_pln_bgn", pln_bgn, err)
      call adios_write (buf_id, "density", cce_density, err)
      call adios_close (buf_id, err)
    endif
  end subroutine cce_send_density

  subroutine cce_send_field(pot0, dpot0, dpot1, nd_bgn, nd_cnt, pln_bgn, pln_cnt, comm)
    implicit none

    real(8), dimension(:), intent(in) :: pot0
    real(8), dimension(:), intent(in) :: dpot0
    real(8), dimension(:), intent(in) :: dpot1
    integer, intent(in) :: nd_bgn
    integer, intent(in) :: nd_cnt
    integer, intent(in) :: pln_bgn
    integer, intent(in) :: pln_cnt
    integer, intent(in) :: comm

    integer :: nd_end
    integer :: pln_end
    character(512) :: cce_filename

    integer(8) :: buf_id
    integer(8) :: buf_grp_sz
    integer(8) :: buf_sz
    integer :: err

    nd_end = nd_bgn + nd_cnt
    pln_end = pln_bgn + pln_cnt

    if (cce_side == 1 .and. cce_comm_field_mode > 0) then
      cce_filename = trim(cce_output_dir)//'/potential_'//trim(cce_my_side)//'_'//trim(cce_density_step_string)//'.bp'

      cce_pot0(:,:)  = pot0(nd_bgn:nd_end, pln_bgn:pln_end)
      cce_dpot0(:,:) = dpot0(nd_bgn:nd_end, pln_bgn:pln_end)
      cce_dpot1(:,:) = dpot1(nd_bgn:nd_end, pln_bgn:pln_end)

      call adios_open(buf_id,'coupling_fields',cce_filename,'w',MPI_COMM_SELF,err)
      buf_grp_sz = (24 + 24 * cce_field_node_count)
      call adios_group_size(buf_id, buf_sz, buf_sz, err)
      call adios_write(buf_id, 'gbl_nd_cnt', cce_field_node_count, err)
      call adios_write(buf_id, 'lcl_nd_cnt', nd_cnt, err)
      call adios_write(buf_id, 'lcl_nd_bgn', nd_bgn, err)
      call adios_write(buf_id, 'gbl_pln_cnt', cce_field_surface_count, err)
      call adios_write(buf_id, 'lcl_pln_cnt', pln_cnt, err)
      call adios_write(buf_id, 'lcl_pln_bgn', pln_bgn, err)
      call adios_write(buf_id, 'pot0', cce_pot0, err)
      call adios_write(buf_id, 'dpot0', cce_dpot0, err)
      call adios_write(buf_id, 'dpot1', cce_dpot1, err)
      call adios_close(buf_id, err)
    endif
  end subroutine cce_send_field

  subroutine cce_receive_density(density, nd_bgn, nd_cnt, pln_bgn, pln_cnt, comm)
    implicit none

    real, dimension(:,:), intent(inout) :: density
    integer, intent(in) :: nd_bgn
    integer, intent(in) :: nd_cnt
    integer, intent(in) :: pln_id
    integer, intent(in) :: pln_cnt
    integer, intent(in) :: comm

    character(512) :: cce_filename
    integer(8) :: buf_id
    integer(8) :: err
    integer(8) :: bb_sel
    integer :: stat
    integer(8), dimension(2) :: bnds
    integer(8), dimension(2) :: cnts

    if (cce_comm_density_mode == 2.or.cce_comm_density_mode == 3) then
      write(cce_plane_string,'(I0.5)') pln_id
      cce_filename=trim(cce_output_dir)//'/density_'//trim(cce_other_side)//'_'//trim(cce_plane_string)//'_'//trim(cce_density_step_string)//'.bp'

      bnds(1) = int(nd_bgn, kind=8)
      bnds(2) = int(pln_bgn, kind=8)
      cnts(1) = int(nd_cnt, kind=8)
      cnts(2) = int(pln_cnt, kind=8)

      call adios_selection_boundingbox (bb_sel, 1, bnds, cnts)
      call adios_read_open_file (buf_id, cce_filename, ADIOS_READ_METHOD_BP, comm, err)
      call adios_schedule_read (buf_id, bb_sel, 'density', 0, 1, cce_density, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
    endif
  end subroutine cce_receive_density

  subroutine cce_receive_field(pot0, dpot0, dpot1, flag_pot0, nd_bgn, nd_cnt, pln_bgn, pln_cnt, comm)
    implicit none
    integer, intent(in) :: flag_pot0
    real(8), dimension(:), intent(in)

    character(512) :: cce_filename

    ! adios
    integer(8) :: buf_id
    integer :: err
    integer :: bb_sel
    integer(8), dimension(2) :: bnds
    integer(8), dimension(2) :: cnts

    if (cce_side == 0 .and. cce_comm_field_mode > 1) then

      bnds(1) = nd_bgn
      cnts(1) = nd_cnt
      bnds(2) = pln_bgn
      cnts(2) = pln_cnt

      cce_filename = trim(cce_output_dir)//'/field_'//trim(cce_other_side)//'.bp'

      call adios_selection_boundingbox(bb_sel, 1, bnds, cnts)
      call adios_read_open_file (buf_id, cce_filename, ADIOS_READ_METHOD_BP, comm, err)
      call adios_schedule_read (buf_id, bb_sel, 'pot0', cce_field_step, 1, cce_pot0, err)
      call adios_schedule_read (buf_id, bb_sel, 'dpot0', cce_field_step, 1, cce_dpot0, err)
      call adios_schedule_read (buf_id, bb_sel, 'dpot1', cce_field_step, 1, cce_dpot1, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
    endif
  end subroutine cce_receive_field

  ! i *think* this is taking the (gyro averaged) ion density for (in XGC) a single
  !  poloidal plane (since each process has only 2), and for a single flux surface,
  !  but for all nodes on the flux surface, and projecting those values onto the coupling grid
  !  but is seems to do this for all flux surfaces, rather than some subset of the flux surfaces
  !  which I might expect for a general coupling which has a limited overlap region

  ! oh wait, this might by APPLYING the coupling, since it seems to ZERO out those
  ! nodes corresponding to the overlap region (in some cases) and then
  ! add values to those locations, but I don't see the density field in here

  ! it could also just be setting some weighting values for use by the actual gyro-averagin function...
  subroutine cce_varpi_grid(rho_ff)
    ! XGC VERSION
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
