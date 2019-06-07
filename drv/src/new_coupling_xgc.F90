module new_coupling_xgc
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
    cce_folder, &
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
    cce_last_surface_coupling_axis, &
    cce_nzeromode

  namelist /surfaces/ cce_surface_first_node, cce_surface_last_node
  namelist /varpi/ cce_varpi, cce_psi

contains

subroutine cce_set_defaults
    implicit none
    cce_input_file = 'coupling.in'
    cce_adios_group = 'coupling'
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
    cce_nzeromode = .false.
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



subroutine cce_destroy

end subroutine cce_destroy


subroutine cce_advance_density

end subroutine cce_advance_density


subroutine cce_advance_field

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

subroutine cce_send_density(density)
implicit none
    real(8), dimension(:), intent(in) :: density
end subroutine 



subroutine cce_send_field(pot0, dpot0, dpot1, comm)
    implicit none

    real(8), dimension(:), intent(in) :: pot0
    real(8), dimension(:), intent(in) :: dpot0
    real(8), dimension(:), intent(in) :: dpot1
    integer, intent(in) :: comm
 
end subroutine



subroutine cce_receive_density()

end subroutine cce_receive_density        

subroutine cce_process_density(density)
        implicit none
        real(8), dimension(:), intent(inout) :: density

end subroutine cce_process_density

subroutine cce_receive_field(flag_pot0)
  implicit none
  integer, intent(in) :: flag_pot0

end subroutine cce_receive_field

subroutine cce_process_field(dpot0, dpot1, pot0, flag_pot0)
    implicit none

    integer, intent(in) :: flag_pot0
    real(8), dimension(:), intent(inout) :: pot0
    real(8), dimension(:), intent(inout) :: dpot0
    real(8), dimension(:), intent(inout) :: dpot1

end subroutine cce_process_field


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

end module new_coupling_xgc
