module coupling_base_setup
#ifdef ADIOS
use adios_read_mod
#endif
use mpi

implicit none


! 0: core
! 1: edge
integer :: cce_side=0

integer, dimension(:), allocatable :: cce_surface_first_node
integer, dimension(:), allocatable :: cce_surface_last_node

integer :: cce_density_model
integer :: cce_first_surface
integer :: cce_last_surface
integer :: cce_surface_number
integer :: cce_all_surface_number
integer :: cce_first_node,cce_last_node
integer :: cce_node_number
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
integer :: cce_npsi

real(8) :: cce_alpha
real(8) :: cce_dt
real(8), dimension(:), allocatable :: cce_varpi
real(8), dimension(:), allocatable :: cce_psi
real, dimension(:), allocatable :: cce_density
real, dimension(:), allocatable :: cce_pot0
real, dimension(:), allocatable :: cce_dpot0
real, dimension(:), allocatable :: cce_dpot1
character(256) :: cce_folder
character(256) :: cce_filename
character(256) :: cce_lock_filename

character(5) :: cce_my_side
character(5) :: cce_other_side

logical :: cce_bcast_dpot
logical :: cce_bcast_pot0
logical :: cce_dpot_index0


namelist /surfaces/ cce_surface_first_node, cce_surface_last_node
namelist /varpi/ cce_varpi, cce_psi

namelist /coupling/&
  cce_all_surface_number, &
  cce_alpha, &
  cce_comm_density_mode, &
  cce_comm_field_mode, &
  cce_density_model, &
  cce_field_model, &
  cce_first_node, &
  cce_first_surface, &
  cce_first_surface_coupling, &
  cce_first_surface_coupling_axis, &
  cce_first_surface_field, &
  cce_folder, &
  cce_last_node, &
  cce_last_surface, &
  cce_last_surface_coupling, &
  cce_last_surface_coupling_axis, &
  cce_last_surface_field, &
  cce_npsi, &
  cce_side, &
  cce_dpot_index0, &
  cce_dt

contains


  subroutine set_default_coupler
     cce_alpha = 0.5D0
     cce_density_model = 0
     cce_step = 0 !In case of restart, one might want to change this
     cce_field_step = 0
     cce_first_node = 10
     cce_last_node = 0
     cce_comm_density_mode = 2
     cce_all_surface_number = 1
     cce_first_surface_coupling = -1
     cce_last_surface_coupling = -1
     cce_field_model = 0
     cce_comm_field_mode = 0
     cce_npsi = -1
     cce_dpot_index0=.false.
     cce_dt = -1
  end subroutine set_default_coupler


  subroutine cce_destroy()
    if (allocated(cce_surface_first_node)) deallocate(cce_surface_first_node)
    if (allocated(cce_surface_last_node)) deallocate(cce_surface_last_node)
    if (allocated(cce_density)) deallocate(cce_density)
    if (allocated(cce_dpot0)) deallocate(cce_dpot0)
    if (allocated(cce_dpot1)) deallocate(cce_dpot1)
    if (allocated(cce_pot0)) deallocate(cce_pot0)
  end subroutine cce_destroy

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

end module coupling_base_setup
