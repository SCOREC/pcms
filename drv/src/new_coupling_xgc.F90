module new_coupling_xgc
  use adios_read_mod
!  use adios_write_mod
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
  integer :: cce_node_number                    !added for variable signature consistency
  integer :: cce_first_node                     !added for variable signature consistency
  integer :: cce_last_node                      !added for variable signature consistency
  integer :: staging_read_method = ADIOS_READ_METHOD_BP

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
subroutine cce_initialize()
  use sim_param                 ! substitute routine to sml_module
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



subroutine cce_destroy

  if(allocated(cce_surface_first_node)) deallocate(cce_surface_first_node)
  if(allocated(cce_surface_last_node)) deallocate(cce_surface_last_node)
  if(allocated(cce_density)) deallocate(cce_density)
  if(allocated(cce_dpot0)) deallocate(cce_dpot0)
  if(allocated(cce_dpot1)) deallocate(cce_dpot1)
  if(allocated(cce_pot0)) deallocate(cce_pot0)
 
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
  !use adios_write_mod   
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

end subroutine 


subroutine cce_send_field(pot0, dpot0, dpot1, flag_pot0)
  use sim_param
  implicit none

  real(8), dimension(:), intent(in) :: pot0
  real(8), dimension(:), intent(in) :: dpot0
  real(8), dimension(:), intent(in) :: dpot1
  integer, intent(in) :: flag_pot0

character(5)::cce_stepstr,planestr
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

    allocate(arrtmp(cce_field_node_number))

#ifdef SPECIFIC_GENEXGC
    if(cce_dpot_index0)then
      arrtmp(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
      if(cce_nzeromode) arrtmp(1:cce_field_node_number)=arrtmp(1:cce_field_node_number)+pot0(cce_field_first_node:cce_field_last_node)
#else
    arrtmp(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
#endif
    if (staging_read_method.gt.0) then
      cce_filename='dpot0_'//trim(cce_my_side)//'.bp'
    else
      cce_filename=trim(cce_folder)//'/dpot0_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
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

#ifdef SPECIFIC_GENEXGC
    else  !else of "if(cce_dpot_index0)then"
      arrtmp(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
      if(cce_nzeromode) arrtmp(1:cce_field_node_number)=arrtmp(1:cce_field_node_number)+pot0(cce_field_first_node:cce_field_last_node)
#else
    arrtmp(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
#endif
    if (staging_read_method.gt.0) then
      cce_filename='dpot1_'//trim(cce_my_side)//'.bp'
    else
      cce_filename=trim(cce_folder)//'/dpot1_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
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
#ifdef SPECIFIC_GENEXGC
    endif !endif of "if(cce_dpot_index0)then"
#else
    if(flag_pot0.eq.0)then
      arrtmp(1:cce_field_node_number)=pot0(cce_field_first_node:cce_field_last_node)
      if (staging_read_method.gt.0) then
        cce_filename='pot0_'//trim(cce_my_side)//'.bp'
      else
        cce_filename=trim(cce_folder)//'/pot0_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
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
    endif
    deallocate(arrtmp)
    !Create an unlock file
#endif

    !Create an unlock file: one file per step
    call mpi_barrier(sml_intpl_comm, err)
    if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
      cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)
    endif

  endif

!  call t_stopf("CCE_SEND_FIELD")
end subroutine



subroutine cce_receive_density()
  use adios_read_mod
  use sim_param                 !  use sml_module
  implicit none

  logical :: ex
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
  integer(8) :: buf_id, buf_size, total_size
  integer :: err
  integer(8) :: sel2, start2(2), count2(2)
  real(8), dimension(:),allocatable :: arrtmp
!  call t_startf("CCE_RECEIVE_DENSITY")

  cce_density=0D0

  !if(cce_side.eq.1.and.cce_comm_density_mode.GT.1)then
  if(cce_comm_density_mode.eq.2.or.cce_comm_density_mode.eq.3)then
    ex=.false.

    write(cce_stepstr,'(I0.5)') cce_step
    write(planestr,'(I0.5)') sml_intpl_mype

#ifndef CCE_TEST_1

    cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.unlock'
!    print *,'Wait unlock filename',cce_filename
    do while((.not.ex).and.(staging_read_method.eq.0))
      inquire(file=cce_filename,EXIST=ex)
    end do

    cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
    if (staging_read_method.gt.0) then
      cce_filename='density_'//trim(cce_other_side)//'.bp'
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
      !stop
    endif

    allocate(arrtmp(cce_node_number))
    arrtmp=0D0
    start2(1) = 0
    count2(1) = cce_node_number
    start2(2) = sml_intpl_mype
    count2(2) = 1
    call adios_selection_boundingbox(sel2, 2, start2, count2)
    call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)

#else

    allocate(arrtmp(cce_node_number))
    arrtmp=0D0
    if(sml_intpl_mype.eq.0) arrtmp(90000)=1D0

    cce_filename=trim(cce_folder)//'/density_test_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'

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

#endif

    cce_density(cce_first_node:cce_last_node)=arrtmp(1:cce_node_number)

    !print *,trim(cce_my_side)//'density_'//trim(planestr)//'_'//trim(cce_stepstr)//'_R',arrtmp(1),arrtmp(cce_node_number)
    deallocate(arrtmp)

  endif
!  call t_stopf("CCE_RECEIVE_DENSITY")

end subroutine cce_receive_density        

subroutine cce_process_density(density)
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
