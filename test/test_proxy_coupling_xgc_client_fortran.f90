module overlap
    use, intrinsic :: ISO_C_BINDING
    implicit none
    public
contains
    function in_overlap(dim, id) bind(C) result(fresult)
        use, intrinsic :: ISO_C_BINDING
        implicit none
        integer(C_INT), intent(in), value :: id
        integer(C_INT), intent(in), value :: dim
        integer(C_INT8_T) :: fresult
        if(id >= 22 .AND. id <= 34) then
            if(dim == 2 .OR. dim == 1 .OR. dim == 0) then
                fresult = 1
                RETURN
            end if
        end if
        fresult = 0
    end function in_overlap
end module overlap

Program main
    use pcms
    use mpi
    use overlap
    use, intrinsic :: ISO_C_BINDING
    implicit none
    !    include "mpif.h"
    integer, parameter :: nplanes = 2

    integer :: num_args, ix, nverts, ierror
    type(SWIGTYPE_p_WdmCplClientHandle) :: client
    type(SWIGTYPE_p_WdmCplReverseClassificationHandle) :: reverse_classification
    type(SWIGTYPE_p_WdmCplFieldHandle), dimension(2) :: fields
    type(SWIGTYPE_p_WdmCplFieldAdapterHandle), dimension(2) :: adapters
    character(len = 80), dimension(:), allocatable :: args
    character(len = :), allocatable :: rc_file
    integer(C_LONG), allocatable, target :: data(:)
    integer(C_LONG), pointer :: data_pointer(:)
    integer :: client_comm, plane_comm
    integer :: world_rank, world_size, plane_rank, plane_size, client_rank, client_size
    integer :: plane
    character(10) :: plane_str
    character(100) :: name

    call MPI_INIT(ierror)
    call pcms_kokkos_initialize_without_args()
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, ierror)
    num_args = command_argument_count()
    if (num_args /= 1) then
        stop 1
    end if
    allocate(args(num_args))
    do ix = 1, num_args
        call get_command_argument(ix, args(ix))
    end do
    rc_file = trim(args(1))

    plane = mod(world_rank, plane)
    call MPI_Comm_split(MPI_COMM_WORLD, plane, world_rank, plane_comm, ierror)
    call MPI_Comm_rank(plane_comm, plane_rank, ierror)
    call MPI_Comm_size(plane_comm, plane_size, ierror)
    if (plane_rank .eq. 0) then
        call MPI_Comm_split(MPI_COMM_WORLD, 0, world_rank, client_comm, ierror)
        call MPI_Comm_rank(client_comm, client_rank, ierror)
        call MPI_Comm_size(client_comm, client_size, ierror)
    else
        call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, client_comm, ierror)
        client_rank = -1
        client_size = -1
    end if

    client = pcms_create_client("proxy_couple", client_comm)
    reverse_classification = pcms_load_reverse_classification(rc_file, MPI_COMM_WORLD)
    nverts = pcms_reverse_classification_count_verts(reverse_classification)
    allocate(data(nverts))
    data_pointer => data

    do ix = 0, nplanes-1
        write(plane_str, '(i0)') ix
        name = 'xgc_gids_plane_' // trim(plane_str)
        if (plane .eq. ix) then
            adapters(ix+1) = pcms_create_xgc_field_adapter(trim(name), MPI_COMM_SELF, c_loc(data_pointer), &
                    size(data_pointer), WDMCPL_LONG_INT, reverse_classification, in_overlap)
        else
            adapters(ix+1) = pcms_create_dummy_field_adapter()
        end if
        if ((plane_rank .eq. 0) .and.  (ix .eq. plane)) then
            ! field participates in the communication
            fields(ix+1) = pcms_add_field(client, trim(name), adapters(ix+1), 1)
        else
            fields(ix+1) = pcms_add_field(client, trim(name), adapters(ix+1), 0)
        end if
    end do
    if(plane_rank .eq. 0) then
        do ix = 1, nverts
            data(ix) = ix
        end do
    endif
    call pcms_begin_send_phase(client)
    call pcms_send_field(fields(plane+1))
    call pcms_end_send_phase(client)
    call pcms_begin_receive_phase(client)
    call pcms_receive_field(fields(plane+1))
    call pcms_end_receive_phase(client)
    do ix = 1, nverts
        if(data(ix) /= ix) then
            print*, "ERROR (1): data[", ix, "] should be ", ix, " not ", data(ix), "."
            STOP 1
        end if
    end do
    if (plane_rank .eq. 0) then
        data = data * 2
    end if
    call pcms_begin_send_phase(client)
    call pcms_send_field(fields(plane+1))
    call pcms_end_send_phase(client)
    call pcms_begin_receive_phase(client)
    call pcms_receive_field(fields(plane+1))
    call pcms_end_receive_phase(client)
    do ix = 1, nverts
        if(data(ix) /= 2 * ix) then
            print*, "ERROR (2): data[", ix, "] should be ", 2 * ix, " not ", data(ix), "."
            STOP 1
        end if
    end do
    do ix = 1, nplanes
        call pcms_destroy_field_adapter(adapters(ix))
    end do
    call pcms_destroy_reverse_classification(reverse_classification)
    call pcms_destroy_client(client)
    if (client_comm .ne. MPI_COMM_NULL) then
        call MPI_Comm_free(client_comm, ierror)
    end if
    call MPI_Comm_free(plane_comm, ierror)
    call pcms_kokkos_finalize()
    call MPI_FINALIZE(ierror)

End Program main
