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
    use wdmcpl
    use mpi
    use flcl_util_kokkos_mod
    use overlap
    use, intrinsic :: ISO_C_BINDING
    implicit none
    !    include "mpif.h"

    integer :: num_args, ix, nverts, rank, size, ierror
    type(SWIGTYPE_p_WdmCplClientHandle) :: client
    type(SWIGTYPE_p_WdmCplReverseClassificationHandle) :: reverse_classification
    type(SWIGTYPE_p_WdmCplFieldHandle) :: field
    type(SWIGTYPE_p_WdmCplFieldAdapterHandle) :: xgc_adapter
    character(len=80), dimension(:), allocatable :: args
    character(len=:), allocatable :: rc_file
    integer(C_LONG), allocatable,target :: data(:)
    integer(C_LONG), pointer :: data_pointer(:)


    call MPI_INIT(ierror)
    call kokkos_initialize_without_args();
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

    num_args = command_argument_count()
    if (num_args /= 1) then
        stop 1
    end if
    allocate(args(num_args))
    do ix = 1, num_args
        call get_command_argument(ix, args(ix))
    end do
    rc_file = trim(args(1))

    client = wdmcpl_create_client("Fortran", MPI_COMM_WORLD)
    reverse_classification = wdmcpl_load_reverse_classification(rc_file, MPI_COMM_WORLD)
    nverts = wdmcpl_reverse_classification_count_verts(reverse_classification)
    allocate(data(nverts))
    data_pointer => data
    do ix = 1, nverts
        data(ix) = ix
    end do

    xgc_adapter = wdmcpl_create_xgc_field_adapter("a1", c_loc(data_pointer),&
            nverts, WDMCPL_LONG_INT, reverse_classification, in_overlap)
    field = wdmcpl_add_field(client, "xgc_gids", xgc_adapter)
    call wdmcpl_send_field(field)
    call wdmcpl_receive_field(field)
    do ix = 1, nverts
        if(data(ix) /= ix) then
            print*, "ERROR (1): data[",ix,"] should be ",ix," not ",data(ix),"."
            STOP 1
        end if
        data(ix) = 2*ix
    end do
    call wdmcpl_send_field_name(client,"xgc_gids")
    call wdmcpl_receive_field_name(client,"xgc_gids")
    do ix = 1, nverts
        if(data(ix) /= 2*ix) then
            print*, "ERROR (2): data[",ix,"] should be ",2*ix," not ",data(ix),"."
            STOP 1
        end if
    end do

    call wdmcpl_destroy_field_adapter(xgc_adapter)
    call wdmcpl_destroy_reverse_classification(reverse_classification)
    call wdmcpl_destroy_client(client)
    call kokkos_finalize()
    call MPI_FINALIZE(ierror)

End Program main
