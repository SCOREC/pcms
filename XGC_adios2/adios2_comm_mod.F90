
module adios2_comm_module
    use sml_module
#ifdef ADIOS2
    use adios2
#endif
    implicit none
#ifdef ADIOS2
    type(adios2_adios) :: adios2obj
    type(adios2_engine), allocatable :: list_engines(:)
    integer :: n_engines

    !! very simple single timer
    character(len=128) :: timer_name
    integer :: timer_comm
    integer :: timer_index
    real(kind=8) :: t_start

    interface adios2_comm_get_type
        module procedure adios2_comm_get_type_dp
        module procedure adios2_comm_get_type_dp_arr1d
        module procedure adios2_comm_get_type_dp_arr2d
        module procedure adios2_comm_get_type_dp_arr3d
        module procedure adios2_comm_get_type_dp_arr4d
        module procedure adios2_comm_get_type_dp_arr5d
        module procedure adios2_comm_get_type_integer4
        module procedure adios2_comm_get_type_integer4_arr1d
        module procedure adios2_comm_get_type_integer4_arr2d
        module procedure adios2_comm_get_type_integer4_arr3d
        module procedure adios2_comm_get_type_integer4_arr4d
        module procedure adios2_comm_get_type_integer4_arr5d
        module procedure adios2_comm_get_type_integer8
        module procedure adios2_comm_get_type_integer8_arr1d
        module procedure adios2_comm_get_type_integer8_arr2d
        module procedure adios2_comm_get_type_integer8_arr3d
        module procedure adios2_comm_get_type_integer8_arr4d
        module procedure adios2_comm_get_type_integer8_arr5d
    end interface

    interface adios2_comm_define_variable   
        module procedure adios2_comm_define_variable_dp
        module procedure adios2_comm_define_variable_dp_arr1d
        module procedure adios2_comm_define_variable_dp_arr2d
        module procedure adios2_comm_define_variable_dp_arr3d
        module procedure adios2_comm_define_variable_dp_arr4d
        module procedure adios2_comm_define_variable_dp_arr5d
        module procedure adios2_comm_define_variable_integer4
        module procedure adios2_comm_define_variable_integer4_arr1d
        module procedure adios2_comm_define_variable_integer4_arr2d
        module procedure adios2_comm_define_variable_integer4_arr3d
        module procedure adios2_comm_define_variable_integer4_arr4d
        module procedure adios2_comm_define_variable_integer4_arr5d
        module procedure adios2_comm_define_variable_integer8
        module procedure adios2_comm_define_variable_integer8_arr1d
        module procedure adios2_comm_define_variable_integer8_arr2d
        module procedure adios2_comm_define_variable_integer8_arr3d
        module procedure adios2_comm_define_variable_integer8_arr4d
        module procedure adios2_comm_define_variable_integer8_arr5d
    end interface

    private :: list_engines
    private :: n_engines
    private :: timer_name, timer_comm, timer_index, t_start
    
contains
    subroutine adios2_comm_init(initfile)
        use input_module
        implicit none
        character(len=*), intent(in) :: initfile
        integer :: ierr

        call adios2_init_config(adios2obj, initfile, sml_comm, .true., ierr)
        allocate(list_engines(16))
        n_engines = 0
    end subroutine adios2_comm_init

    subroutine adios2_comm_finalize()
        implicit none
        integer :: ierr
        integer :: i

        do i = 1, n_engines
            if (sml_mype.eq.0) print *, 'ADIOS2: close output ', trim(list_engines(i)%name)
            call adios2_close(list_engines(i), ierr)
        enddo
        call adios2_finalize(adios2obj, ierr)
    end subroutine adios2_comm_finalize

    subroutine adios2_comm_engine_push(engine)
        implicit none
        type(adios2_engine), intent(in) :: engine
        type(adios2_engine), allocatable :: tmp(:)

        if (n_engines.ge.size(list_engines)) then
            if (sml_mype.eq.0) print *, 'Increasing the size of list_engines to ', size(list_engines)*2
            allocate(tmp(size(list_engines)*2))
            tmp(1:n_engines) = list_engines(1:n_engines)
            deallocate(list_engines)
            call move_alloc(tmp,list_engines)
        endif
        n_engines = n_engines+1
        list_engines(n_engines) = engine
        if (sml_mype.eq.0) print *, 'ADIOS2: push to close on finalizing ', trim(list_engines(n_engines)%name)
    end subroutine

    subroutine adios2_comm_time_start(name, index, comm)
        use mpi
        implicit none
        character(len=*), intent(in) :: name
        integer, intent(in) :: index
        integer, intent(in) :: comm
        integer :: rank, err

        timer_name = name
        timer_index = index
        timer_comm = comm

        call mpi_barrier(comm, err);
        t_start = MPI_Wtime();
    end subroutine

    subroutine adios2_comm_time_end(name)
        use mpi
        implicit none
        character(len=*), intent(in) :: name
        real(kind=8) :: t_end, t_elap, t_sum, buf_sum
        real(kind=8), allocatable :: t_elapall(:)
        integer :: comm_rank, comm_size, i, err
        character(len=128) :: filename

        if (timer_name.eq.trim(name)) then
            call mpi_barrier(timer_comm, err);
            t_end = MPI_Wtime();
            t_elap = t_end - t_start

            call mpi_comm_size(timer_comm,comm_size,err)
            call mpi_comm_rank(timer_comm,comm_rank,err)
            !print *, comm_rank, 't_elap', t_elap
            allocate(t_elapall(comm_size))

            call mpi_gather(t_elap, 1, MPI_REAL8, t_elapall, 1, MPI_REAL8, 0, timer_comm, err)
            !call mpi_gather(buf_size, 1, MPI_INTEGER8, buf_sizeall, 1, MPI_INTEGER8, 0, timer_comm, err)

            if (comm_rank.eq.0) then
                write(filename,'("timing_",a,".",i5.5,".dat")') trim(timer_name), timer_index
                open(unit=99, file=filename, status='replace')
                do i=1,comm_size
                    write (99,'(I8,1X,F15.3)') i, t_elapall(i)
                enddo
                close(99)

                print '(1X, "timing: ", A, " (min,max,avg)=", 3ES15.3, 3ES15.3, 3ES15.3)', &
                    trim(timer_name), minval(t_elapall), maxval(t_elapall), sum(t_elapall)/size(t_elapall)
            endif
            deallocate(t_elapall)
        end if
    end subroutine

    function adios2_comm_get_type_dp(x) result(y)
        implicit none
        real(8), intent(in) :: x
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_dp_arr1d(x) result(y)
        implicit none
        real(8), intent(in) :: x(:)
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_dp_arr2d(x) result(y)
        implicit none
        real(8), intent(in) :: x(:,:)
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_dp_arr3d(x) result(y)
        implicit none
        real(8), intent(in) :: x(:,:,:)
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_dp_arr4d(x) result(y)
        implicit none
        real(8), intent(in) :: x(:,:,:,:)
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_dp_arr5d(x) result(y)
        implicit none
        real(8), intent(in) :: x(:,:,:,:,:)
        integer :: y
        
        y = adios2_type_dp
    end function

    function adios2_comm_get_type_integer4(x) result(y)
        implicit none
        integer(4), intent(in) :: x
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer4_arr1d(x) result(y)
        implicit none
        integer(4), intent(in) :: x(:)
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer4_arr2d(x) result(y)
        implicit none
        integer(4), intent(in) :: x(:,:)
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer4_arr3d(x) result(y)
        implicit none
        integer(4), intent(in) :: x(:,:,:)
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer4_arr4d(x) result(y)
        implicit none
        integer(4), intent(in) :: x(:,:,:,:)
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer4_arr5d(x) result(y)
        implicit none
        integer(4), intent(in) :: x(:,:,:,:,:)
        integer :: y
        
        y = adios2_type_integer4
    end function

    function adios2_comm_get_type_integer8(x) result(y)
        implicit none
        integer(8), intent(in) :: x
        integer :: y
        
        y = adios2_type_integer8
    end function

    function adios2_comm_get_type_integer8_arr1d(x) result(y)
        implicit none
        integer(8), intent(in) :: x(:)
        integer :: y
        
        y = adios2_type_integer8
    end function

    function adios2_comm_get_type_integer8_arr2d(x) result(y)
        implicit none
        integer(8), intent(in) :: x(:,:)
        integer :: y
        
        y = adios2_type_integer8
    end function

    function adios2_comm_get_type_integer8_arr3d(x) result(y)
        implicit none
        integer(8), intent(in) :: x(:,:,:)
        integer :: y
        
        y = adios2_type_integer8
    end function

    function adios2_comm_get_type_integer8_arr4d(x) result(y)
        implicit none
        integer(8), intent(in) :: x(:,:,:,:)
        integer :: y
        
        y = adios2_type_integer8
    end function

    function adios2_comm_get_type_integer8_arr5d(x) result(y)
        implicit none
        integer(8), intent(in) :: x(:,:,:,:,:)
        integer :: y
        
        y = adios2_type_integer8
    end function

    !! adios2_comm_define_variable
    subroutine adios2_comm_define_variable_dp(variable, io, name, x, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
    end subroutine

    subroutine adios2_comm_define_variable_dp_arr1d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x(:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_dp_arr2d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x(:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_dp_arr3d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x(:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_dp_arr4d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x(:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_dp_arr5d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(8), intent(in) :: x(:,:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4(variable, io, name, x, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4_arr1d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x(:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4_arr2d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x(:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4_arr3d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x(:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4_arr4d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x(:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer4_arr5d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(4), intent(in) :: x(:,:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine
    
    subroutine adios2_comm_define_variable_integer8(variable, io, name, x, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer8_arr1d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x(:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer8_arr2d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x(:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer8_arr3d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x(:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer8_arr4d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x(:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine

    subroutine adios2_comm_define_variable_integer8_arr5d(&
        variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(8), intent(in) :: x(:,:,:,:,:)
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
            size(shape_dims), shape_dims, start_dims, count_dims, &
            is_constant_dims, ierr)
    end subroutine
#endif
end module adios2_comm_module
