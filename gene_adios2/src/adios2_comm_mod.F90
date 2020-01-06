module adios2_comm_module
#ifdef ADIOS2
  use discretization, only: mype
  use adios2
#endif
  implicit none
#ifdef ADIOS2

  public :: adios_cfg
  public :: adios2_comm_init, adios2_comm_finalize
  public :: engines, adios2obj

  private

  interface adios2_comm_get_type
     module procedure adios2_comm_get_type_dp
     module procedure adios2_comm_get_type_dp_arr1d
     module procedure adios2_comm_get_type_dp_arr2d
     module procedure adios2_comm_get_type_dp_arr3d
     module procedure adios2_comm_get_type_dp_arr4d
     module procedure adios2_comm_get_type_dp_arr5d
     module procedure adios2_comm_get_type_cdp
     module procedure adios2_comm_get_type_cdp_arr1d
     module procedure adios2_comm_get_type_cdp_arr2d
     module procedure adios2_comm_get_type_cdp_arr3d
     module procedure adios2_comm_get_type_cdp_arr4d
     module procedure adios2_comm_get_type_cdp_arr5d
     module procedure adios2_comm_get_type_c
     module procedure adios2_comm_get_type_c_arr1d
     module procedure adios2_comm_get_type_c_arr2d
     module procedure adios2_comm_get_type_c_arr3d
     module procedure adios2_comm_get_type_c_arr4d
     module procedure adios2_comm_get_type_c_arr5d
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
     module procedure adios2_comm_define_variable_cdp
     module procedure adios2_comm_define_variable_cdp_arr1d
     module procedure adios2_comm_define_variable_cdp_arr2d
     module procedure adios2_comm_define_variable_cdp_arr3d
     module procedure adios2_comm_define_variable_cdp_arr4d
     module procedure adios2_comm_define_variable_cdp_arr5d
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

  type(adios2_adios) :: adios2obj
  type(adios2_engine), allocatable :: engines(:)
  integer :: n_engines
  character(len=128):: adios_cfg='adioscfg.xml'
contains
  subroutine adios2_comm_init(initfile, comm)
    character(len=*), intent(in) :: initfile
    integer, intent(in):: comm
    integer :: ierr

    call adios2_init_config(adios2obj, initfile, comm, .true., ierr)
    
    n_engines=2
    !first write second read
    allocate(engines(n_engines))

  end subroutine adios2_comm_init

  subroutine adios2_comm_finalize()
    implicit none
    integer :: ierr
    integer :: i

    do i = 1, n_engines
!       if (mype.eq.0) print *, 'ADIOS2: close output ', trim(list_engines(i)%name)
       call adios2_close(engines(i), ierr)
    enddo
    call adios2_finalize(adios2obj, ierr)

  end subroutine adios2_comm_finalize

  function adios2_comm_get_type_dp(x) result(y)
    implicit none
    real(8), intent(in) :: x
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp

  function adios2_comm_get_type_dp_arr1d(x) result(y)
    implicit none
    real(8), intent(in) :: x(:)
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp_arr1d

  function adios2_comm_get_type_dp_arr2d(x) result(y)
    implicit none
    real(8), intent(in) :: x(:,:)
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp_arr2d

  function adios2_comm_get_type_dp_arr3d(x) result(y)
    implicit none
    real(8), intent(in) :: x(:,:,:)
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp_arr3d

  function adios2_comm_get_type_dp_arr4d(x) result(y)
    implicit none
    real(8), intent(in) :: x(:,:,:,:)
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp_arr4d

  function adios2_comm_get_type_dp_arr5d(x) result(y)
    implicit none
    real(8), intent(in) :: x(:,:,:,:,:)
    integer :: y

    y = adios2_type_dp
  end function adios2_comm_get_type_dp_arr5d

  function adios2_comm_get_type_cdp(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp

  function adios2_comm_get_type_cdp_arr1d(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x(:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp_arr1d

  function adios2_comm_get_type_cdp_arr2d(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x(:,:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp_arr2d

  function adios2_comm_get_type_cdp_arr3d(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x(:,:,:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp_arr3d

  function adios2_comm_get_type_cdp_arr4d(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x(:,:,:,:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp_arr4d

  function adios2_comm_get_type_cdp_arr5d(x) result(y)
    implicit none
    complex(kind=8), intent(in) :: x(:,:,:,:,:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_cdp_arr5d

  function adios2_comm_get_type_c(x) result(y)
    implicit none
    complex(4), intent(in) :: x
    integer :: y

    y = adios2_type_complex
  end function adios2_comm_get_type_c

  function adios2_comm_get_type_c_arr1d(x) result(y)
    implicit none
    complex(4), intent(in) :: x(:)
    integer :: y

    y = adios2_type_complex
  end function adios2_comm_get_type_c_arr1d

  function adios2_comm_get_type_c_arr3d(x) result(y)
    implicit none
    complex(4), intent(in) :: x(:,:,:)
    integer :: y

    y = adios2_type_complex
  end function adios2_comm_get_type_c_arr3d

  function adios2_comm_get_type_c_arr2d(x) result(y)
    implicit none
    complex(4), intent(in) :: x(:,:)
    integer :: y

    y = adios2_type_complex_dp
  end function adios2_comm_get_type_c_arr2d

  function adios2_comm_get_type_c_arr4d(x) result(y)
    implicit none
    complex(4), intent(in) :: x(:,:,:,:)
    integer :: y

    y = adios2_type_complex
  end function adios2_comm_get_type_c_arr4d

  function adios2_comm_get_type_c_arr5d(x) result(y)
    implicit none
    complex(4), intent(in) :: x(:,:,:,:,:)
    integer :: y

    y = adios2_type_complex
  end function adios2_comm_get_type_c_arr5d

  function adios2_comm_get_type_integer4(x) result(y)
    implicit none
    integer(4), intent(in) :: x
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4

  function adios2_comm_get_type_integer4_arr1d(x) result(y)
    implicit none
    integer(4), intent(in) :: x(:)
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4_arr1d

  function adios2_comm_get_type_integer4_arr2d(x) result(y)
    implicit none
    integer(4), intent(in) :: x(:,:)
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4_arr2d

  function adios2_comm_get_type_integer4_arr3d(x) result(y)
    implicit none
    integer(4), intent(in) :: x(:,:,:)
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4_arr3d

  function adios2_comm_get_type_integer4_arr4d(x) result(y)
    implicit none
    integer(4), intent(in) :: x(:,:,:,:)
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4_arr4d

  function adios2_comm_get_type_integer4_arr5d(x) result(y)
    implicit none
    integer(4), intent(in) :: x(:,:,:,:,:)
    integer :: y

    y = adios2_type_integer4
  end function adios2_comm_get_type_integer4_arr5d

  function adios2_comm_get_type_integer8(x) result(y)
    implicit none
    integer(8), intent(in) :: x
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8

  function adios2_comm_get_type_integer8_arr1d(x) result(y)
    implicit none
    integer(8), intent(in) :: x(:)
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8_arr1d

  function adios2_comm_get_type_integer8_arr2d(x) result(y)
    implicit none
    integer(8), intent(in) :: x(:,:)
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8_arr2d

  function adios2_comm_get_type_integer8_arr3d(x) result(y)
    implicit none
    integer(8), intent(in) :: x(:,:,:)
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8_arr3d

  function adios2_comm_get_type_integer8_arr4d(x) result(y)
    implicit none
    integer(8), intent(in) :: x(:,:,:,:)
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8_arr4d

  function adios2_comm_get_type_integer8_arr5d(x) result(y)
    implicit none
    integer(8), intent(in) :: x(:,:,:,:,:)
    integer :: y

    y = adios2_type_integer8
  end function adios2_comm_get_type_integer8_arr5d

  !! adios2_comm_define_variable
  subroutine adios2_comm_define_variable_dp(variable, io, name, x, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    real(8), intent(in) :: x
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
  end subroutine adios2_comm_define_variable_dp

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
  end subroutine adios2_comm_define_variable_dp_arr1d

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
  end subroutine adios2_comm_define_variable_dp_arr2d

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
  end subroutine adios2_comm_define_variable_dp_arr3d

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
  end subroutine adios2_comm_define_variable_dp_arr4d

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
  end subroutine adios2_comm_define_variable_dp_arr5d

  subroutine adios2_comm_define_variable_cdp(variable, io, name, x, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
  end subroutine adios2_comm_define_variable_cdp

  subroutine adios2_comm_define_variable_cdp_arr1d(&
       variable, io, name, x, &
       shape_dims, start_dims, count_dims, &
       is_constant_dims, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x(:)
    integer(kind=8), dimension(:), intent(in) :: shape_dims
    integer(kind=8), dimension(:), intent(in) :: start_dims
    integer(kind=8), dimension(:), intent(in) :: count_dims
    logical, intent(in) :: is_constant_dims
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
         size(shape_dims), shape_dims, start_dims, count_dims, &
         is_constant_dims, ierr)
  end subroutine adios2_comm_define_variable_cdp_arr1d

  subroutine adios2_comm_define_variable_cdp_arr2d(&
       variable, io, name, x, &
       shape_dims, start_dims, count_dims, &
       is_constant_dims, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x(:,:)
    integer(kind=8), dimension(:), intent(in) :: shape_dims
    integer(kind=8), dimension(:), intent(in) :: start_dims
    integer(kind=8), dimension(:), intent(in) :: count_dims
    logical, intent(in) :: is_constant_dims
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
         size(shape_dims), shape_dims, start_dims, count_dims, &
         is_constant_dims, ierr)
  end subroutine adios2_comm_define_variable_cdp_arr2d

  subroutine adios2_comm_define_variable_cdp_arr3d(&
       variable, io, name, x, &
       shape_dims, start_dims, count_dims, &
       is_constant_dims, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x(:,:,:)
    integer(kind=8), dimension(:), intent(in) :: shape_dims
    integer(kind=8), dimension(:), intent(in) :: start_dims
    integer(kind=8), dimension(:), intent(in) :: count_dims
    logical, intent(in) :: is_constant_dims
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
         size(shape_dims), shape_dims, start_dims, count_dims, &
         is_constant_dims, ierr)
  end subroutine adios2_comm_define_variable_cdp_arr3d

  subroutine adios2_comm_define_variable_cdp_arr4d(&
       variable, io, name, x, &
       shape_dims, start_dims, count_dims, &
       is_constant_dims, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x(:,:,:,:)
    integer(kind=8), dimension(:), intent(in) :: shape_dims
    integer(kind=8), dimension(:), intent(in) :: start_dims
    integer(kind=8), dimension(:), intent(in) :: count_dims
    logical, intent(in) :: is_constant_dims
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
         size(shape_dims), shape_dims, start_dims, count_dims, &
         is_constant_dims, ierr)
  end subroutine adios2_comm_define_variable_cdp_arr4d

  subroutine adios2_comm_define_variable_cdp_arr5d(&
       variable, io, name, x, &
        shape_dims, start_dims, count_dims, &
        is_constant_dims, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    complex(kind=8), intent(in) :: x(:,:,:,:,:)
    integer(kind=8), dimension(:), intent(in) :: shape_dims
    integer(kind=8), dimension(:), intent(in) :: start_dims
    integer(kind=8), dimension(:), intent(in) :: count_dims
    logical, intent(in) :: is_constant_dims
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), &
         size(shape_dims), shape_dims, start_dims, count_dims, &
         is_constant_dims, ierr)
  end subroutine adios2_comm_define_variable_cdp_arr5d

  subroutine adios2_comm_define_variable_integer4(variable, io, name, x, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    integer(4), intent(in) :: x
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
  end subroutine adios2_comm_define_variable_integer4

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
  end subroutine adios2_comm_define_variable_integer4_arr1d

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
  end subroutine adios2_comm_define_variable_integer4_arr2d

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
  end subroutine adios2_comm_define_variable_integer4_arr3d

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
  end subroutine adios2_comm_define_variable_integer4_arr4d

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
  end subroutine adios2_comm_define_variable_integer4_arr5d

  subroutine adios2_comm_define_variable_integer8(variable, io, name, x, ierr)
    type(adios2_variable), intent(out) :: variable
    type(adios2_io), intent(in) :: io
    character*(*), intent(in) :: name
    integer(8), intent(in) :: x
    integer, intent(out) :: ierr

    call adios2_define_variable(variable, io, name, adios2_comm_get_type(x), ierr)
  end subroutine adios2_comm_define_variable_integer8

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
  end subroutine adios2_comm_define_variable_integer8_arr1d

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
  end subroutine adios2_comm_define_variable_integer8_arr2d

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
  end subroutine adios2_comm_define_variable_integer8_arr3d

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
 end subroutine adios2_comm_define_variable_integer8_arr4d

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
  end subroutine adios2_comm_define_variable_integer8_arr5d
#endif
end module adios2_comm_module
