module assert_mod
  use sml_module, only: sml_mype
  implicit none

  public

  interface assert
     module procedure   &
          assert_int1, assert_int2,  &
          assert_dble1, assert_log1, assert_char1
  end interface
  
contains

  subroutine assert_int2(lcond,msg,ival1,ival2)
    implicit none
    logical, intent(in) :: lcond
    character(len=*), intent(in) :: msg
    integer, intent(in) :: ival1,ival2

    if (.not.lcond) then
       write(*,*) sml_mype,msg,ival1,ival2
       stop '** assertion error ** '
    endif
    return
  end subroutine assert_int2

  subroutine assert_int1(lcond,msg,ival)
    implicit none
    logical, intent(in) :: lcond
    character(len=*), intent(in) :: msg
    integer, intent(in) :: ival


    if (.not.lcond) then
       write(*,*) sml_mype, msg, ival
       stop 'assertion error '
    endif
    return
  end subroutine assert_int1


  subroutine assert_dble1(lcond,msg,ival)
    implicit none
    logical, intent(in) :: lcond
    character(len=*), intent(in) :: msg
    doubleprecision, intent(in) :: ival


    if (.not.lcond) then
       write(*,*) sml_mype, msg, ival
       stop 'assertion error '
    endif
    return
  end subroutine assert_dble1


  subroutine assert_log1(lcond,msg,ival)
    implicit none
    logical, intent(in) :: lcond
    character(len=*), intent(in) :: msg
    logical, intent(in) :: ival


    if (.not.lcond) then
       write(*,*) sml_mype, msg, ival
       stop 'assertion error '
    endif
    return
  end subroutine assert_log1

  subroutine assert_char1(lcond,msg,ival)
    implicit none
    logical, intent(in) :: lcond
    character(len=*), intent(in) :: msg
    character(len=*), intent(in) :: ival


    if (.not.lcond) then
       write(*,*) sml_mype, msg, ival
       stop 'assertion error '
    endif
    return
  end subroutine assert_char1

end module assert_mod
