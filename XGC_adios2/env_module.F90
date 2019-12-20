        module env_module
        implicit none

        integer, parameter :: idebug = 0

!       ------------------------------------------
!       module to extract environment variables
!       the original or default value is unchanged
!       if the environment variable is not present
!       ------------------------------------------

        interface get_env_variable
        module procedure                                                   &
     &         get_env_variable_i, get_env_variable_l,                     &
     &         get_env_variable_d, get_env_variable_s
        end interface get_env_variable

        contains

        subroutine get_env_variable_i(cname,value)
        implicit none
        character*(*) , intent(in) :: cname
        integer, intent(inout) :: value

        integer :: original_value

        logical :: isok
        integer :: ier,ioerr
        character*255 :: str

        original_value = value

        call get_environment_variable(name=cname,value=str,status=ier)

        isok = (ier.eq.0)
        if (isok) then
          read(str,*,iostat=ioerr) value
          isok = (ioerr.eq.0)
        endif

        if (.not.isok) then

           if (idebug.ge.1) then
             write(*,*) 'name=',trim(cname)
             write(*,*) 'str=',trim(str)
             write(*,*) 'value=',value,' original_value',original_value
          endif

           value = original_value

        endif
        return
        end subroutine get_env_variable_i



        subroutine get_env_variable_l(cname,value)
        implicit none
        character*(*) , intent(in) :: cname
        logical, intent(inout) :: value

        logical :: original_value

        logical :: isok
        integer :: ier,ioerr
        character*255 :: str

        original_value = value

        call get_environment_variable(name=cname,value=str,status=ier)

        isok = (ier.eq.0)
        if (isok) then
          value = (str(1:1) .eq. 'T').or.(str(1:1).eq.'t')
        endif

        if (.not.isok) then
!          -----------------------------------------------
!          true if value='t' or value='T', false otherwise
!          -----------------------------------------------

           if (idebug.ge.1) then
             write(*,*) 'name=',trim(cname)
             write(*,*) 'str=',trim(str)
             write(*,*) 'value=',value,' original_value',original_value
          endif

           value = original_value
        endif
        return
        end subroutine get_env_variable_l



        subroutine get_env_variable_d(cname,value)
        implicit none
        character*(*) , intent(in) :: cname
        integer, parameter :: dp = kind(1.0d0)
        real(kind=dp), intent(inout) :: value

        real(kind=dp) :: original_value

        logical :: isok
        integer :: ier,ioerr
        character*255 :: str

        original_value = value

        call get_environment_variable(name=cname,value=str,status=ier)

        isok = (ier.eq.0)
        if (isok) then
          read(str,*,iostat=ioerr) value
          isok = (ioerr.eq.0)
        endif

        if (.not.isok) then

           if (idebug.ge.1) then
             write(*,*) 'name=',trim(cname)
             write(*,*) 'str=',trim(str)
             write(*,*) 'value=',value,' original_value',original_value
          endif

           value = original_value
        endif
        return
        end subroutine get_env_variable_d




        subroutine get_env_variable_s(cname,value)
        implicit none
        character*(*) , intent(in) :: cname
        integer, parameter :: sp = kind(1.0)
        real(kind=sp), intent(inout) :: value

        real(kind=sp) :: original_value

        logical :: isok
        integer :: ier,ioerr
        character*255 :: str

        original_value = value

        call get_environment_variable(name=cname,value=str,status=ier)

        isok = (ier.eq.0)
        if (isok) then
          read(str,*,iostat=ioerr) value
          isok = (ioerr.eq.0)
        endif

        if (.not.isok) then

           if (idebug.ge.1) then
             write(*,*) 'name=',trim(cname)
             write(*,*) 'str=',trim(str)
             write(*,*) 'value=',value,' original_value',original_value
          endif

           value = original_value
        endif
        return
        end subroutine get_env_variable_s

        end module env_module
