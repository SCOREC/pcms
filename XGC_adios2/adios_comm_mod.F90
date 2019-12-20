module adios_comm_module
    use sml_module
    use adios_read_mod
    implicit none
    character (len=10) :: staging_read_method_name
    character(len=128) :: staging_read_params = "verbose=2"
    
    integer :: staging_read_method = ADIOS_READ_METHOD_BP !! internal value

    namelist /adios_param/ staging_read_method_name, staging_read_params
contains
    subroutine adios_comm_init()
        use input_module
        implicit none
        integer :: ierr

        staging_read_method_name = 'BP'//char(0)

        open(unit=20,file='adios.in', status='old',action='read')
        READ(20, NML=adios_param)
        close(unit=20)

        staging_read_method = adios_comm_get_read_method (trim(staging_read_method_name))
        if(sml_mype==0) print *, 'staging_read_method_name=', trim(staging_read_method_name)
        if(sml_mype==0) print *, 'staging_read_method=', staging_read_method
        call adios_init('adioscfg.xml'//char(0), sml_comm, ierr)
        call adios_read_init_method(staging_read_method, sml_comm, trim(staging_read_params), ierr)
    end subroutine adios_comm_init

    subroutine adios_comm_finalize()
        implicit none
        integer :: ierr
        call adios_finalize(sml_mype,ierr)
        call adios_read_finalize_method (staging_read_method, ierr)
    end subroutine adios_comm_finalize

    function adios_comm_get_read_method(method)
        implicit none
        integer :: adios_comm_get_read_method
        character(len=*), intent(in) :: method

        adios_comm_get_read_method = ADIOS_READ_METHOD_BP
        if (trim(method) .eq. 'FLEXPATH') then
            adios_comm_get_read_method = ADIOS_READ_METHOD_FLEXPATH
        else if (trim(method) .eq. 'DATASPACES') then
            adios_comm_get_read_method = ADIOS_READ_METHOD_DATASPACES
        else if (trim(method) .eq. 'DIMES') then
            adios_comm_get_read_method = ADIOS_READ_METHOD_DIMES
        end if
    end function adios_comm_get_read_method
end module adios_comm_module
