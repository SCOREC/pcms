module adios_comm_module

    use adios_read_mod
	use mpi

    implicit none

    character (len=10) :: staging_read_method_name
    integer :: staging_read_method = ADIOS_READ_METHOD_BP !! internal value
	character(len=128) :: staging_read_params = "verbose=2"


    namelist /adios_param/ staging_read_method_name, staging_read_params

	contains

		subroutine adios_comm_get_read_method(rank, fileunit)
			integer, intent(in), optional :: rank
			integer, intent(in), optional :: fileunit
			integer :: ierr, funit

			! Default value
			staging_read_method_name = 'BP'//char(0)

			! Use a default fileunit, but allow argument for it
			if (present(fileunit)) then
				funit = fileunit
			else
				funit = 100
			end if

			! Read from file
			open(unit=funit, file='adios.in', status='old', action='read')
			read(funit, nml=adios_param)
			close(unit=funit)

			! Change string input to ADIOS integer type
			if (trim(staging_read_method_name) .eq. 'FLEXPATH') then
				staging_read_method = ADIOS_READ_METHOD_FLEXPATH
			else if (trim(staging_read_method_name) .eq. 'DATASPACES') then
				staging_read_method = ADIOS_READ_METHOD_DATASPACES
			else if (trim(staging_read_method_name) .eq. 'DIMES') then
				staging_read_method = ADIOS_READ_METHOD_DIMES
			end if

			! Debug output
			if (present(rank)) then
				if (rank == 0) then
					print *, 'staging_read_method_name=', trim(staging_read_method_name)
					print *, 'staging_read_method=', staging_read_method
				end if
			end if

		end subroutine adios_comm_get_read_method


		subroutine adios_comm_init(xmlfile, comm)
			character(len=*), intent(in) :: xmlfile
			integer, intent(in) :: comm
			integer :: ierr

			call adios_comm_get_read_method(comm)
			call adios_init(trim(xmlfile), comm, ierr)
			call adios_read_init_method(staging_read_method, comm, 'verbose=3', ierr)

		end subroutine adios_comm_init


		subroutine adios_comm_finalize(comm)
			integer, intent(in) :: comm
			integer :: ierr, rank

			call mpi_comm_rank(comm, rank, ierr)
			call adios_finalize(rank, ierr)
			call adios_read_finalize_method(staging_read_method, ierr)
		end subroutine adios_comm_finalize


end module adios_comm_module
