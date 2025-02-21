! Created by Fuad Hasan on 2/21/25.

program test_interpolation
    use pcms_interpolator
    implicit none

    type(PcmsInterpolatorHandle) :: interpolator
    type(PcmsInterpolatorOHMeshHandle) :: mesh
    character(len=100) :: filename
    real(8) :: radius
    radius = 0.1d0

    ! Read the mesh name from the command line
    if (command_argument_count() /= 1) then
        print *, "Usage: test_interpolation <mesh_filename>"
        stop 1
    end if

    call get_command_argument(1, filename)
    filename = trim(filename)

    print *, "Reading mesh from file: ", filename


    call pcms_kokkos_initialize_without_args()

    mesh = read_oh_mesh(filename)
    interpolator = pcms_create_interpolator(mesh, radius)


    call pcms_destroy_interpolator(interpolator)
    call release_oh_mesh(mesh)

    call pcms_kokkos_finalize()

end program test_interpolation