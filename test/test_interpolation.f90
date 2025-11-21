! Created by Fuad Hasan on 2/21/25.

program test_interpolation
    use pcms_interpolator
    use iso_c_binding
    implicit none

    !!!!!!!!!!!!!! Declare the variables !!!!!!!!!!!!!!
    type(PcmsInterpolatorHandle) :: interpolator, point_cloud_interpolator
    type(PcmsInterpolatorOHMeshHandle) :: mesh
    character(len=100) :: filename, num_faces_str, num_vertices_str
    real(8) :: radius
    integer :: num_faces, num_vertices
    integer :: i
    ! Didn't use real(c_double) to show that it works with real(8) as well
    real(8), allocatable, target :: source_at_face(:), target_at_vertex(:)

    real(8), allocatable, target :: point_cloud_source_points(:)
    real(8), allocatable, target :: point_cloud_target_points(:)
    real(8), allocatable, target :: point_cloud_source_values(:)
    real(8), allocatable, target :: point_cloud_target_values(:)
    integer :: num_point_cloud_source_points
    integer :: num_point_cloud_target_points

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!! Read Arguments !!!!!!!!!!!!!!!!!!!!!!!!!

    radius = 0.1d0
    ! number of faces and vertices

    ! Read the mesh name from the command line
    if (command_argument_count() /= 3) then
        print *, "Usage: test_interpolation <mesh_filename> num_faces num_vertices"
        stop 1
    end if

    call get_command_argument(1, filename)
    filename = trim(filename)
    call get_command_argument(2, num_faces_str)
    call get_command_argument(3, num_vertices_str)
    read(num_faces_str, *) num_faces
    read(num_vertices_str, *) num_vertices

    print *, "Reading mesh from file: ", filename
    print *, "Number of faces:        ", num_faces
    print *, "Number of vertices:     ", num_vertices

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!! Initialize !!!!!!!!!!!!!!!!!
    call pcms_kokkos_initialize_without_args()

    mesh = read_oh_mesh(filename)
    interpolator = pcms_create_interpolator(mesh, radius)

    allocate(source_at_face(num_faces))
    allocate(target_at_vertex(num_vertices))

    do i = 1, num_faces
        source_at_face(i) = 2.0d0
    end do

    do i = 1, num_vertices
        target_at_vertex(i) = 0.0d0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!! Interpolate !!!!!!!!!!!!!!!!!!!!
    call pcms_interpolate(interpolator, c_loc(source_at_face), num_faces, c_loc(target_at_vertex), num_vertices)

    !!!!!!!!!!!!!!! Test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This is not checking the correctness of the interpolation
    ! It only shows the functionality of the interpolation fortran API

    ! print the interpolated values
    do i = 1, num_vertices
        print *, "(", i, ", ", target_at_vertex(i), ")"
        ! if values are not close to 2.0, then the interpolation is not working; exit with error
        if (abs(target_at_vertex(i) - 2.0d0) > 1.0d-4) then
            print *, "Interpolation failed, expected about 2.0, got ", target_at_vertex(i)
            stop 1
        end if
    end do

    print *, "Mesh Node to Cell Centroid Interpolation successful!"

    !!!!!!!!!!!!!!!!! 2D Point Cloud Interpolator Test !!!!!!!!!!!!!!!!!!!!
    num_point_cloud_source_points = 9 ! first order interpolation needs at least 6 points
    num_point_cloud_target_points = 4
    allocate(point_cloud_source_points(2 * num_point_cloud_source_points))
    allocate(point_cloud_target_points(2 * num_point_cloud_target_points))
    allocate(point_cloud_source_values(num_point_cloud_source_points))
    allocate(point_cloud_target_values(num_point_cloud_target_points))

    !!!!!! Geometry of the Point Clouds !!!!!!
    !       *(0,1)            *(0.5,1)            *(1,1)
    !         |                                     |
    !         |    @(.25,.75)           @(.75,.75)  |
    !       *(0,0.5)          *(0.5,0.5)          *(1,0.5)
    !         |    @(.25,.25)           @(.75,.25)  |
    !         |                                     |
    !       *(0,0)            *(0.5,0)            *(1,0)
    !     ----------------------------------------------------> X
    !    * -> Source Points
    !    @ -> Target Points
    !    All source points have value 2.0

    point_cloud_source_points = [0.0d0, 0.0d0, &
                                 1.0d0, 0.0d0, &
                                 0.0d0, 1.0d0, &
                                 1.0d0, 1.0d0, &
                                 0.5d0, 0.5d0, &
                                 0.0d0, 0.5d0, &
                                 1.0d0, 0.5d0, &
                                 0.5d0, 0.0d0, &
                                 0.5d0, 1.0d0]
    point_cloud_source_values = [2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0]
    point_cloud_target_points = [0.25d0, 0.25d0, &
                                 0.75d0, 0.25d0, &
                                 0.25d0, 0.75d0, &
                                 0.75d0, 0.75d0]

    point_cloud_interpolator = pcms_create_point_based_interpolator(c_loc(point_cloud_source_points), &
                                                                        num_point_cloud_source_points*2, &
                                                                        c_loc(point_cloud_target_points), &
                                                                        num_point_cloud_target_points*2, &
                                                                        1.0d0, & ! radius
                                                                        1, & ! degree
                                                                        6, &! min neighbors
                                                                        0.0d0, & ! lambda
                                                                        5.0d0) ! decay factor

    call pcms_interpolate(point_cloud_interpolator, c_loc(point_cloud_source_values), &
                          num_point_cloud_source_points, &
                          c_loc(point_cloud_target_values), &
                          num_point_cloud_target_points)

    ! print the interpolated values
    do i = 1, num_point_cloud_target_points
        print *, "Point Cloud Target Point ", i, " value: ", point_cloud_target_values(i)
        ! if values are not close to 2.0, then the interpolation is not working; exit with error
        if (abs(point_cloud_target_values(i) - 2.0d0) > 1.0d-6) then
            print *, "Point Cloud Interpolation failed, expected about 2.0, got ", point_cloud_target_values(i)
        end if
    end do
    print *, "2D Point Cloud Interpolation successful!"

    !!!!!!!!!!!!!!!! Destroy !!!!!!!!!!!!!!!!!!!!!!!!
    call pcms_destroy_interpolator(interpolator)
    call pcms_destroy_interpolator(point_cloud_interpolator)
    call release_oh_mesh(mesh)

    call pcms_kokkos_finalize()

    deallocate(source_at_face)
    deallocate(target_at_vertex)

end program test_interpolation
