find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)
message(STATUS "Python_EXECUTABLE ${Python_EXECUTABLE}")

pybind11_add_module(pcms 
  pythonapi.cpp
  bind_field_base.cpp
  bind_field_layout.cpp
  bind_omega_h_field2.cpp
  bind_omega_h_field_layout.cpp
  bind_uniform_grid_field.cpp
  bind_uniform_grid_field_layout.cpp
  bind_omega_h.cpp
  bind_transfer_field2.cpp
        bind_mls_interpolation.cpp
        bind_mesh_utilities.cpp
)
set_target_properties(pcms PROPERTIES CXX_STANDARD 17)
target_link_libraries(pcms PRIVATE Kokkos::kokkos pybind11::module MPI::MPI_C pcms::core pcms::interpolator)



install(
        TARGETS pcms
        EXPORT pcms_python-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/pythonapi/)
install(
        EXPORT pcms_python-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)