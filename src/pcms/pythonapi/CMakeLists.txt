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
target_link_libraries(pcms PRIVATE Kokkos::kokkos pybind11::module MPI::MPI_C pcms::core pcms::transfer)

# if we are building as a python package from the pyproject.toml
# scikit-build-core will set this variable
if(SKBUILD)
  set(PCMS_PYTHON_INSTALL_LOCATION "${SKBUILD_PLATLIB_DIR}")
else()
  # this is useful for direct CMake installation into python.
  # This enables spack environment installation to directly load the PyOmega_h module
  set(PCMS_PYTHON_INSTALL_LOCATION "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
  # Ensure the extension can resolve pcms shared libraries from a standard
  # prefix install layout (e.g. Spack), where the module lives under
  # <prefix>/<libdir>/pythonX.Y/site-packages and core libs live in <libdir>.
  set_target_properties(pcms PROPERTIES INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../../../")
endif()


install(
        TARGETS pcms
	RUNTIME DESTINATION ${PCMS_PYTHON_INSTALL_LOCATION}
	LIBRARY DESTINATION ${PCMS_PYTHON_INSTALL_LOCATION})
