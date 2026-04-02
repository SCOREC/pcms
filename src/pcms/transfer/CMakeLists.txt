find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_TRANSFER_HEADERS
    pcms_interpolator_aliases.hpp
    adj_search.hpp
    mls_interpolation_impl.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
    mls_interpolation.hpp
    pcms_interpolator_view_utils.hpp
    pcms_interpolator_logger.hpp
    interpolation_base.h
    interpolation_helpers.h
    spline_interpolator.hpp
    load_vector_integrator.hpp
    mass_matrix_integrator.hpp
    mesh_intersection.hpp
)


set(PCMS_FIELD_TRANSFER_SOURCES
    adj_search.cpp
    load_vector_integrator.cpp
    mesh_intersection.cpp
    mls_interpolation.cpp
    interpolation_base.cpp)

if(PCMS_ENABLE_PETSC)
  list(APPEND PCMS_FIELD_TRANSFER_HEADERS
      calculate_load_vector.hpp
      calculate_mass_matrix.hpp
      conservative_projection_solver.hpp
      coo_assembly_utils.hpp
      petsc_utils.hpp)
  list(APPEND PCMS_FIELD_TRANSFER_SOURCES
      calculate_load_vector.cpp
      calculate_mass_matrix.cpp
      conservative_projection_solver.cpp
      petsc_utils.cpp
      coo_assembly_utils.cpp)
endif()

add_library(pcms_transfer ${PCMS_FIELD_TRANSFER_SOURCES})
set_target_properties(pcms_transfer PROPERTIES
        OUTPUT_NAME pcmstransfer
        EXPORT_NAME transfer)
target_sources(pcms_transfer PUBLIC
        FILE_SET transfer
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_FIELD_TRANSFER_HEADERS})
add_library(pcms::transfer ALIAS pcms_transfer)
target_compile_features(pcms_transfer PUBLIC cxx_std_17)

target_link_libraries(pcms_transfer PUBLIC
    pcms::core
    Omega_h::omega_h
    PRIVATE Kokkos::kokkoskernels)

target_link_libraries(pcms_transfer PUBLIC meshfields::meshfields)

if(PCMS_ENABLE_PETSC)
  target_link_libraries(pcms_transfer PRIVATE PETSc::PETSc)
endif()

target_compile_definitions(pcms_transfer PUBLIC R3D_USE_KOKKOS)

target_include_directories(pcms_transfer INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
  TARGETS pcms_transfer
  EXPORT pcms_transfer-targets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/transfer
  FILE_SET transfer DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    EXPORT pcms_transfer-targets
    NAMESPACE pcms::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
