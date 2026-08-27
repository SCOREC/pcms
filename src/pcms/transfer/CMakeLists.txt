set(PCMS_FIELD_TRANSFER_HEADERS
    linear_interpolant.hpp
    multidimarray.hpp
    interpolation_base.h
    interpolation_helpers.h
    mesh_intersection.hpp
    interpolator.h
    copy.h
    transfer_operator.hpp
    mass_matrix_type.hpp
    monte_carlo_sampling.hpp
)


set(PCMS_FIELD_TRANSFER_SOURCES
    mesh_intersection.cpp
    interpolation_base.cpp)


if(PCMS_ENABLE_MESHFIELDS)
  list(APPEND PCMS_FIELD_TRANSFER_HEADERS
      mass_matrix_integrator.hpp)
endif()

if(PCMS_ENABLE_PETSC AND PCMS_ENABLE_MESHFIELDS)
  list(APPEND PCMS_FIELD_TRANSFER_HEADERS
      linear_form_integrator.hpp
      bilinear_form_integrator.hpp
      conservative_projection_solver.hpp
      omega_h_conservative_projection.hpp
      omega_h_intersection_rhs_integrator.hpp
      omega_h_mass_integrator.hpp
      omega_h_mc_rhs_integrator.hpp
      omega_h_control_variate_projection.hpp
      coo_assembly_utils.hpp
      petsc_utils.hpp
      omega_h_form_integrator_utils.hpp)
  list(APPEND PCMS_FIELD_TRANSFER_SOURCES
      conservative_projection_solver.cpp
      omega_h_conservative_projection.cpp
      omega_h_intersection_rhs_integrator.cpp
      omega_h_mass_integrator.cpp
      omega_h_mc_rhs_integrator.cpp
      omega_h_control_variate_projection.cpp
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
target_compile_features(pcms_transfer PUBLIC cxx_std_20)

target_link_libraries(pcms_transfer PUBLIC
    pcms::core
    pcms::field
    pcms::localization
    Omega_h::omega_h)

if(PCMS_ENABLE_MESHFIELDS)
  target_link_libraries(pcms_transfer PRIVATE meshfields::meshfields)
endif()

if(PCMS_ENABLE_PETSC AND PCMS_ENABLE_MESHFIELDS)
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
