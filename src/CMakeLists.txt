# TODO split out the field transfer library
set(PCMS_HEADERS
        pcms.h
        pcms/arrays.h
        pcms/assert.h
        pcms/bounding_box.h
        pcms/common.h
        pcms/coordinate.h
        pcms/coordinate_systems.h
        pcms/coordinate_transform.h
        pcms/field.h
        pcms/create_field.h
        pcms/field_communicator.h
        pcms/field_communicator2.h
        pcms/field_evaluation_methods.h
        pcms/memory_spaces.h
        pcms/types.h
        pcms/array_mask.h
        pcms/inclusive_scan.h
        pcms/profile.h
        pcms/print.h
        pcms/partition.h
        pcms/coupler.h
        pcms/coordinate_system.h
        pcms/field_layout.h
        pcms/adapter/point_cloud/point_cloud_layout.h
        pcms/adapter/point_cloud/point_cloud.h
        pcms/adapter/omega_h/omega_h_field_layout.h
        pcms/adapter/omega_h/omega_h_field2.h
        )

set(PCMS_SOURCES
        pcms.cpp
        pcms/assert.cpp
        pcms/print.cpp
        pcms/create_field.cpp
        pcms/adapter/xgc/xgc_field_adapter.h
        pcms/adapter/point_cloud/point_cloud_layout.cpp
        pcms/adapter/point_cloud/point_cloud.cpp
        pcms/adapter/omega_h/omega_h_field_layout.cpp
        pcms/adapter/omega_h/omega_h_field2.cpp
        pcms/adapter/xgc/xgc_field_adapter.h)

configure_file(pcms/version.h.in pcms/version.h)

if(PCMS_ENABLE_XGC)
  list(APPEND PCMS_SOURCES pcms/adapter/xgc/xgc_reverse_classification.cpp)
  list(APPEND PCMS_HEADERS pcms/adapter/xgc/xgc_reverse_classification.h)
endif()
if(PCMS_ENABLE_OMEGA_H)
  list(APPEND PCMS_SOURCES pcms/point_search.cpp)
  list(APPEND PCMS_HEADERS
          pcms/adapter/omega_h/omega_h_field.h
          pcms/transfer_field.h
          pcms/transfer_field2.h
          pcms/uniform_grid.h
          pcms/point_search.h)
endif ()

find_package(Kokkos REQUIRED)
find_package(perfstubs REQUIRED)

add_library(pcms_core ${PCMS_SOURCES})
set_target_properties(pcms_core PROPERTIES OUTPUT_NAME pcmscore EXPORT_NAME
                                                                core)
add_library(pcms::core ALIAS pcms_core)
target_compile_features(pcms_core PUBLIC cxx_std_17)
target_link_libraries(pcms_core PUBLIC meshfields::meshfields redev::redev MPI::MPI_CXX Kokkos::kokkos perfstubs)
if(PCMS_ENABLE_OMEGA_H)
  target_link_libraries(pcms_core PUBLIC Omega_h::omega_h)
  target_compile_definitions(pcms_core PUBLIC -DPCMS_HAS_OMEGA_H)
endif()
if(PCMS_ENABLE_SERVER)
  target_compile_definitions(pcms_core PUBLIC -DPCMS_HAS_SERVER)
endif()

if(PCMS_HAS_ASAN)
  target_compile_options(pcms_core PRIVATE -fsanitize=address
                                           -fno-omit-frame-pointer)
endif()

if(PCMS_PRINT_ENABLED)
  add_definitions(-DPCMS_PRINT_ENABLED)
  target_compile_definitions(pcms_core INTERFACE -DPCMS_PRINT_ENABLED)
endif()

if(spdlog_FOUND)
  add_definitions(-DPCMS_SPDLOG_ENABLED)
  target_compile_definitions(pcms_core INTERFACE -DPCMS_SPDLOG_ENABLED)
  target_link_libraries(pcms_core PUBLIC spdlog::spdlog)
endif()

# export the library
set_target_properties(pcms_core PROPERTIES PUBLIC_HEADER "${PCMS_HEADERS}")
target_include_directories(
  pcms_core
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
         "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
install(
  TARGETS pcms_core
  EXPORT pcms_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/config.cmake.in" "${CMAKE_BINARY_DIR}/pcms-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

install(FILES "${PROJECT_BINARY_DIR}/pcms-config.cmake"
              "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  EXPORT pcms_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

add_library(pcms_pcms INTERFACE)
target_link_libraries(pcms_pcms INTERFACE pcms::core)
set_target_properties(pcms_pcms PROPERTIES EXPORT_NAME pcms)
if(PCMS_ENABLE_C)
  add_subdirectory(pcms/capi)
  target_link_libraries(pcms_pcms INTERFACE pcms::capi)
endif()
if(PCMS_ENABLE_Fortran)
  add_subdirectory(pcms/fortranapi)
  target_link_libraries(pcms_pcms INTERFACE pcms::fortranapi)
endif()

add_subdirectory(pcms/interpolator)

install(
  TARGETS pcms_pcms
  EXPORT pcms-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)
# install external headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pcms/external/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/external)
install(
  EXPORT pcms-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
