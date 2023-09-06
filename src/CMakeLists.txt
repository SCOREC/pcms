# TODO split out the field transfer library
set(WDMCPL_HEADERS
        pcms.h
        pcms/arrays.h
        pcms/assert.h
        pcms/bounding_box.h
        pcms/common.h
        pcms/coordinate.h
        pcms/coordinate_systems.h
        pcms/coordinate_transform.h
        pcms/field.h
        pcms/field_communicator.h
        pcms/field_evaluation_methods.h
        pcms/memory_spaces.h
        pcms/types.h
        pcms/array_mask.h
        pcms/inclusive_scan.h
        pcms/profile.h
        )

set(WDMCPL_SOURCES
        pcms.cpp
        pcms/assert.cpp
        pcms/xgc_field_adapter.h)
set(WDMCPL_SOURCES pcms.cpp pcms/assert.cpp)
if(WDMCPL_ENABLE_XGC)
  list(APPEND WDMCPL_SOURCES  pcms/xgc_reverse_classification.cpp)
  list(APPEND WDMCPL_HEADERS pcms/xgc_reverse_classification.h)
endif()
if (WDMCPL_ENABLE_OMEGA_H)
  list(APPEND WDMCPL_SOURCES pcms/point_search.cpp)
  list(APPEND WDMCPL_HEADERS
          pcms/omega_h_field.h
          pcms/transfer_field.h
          pcms/uniform_grid.h
          pcms/point_search.h)
endif ()
if (WDMCPL_ENABLE_SERVER)
  if (NOT WDMCPL_ENABLE_OMEGA_H)
    message(ERROR "WDMCPL_ENABLE_OMEGA_H is required for server implementation")
  endif ()
  list(APPEND WDMCPL_HEADERS pcms/server.h)
endif ()

find_package(Kokkos REQUIRED)
find_package(perfstubs REQUIRED)

add_library(pcms_core ${WDMCPL_SOURCES})
set_target_properties(pcms_core PROPERTIES
        OUTPUT_NAME pcmscore
        EXPORT_NAME core)
add_library(pcms::core ALIAS pcms_core)
target_compile_features(pcms_core PUBLIC cxx_std_17)
target_link_libraries(pcms_core PUBLIC redev::redev MPI::MPI_CXX Kokkos::kokkos perfstubs)
if(WDMCPL_ENABLE_OMEGA_H)
  target_link_libraries(pcms_core PUBLIC Omega_h::omega_h)
  target_compile_definitions(pcms_core PUBLIC -DWDMCPL_HAS_OMEGA_H)
endif()
if(WDMCPL_ENABLE_SERVER)
  target_compile_definitions(pcms_core PUBLIC -DWDMCPL_HAS_SERVER)
endif()
if(WDMCPL_ENABLE_CLIENT)
  list(APPEND WDMCPL_HEADERS pcms/client.h)
  target_compile_definitions(pcms_core PUBLIC -DWDMCPL_HAS_CLIENT)
endif()

if(WDMCPL_HAS_ASAN)
  target_compile_options(pcms_core PRIVATE -fsanitize=address -fno-omit-frame-pointer)
endif()


## export the library
set_target_properties(pcms_core PROPERTIES
  PUBLIC_HEADER "${WDMCPL_HEADERS}")
target_include_directories(pcms_core
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
install(
  TARGETS pcms_core
  EXPORT pcms_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_BINARY_DIR}/pcms-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/pcms-config.cmake"
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  EXPORT pcms_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

add_library(pcms_pcms INTERFACE)
target_link_libraries(pcms_pcms INTERFACE pcms::core)
set_target_properties(pcms_pcms PROPERTIES
        EXPORT_NAME pcms)
if (WDMCPL_ENABLE_C)
  add_subdirectory(pcms/capi)
  target_link_libraries(pcms_pcms INTERFACE pcms::capi)
endif ()
if(WDMCPL_ENABLE_Fortran)
  add_subdirectory(pcms/fortranapi)
  target_link_libraries(pcms_pcms INTERFACE pcms::fortranapi)
endif()

install(
        TARGETS pcms_pcms
        EXPORT pcms-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)
# install external headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pcms/external/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/external)
install(
  EXPORT pcms-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
