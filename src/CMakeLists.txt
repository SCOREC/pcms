

# TODO split out the field transfer library
set(
  PCMS_HEADERS
  pcms.h
  pcms/bounding_box.h
  pcms/coordinate.h
  pcms/coordinate_systems.h
  pcms/coordinate_transform.h
  pcms/field.h
  pcms/create_field.h
  pcms/field_communicator.h
  pcms/global_communicator.h
  pcms/field_communicator2.h
  pcms/field_evaluation_methods.h
  pcms/partition.h
  pcms/coupler.h
  pcms/coordinate_system.h
  pcms/field_layout.h
)

set(
  PCMS_SOURCES
  pcms.cpp
  pcms/create_field.cpp
  pcms/adapter/point_cloud/point_cloud_layout.cpp
  pcms/adapter/point_cloud/point_cloud.cpp
  pcms/adapter/omega_h/omega_h_field_layout.cpp
  pcms/adapter/omega_h/omega_h_field2.cpp
)
set(
  ADAPTER_HEADERS
  ${CMAKE_SOURCE_DIR}/src/pcms/adapter/point_cloud
)
configure_file(pcms/version.h.in pcms/version.h)
configure_file(pcms/configuration.h.in pcms/configuration.h)
list(APPEND PCMS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/pcms/version.h ${CMAKE_CURRENT_BINARY_DIR}/pcms/configuration.h)

add_subdirectory(pcms/utility)

if(PCMS_ENABLE_XGC)
  list(APPEND PCMS_SOURCES pcms/adapter/xgc/xgc_reverse_classification.cpp)
  list(APPEND ADAPTER_HEADERS ${CMAKE_SOURCE_DIR}/src/pcms/adapter/xgc)
endif()
if(PCMS_ENABLE_OMEGA_H)
  list(APPEND PCMS_SOURCES pcms/point_search.cpp)
  list(
    APPEND
    PCMS_HEADERS
    pcms/transfer_field.h
    pcms/transfer_field2.h
    pcms/uniform_grid.h
    pcms/point_search.h
    #install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/pcms/adapter/omega_h/omega_h_field.h
    #DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/adapter/omega_h)
  )
  list(APPEND ADAPTER_HEADERS ${CMAKE_SOURCE_DIR}/src/pcms/adapter/omega_h)
endif()

find_package(Kokkos REQUIRED)
find_package(perfstubs REQUIRED)

add_library(pcms_core ${PCMS_SOURCES})
set_target_properties(
  pcms_core PROPERTIES OUTPUT_NAME pcmscore EXPORT_NAME
  core
)
add_library(pcms::core ALIAS pcms_core)
target_compile_features(pcms_core PUBLIC cxx_std_17)
target_link_libraries(
  pcms_core PUBLIC meshfields::meshfields redev::redev
  MPI::MPI_CXX Kokkos::kokkos perfstubs pcms::utility
)
if(PCMS_ENABLE_OMEGA_H)
  target_link_libraries(pcms_core PUBLIC Omega_h::omega_h)
endif()

if(PCMS_HAS_ASAN)
  target_compile_options(
    pcms_core PRIVATE -fsanitize=address
    -fno-omit-frame-pointer
  )
endif()

# export the library
set_target_properties(pcms_core PROPERTIES PUBLIC_HEADER "${PCMS_HEADERS}")
target_include_directories(
  pcms_core
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
install(
  TARGETS pcms_core
  EXPORT pcms_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/config.cmake.in" "${CMAKE_BINARY_DIR}/pcms-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  COMPATIBILITY AnyNewerVersion
)

install(
  FILES "${PROJECT_BINARY_DIR}/pcms-config.cmake"
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)

install(
  EXPORT pcms_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)

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
target_link_libraries(pcms_pcms INTERFACE pcms::interpolator)

install(
  TARGETS pcms_pcms
  EXPORT pcms-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

# Mirror the src tree for adapter headers
foreach(hdr IN LISTS ADAPTER_HEADERS)
  install(DIRECTORY ${hdr}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/adapter
          FILES_MATCHING PATTERN "*.h")
endforeach ()

# install external headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pcms/external/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/external
)
install(
  EXPORT pcms-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
