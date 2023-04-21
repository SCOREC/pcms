# TODO split out the field transfer library
set(WDMCPL_HEADERS
        wdmcpl.h
        wdmcpl/arrays.h
        wdmcpl/assert.h
        wdmcpl/bounding_box.h
        wdmcpl/common.h
        wdmcpl/coordinate.h
        wdmcpl/coordinate_systems.h
        wdmcpl/coordinate_transform.h
        wdmcpl/field.h
        wdmcpl/field_communicator.h
        wdmcpl/field_evaluation_methods.h
        wdmcpl/memory_spaces.h
        wdmcpl/types.h
        wdmcpl/array_mask.h
        )

set(WDMCPL_SOURCES
        wdmcpl.cpp
        wdmcpl/assert.cpp
        wdmcpl/xgc_field_adapter.h)
set(WDMCPL_SOURCES wdmcpl.cpp wdmcpl/assert.cpp)
if(WDMCPL_ENABLE_XGC)
  list(APPEND WDMCPL_SOURCES  wdmcpl/xgc_reverse_classification.cpp)
  list(APPEND WDMCPL_HEADERS wdmcpl/xgc_reverse_classification.h)
endif()
if (WDMCPL_ENABLE_OMEGA_H)
  list(APPEND WDMCPL_SOURCES wdmcpl/point_search.cpp)
  list(APPEND WDMCPL_HEADERS
          wdmcpl/omega_h_field.h
          wdmcpl/transfer_field.h
          wdmcpl/uniform_grid.h
          wdmcpl/point_search.h)
endif ()
if (WDMCPL_ENABLE_SERVER)
  if (NOT WDMCPL_ENABLE_OMEGA_H)
    message(ERROR "WDMCPL_ENABLE_OMEGA_H is required for server implementation")
  endif ()
  list(APPEND WDMCPL_HEADERS wdmcpl/server.h)
endif ()

find_package(Kokkos REQUIRED)

add_library(wdmcpl_core ${WDMCPL_SOURCES})
set_target_properties(wdmcpl_core PROPERTIES
        OUTPUT_NAME wdmcplcore
        EXPORT_NAME core)
add_library(wdmcpl::core ALIAS wdmcpl_core)
target_compile_features(wdmcpl_core PUBLIC cxx_std_17)
target_link_libraries(wdmcpl_core PUBLIC redev::redev MPI::MPI_CXX Kokkos::kokkos)
if(WDMCPL_ENABLE_OMEGA_H)
  target_link_libraries(wdmcpl_core PUBLIC Omega_h::omega_h)
  target_compile_definitions(wdmcpl_core PUBLIC -DWDMCPL_HAS_OMEGA_H)
endif()
if(WDMCPL_ENABLE_SERVER)
  target_compile_definitions(wdmcpl_core PUBLIC -DWDMCPL_HAS_SERVER)
endif()
if(WDMCPL_ENABLE_CLIENT)
  list(APPEND WDMCPL_HEADERS wdmcpl/client.h)
  target_compile_definitions(wdmcpl_core PUBLIC -DWDMCPL_HAS_CLIENT)
endif()

if(WDMCPL_HAS_ASAN)
  target_compile_options(wdmcpl_core PRIVATE -fsanitize=address -fno-omit-frame-pointer)
endif()


## export the library
set_target_properties(wdmcpl_core PROPERTIES
  PUBLIC_HEADER "${WDMCPL_HEADERS}")
target_include_directories(wdmcpl_core
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
install(
  TARGETS wdmcpl_core
  EXPORT wdmcpl_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wdmcpl)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_BINARY_DIR}/wdmcpl-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/wdmcpl-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/wdmcpl-config.cmake"
  "${PROJECT_BINARY_DIR}/wdmcpl-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)

install(
  EXPORT wdmcpl_core-targets
  NAMESPACE wdmcpl::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)

add_library(wdmcpl_wdmcpl INTERFACE)
target_link_libraries(wdmcpl_wdmcpl INTERFACE wdmcpl::core)
set_target_properties(wdmcpl_wdmcpl PROPERTIES
        EXPORT_NAME wdmcpl)
if (WDMCPL_ENABLE_C)
  add_subdirectory(wdmcpl/capi)
  target_link_libraries(wdmcpl_wdmcpl INTERFACE wdmcpl::capi)
endif ()
if(WDMCPL_ENABLE_Fortran)
  add_subdirectory(wdmcpl/fortranapi)
  target_link_libraries(wdmcpl_wdmcpl INTERFACE wdmcpl::fortranapi)
endif()

install(
        TARGETS wdmcpl_wdmcpl
        EXPORT wdmcpl-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wdmcpl)
install(
  EXPORT wdmcpl-targets
  NAMESPACE wdmcpl::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)