set(WDMCPL_HEADERS
  wdmcpl.h
  )

set(WDMCPL_SOURCES  wdmcpl.cpp)

add_library(wdmcpl ${WDMCPL_SOURCES})
add_library(wdmcpl::wdmcpl ALIAS wdmcpl)
target_compile_features(wdmcpl PUBLIC cxx_std_20)
target_link_libraries(wdmcpl PUBLIC Omega_h::omega_h redev::redev MPI::MPI_CXX)

find_package(Kokkos REQUIRED)


add_library(wdmcpl_static wdmcpl/wdmcpl_static.cpp wdmcpl/utility.h)
add_library(wdmcpl::wdmcpl_static ALIAS wdmcpl_static)
target_compile_features(wdmcpl_static PUBLIC cxx_std_20)
target_link_libraries(wdmcpl_static PUBLIC Omega_h::omega_h redev::redev MPI::MPI_CXX Kokkos::kokkos)
target_include_directories(wdmcpl_static
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

#add_library(wdmcpl_dynamic wdmcpl/wdmcpl_dynamic.cpp wdmcpl/utility.h wdmcpl/wdmcpl_dynamic.cpp wdmcpl/wdmcpl_dynamic.h)
#add_library(wdmcpl::wdmcpl_dynamic ALIAS wdmcpl_dynamic)
#target_compile_features(wdmcpl_dynamic PUBLIC cxx_std_20)
#target_link_libraries(wdmcpl_dynamic PUBLIC Omega_h::omega_h redev::redev MPI::MPI_CXX)
#target_include_directories(wdmcpl_dynamic
#        PUBLIC
#        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
#        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")


if(WDMCPL_HAS_ASAN)
  target_compile_options(wdmcpl PRIVATE -fsanitize=address -fno-omit-frame-pointer)
endif()

## export the library
set_target_properties(wdmcpl PROPERTIES
  PUBLIC_HEADER "${WDMCPL_HEADERS}")
target_include_directories(wdmcpl
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
install(
  TARGETS wdmcpl
  EXPORT wdmcpl-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
  EXPORT wdmcpl-targets
  NAMESPACE wdmcpl::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)
