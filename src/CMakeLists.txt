# TODO split out the field transfer library
set(
  PCMS_HEADERS
  pcms.h
)

set(
  PCMS_SOURCES
  pcms.cpp
)

configure_file(pcms/version.h.in pcms/version.h)
configure_file(pcms/configuration.h.in pcms/configuration.h)
list(APPEND PCMS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/pcms/version.h ${CMAKE_CURRENT_BINARY_DIR}/pcms/configuration.h)

add_subdirectory(pcms/utility)
add_subdirectory(pcms/transient)
add_subdirectory(pcms/discretization)
add_subdirectory(pcms/localization)
add_subdirectory(pcms/field)
add_subdirectory(pcms/coupler)



find_package(Kokkos REQUIRED)
find_package(perfstubs REQUIRED)

add_library(pcms_core ${PCMS_SOURCES})
set_target_properties(
  pcms_core PROPERTIES OUTPUT_NAME pcmscore EXPORT_NAME
  core
)
add_library(pcms::core ALIAS pcms_core)
target_compile_features(pcms_core PUBLIC cxx_std_20)
target_link_libraries(
  pcms_core PUBLIC redev::redev
  MPI::MPI_CXX Kokkos::kokkos perfstubs pcms::utility pcms::localization pcms::field pcms::coupler
)

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

install(FILES "${CMAKE_SOURCE_DIR}/cmake/FindPETSc.cmake"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  EXPORT pcms_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)

add_library(pcms_pcms INTERFACE)
target_link_libraries(pcms_pcms INTERFACE pcms::core pcms::transient)
set_target_properties(pcms_pcms PROPERTIES EXPORT_NAME pcms)
if (PCMS_ENABLE_Python)
  # Disable LTO/IPO before adding Python subdirectory to avoid fatbinData conflicts
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
  add_subdirectory(pcms/pythonapi)
endif()
if(PCMS_ENABLE_C)
  add_subdirectory(pcms/capi)
  target_link_libraries(pcms_pcms INTERFACE pcms::capi)
endif()
if(PCMS_ENABLE_Fortran)
  add_subdirectory(pcms/fortranapi)
  target_link_libraries(pcms_pcms INTERFACE pcms::fortranapi)
endif()

add_subdirectory(pcms/transfer)
target_link_libraries(pcms_pcms INTERFACE pcms::transfer)

install(
  TARGETS pcms_pcms
  EXPORT pcms-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)


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
