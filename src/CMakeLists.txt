
# The version comes from git describe, so version.hpp is regenerated on every
# build, and once now so it exists after configure. See PcmsVersion.cmake.
find_package(Git QUIET)
set(
  PCMS_VERSION_ARGS
  -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
  -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
  -DPROJECT_VERSION=${PROJECT_VERSION}
  -DVERSION_OVERRIDE=${SKBUILD_PROJECT_VERSION_FULL}
  -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/pcms/version.hpp.in
  -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/pcms/version.hpp
  -P ${PROJECT_SOURCE_DIR}/cmake/PcmsVersion.cmake
)
execute_process(COMMAND ${CMAKE_COMMAND} ${PCMS_VERSION_ARGS})
add_custom_target(
  pcms_version_header
  COMMAND ${CMAKE_COMMAND} ${PCMS_VERSION_ARGS}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/pcms/version.hpp
)
configure_file(pcms/configuration.hpp.in pcms/configuration.hpp)

add_subdirectory(pcms/utility)
# the generated version.hpp is found through pcms::utility's include directories
add_dependencies(pcms_utility pcms_version_header)
add_subdirectory(pcms/transient)
add_subdirectory(pcms/geometry)
add_subdirectory(pcms/localization)
add_subdirectory(pcms/field)
if(PCMS_ENABLE_COUPLER)
  add_subdirectory(pcms/coupler)
endif()
# backend adapters; each adds its files to the core targets above
add_subdirectory(pcms/backends)
add_subdirectory(pcms/eqdsk)

# Top-level and generated headers. They are found through pcms::utility's
# include directories (source and binary src/ directories).
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/pcms/version.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/pcms/configuration.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

# the package config finds the private dependencies only for a static build
if(BUILD_SHARED_LIBS)
  set(PCMS_BUILD_SHARED TRUE)
else()
  set(PCMS_BUILD_SHARED FALSE)
endif()
configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/config.cmake.in" "${CMAKE_BINARY_DIR}/pcms-config.cmake"
  INSTALL_DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  COMPATIBILITY SameMinorVersion
)

install(
  FILES "${PROJECT_BINARY_DIR}/pcms-config.cmake"
  "${PROJECT_BINARY_DIR}/pcms-config-version.cmake"
  DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)

# FindPETSc.cmake stays installed: the package config calls
# find_dependency(PETSc) when PETSc is enabled, which needs it.
install(FILES "${CMAKE_SOURCE_DIR}/cmake/FindPETSc.cmake"
          DESTINATION ${PCMS_INSTALL_CMAKEDIR})


# pcms::pcms links every PCMS library that was built; downstream codes
# should link it.
add_library(pcms_pcms INTERFACE)
add_library(pcms::pcms ALIAS pcms_pcms)
target_link_libraries(
  pcms_pcms INTERFACE pcms::utility pcms::localization
                      pcms::field)
if(TARGET pcms::coupler)
  target_link_libraries(pcms_pcms INTERFACE pcms::coupler)
endif()
target_link_libraries(pcms_pcms INTERFACE pcms::transient)
set_target_properties(pcms_pcms PROPERTIES EXPORT_NAME pcms)
target_sources(
  pcms_pcms INTERFACE FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include
                      FILES ${PROJECT_SOURCE_DIR}/include/pcms/pcms.hpp)
if (PCMS_ENABLE_PYTHON)
  # Disable LTO/IPO before adding Python subdirectory to avoid fatbinData conflicts
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
  add_subdirectory(pcms/pythonapi)
endif()
if(PCMS_ENABLE_C)
  add_subdirectory(pcms/capi)
  target_link_libraries(pcms_pcms INTERFACE pcms::capi)
endif()
if(PCMS_ENABLE_FORTRAN)
  add_subdirectory(pcms/fortranapi)
  target_link_libraries(pcms_pcms INTERFACE pcms::fortranapi)
endif()

add_subdirectory(pcms/transfer)
target_link_libraries(pcms_pcms INTERFACE pcms::transfer)

# Library versions: the SOVERSION is major.minor, since 0.x releases may
# break the ABI in every minor release (find_package uses SameMinorVersion).
foreach(
  lib IN
  ITEMS pcms_utility pcms_geometry pcms_localization
        pcms_field
        pcms_transfer pcms_coupler pcms_transient pcms_capi_core
        pcms_capi_transfer pcms_fortranapi_core pcms_fortranapi_transfer)
  if(TARGET ${lib})
    set_target_properties(
      ${lib} PROPERTIES VERSION ${PROJECT_VERSION}
                        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
  endif()
endforeach()

# Deprecated: pcms::core used to be the coupler plus everything else. Kept
# for one release as an alias of pcms::pcms.
add_library(pcms_core INTERFACE)
add_library(pcms::core ALIAS pcms_core)
target_link_libraries(pcms_core INTERFACE pcms_pcms)
set_target_properties(
  pcms_core PROPERTIES EXPORT_NAME core DEPRECATION
                       "pcms::core is deprecated; link pcms::pcms instead.")
install(TARGETS pcms_core EXPORT pcms_core-targets)
install(
  EXPORT pcms_core-targets
  NAMESPACE pcms::
  DESTINATION ${PCMS_INSTALL_CMAKEDIR})

install(
  TARGETS pcms_pcms
  EXPORT pcms-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILE_SET HEADERS
)


install(
  EXPORT pcms-targets
  NAMESPACE pcms::
  DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
