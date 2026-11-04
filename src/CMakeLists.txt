find_package(perfstubs REQUIRED)

configure_file(pcms/version.hpp.in pcms/version.hpp)
configure_file(pcms/configuration.hpp.in pcms/configuration.hpp)

add_subdirectory(pcms/utility)
if(PCMS_ENABLE_TRANSIENT)
  add_subdirectory(pcms/transient)
endif()
add_subdirectory(pcms/discretization)
add_subdirectory(pcms/localization)
add_subdirectory(pcms/field)
if(PCMS_ENABLE_COUPLER)
  add_subdirectory(pcms/coupler)
endif()
# backend adapters; each adds its files to the core targets above
add_subdirectory(pcms/backends)

# Top-level and generated headers. They are found through pcms::utility's
# include directories (source and binary src/ directories).
install(
  FILES ../include/pcms/pcms.hpp ${CMAKE_CURRENT_BINARY_DIR}/pcms/version.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/pcms/configuration.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

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

install(FILES "${CMAKE_SOURCE_DIR}/cmake/FindPETSc.cmake"
          DESTINATION ${PCMS_INSTALL_CMAKEDIR})


# pcms::pcms links every PCMS library that was built; downstream codes
# should link it.
add_library(pcms_pcms INTERFACE)
add_library(pcms::pcms ALIAS pcms_pcms)
target_link_libraries(
  pcms_pcms INTERFACE pcms::utility pcms::discretization pcms::localization
                      pcms::field)
if(TARGET pcms::coupler)
  target_link_libraries(pcms_pcms INTERFACE pcms::coupler)
endif()
if(PCMS_ENABLE_TRANSIENT)
  target_link_libraries(pcms_pcms INTERFACE pcms::transient)
endif()
set_target_properties(pcms_pcms PROPERTIES EXPORT_NAME pcms)
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
  ITEMS pcms_utility pcms_discretization pcms_localization pcms_field
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
)


# install external headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pcms/external/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/external
)
install(
  EXPORT pcms-targets
  NAMESPACE pcms::
  DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
