# no public headers: the MLS support search is private
set(
        PCMS_LOCALIZATION_HEADERS
        mls_support_helpers.hpp
)

set(
        PCMS_LOCALIZATION_SOURCES
        mls_support_helpers.cpp
)
add_library(pcms_localization ${PCMS_LOCALIZATION_SOURCES})
add_library(pcms::localization ALIAS pcms_localization)
target_include_directories(
        pcms_localization
        PUBLIC
        # include path should be pcms/localization/
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        # include path to search for pcms/config.h
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(
        pcms_localization PUBLIC
        Kokkos::kokkos
        pcms::utility
        pcms::geometry
)
# Omega_h is used by the implementation only (no public header includes it)
if(PCMS_ENABLE_OMEGA_H)
  target_link_libraries(pcms_localization PRIVATE Omega_h::omega_h)
endif()
target_compile_features(pcms_localization PUBLIC cxx_std_20)

set_target_properties(
        pcms_localization PROPERTIES OUTPUT_NAME pcmslocalization EXPORT_NAME
        localization
)

## export the library
target_sources(pcms_localization PRIVATE ${PCMS_LOCALIZATION_HEADERS})

install(
        TARGETS pcms_localization
        EXPORT pcms_localization-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
        EXPORT pcms_localization-targets
        NAMESPACE pcms::
        DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
