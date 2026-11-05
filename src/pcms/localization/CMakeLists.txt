set(
        PCMS_LOCALIZATION_HEADERS
        point_search.hpp
        support_results.hpp
        mls_support_helpers.hpp
        localization_factory.hpp
        ../../../include/pcms/localization/mls_options.hpp
        mesh_localization.hpp
)

set(
        PCMS_LOCALIZATION_SOURCES
        point_search.cpp
        mls_support_helpers.cpp
        mesh_localization.cpp
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
        pcms::discretization
        Omega_h::omega_h
)
target_compile_features(pcms_localization PUBLIC cxx_std_20)

set_target_properties(
        pcms_localization PROPERTIES OUTPUT_NAME pcmslocalization EXPORT_NAME
        localization
)

## export the library
# Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(PCMS_LOCALIZATION_HEADERS_PUBLIC ${PCMS_LOCALIZATION_HEADERS})
list(FILTER PCMS_LOCALIZATION_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(PCMS_LOCALIZATION_HEADERS_PRIVATE ${PCMS_LOCALIZATION_HEADERS})
list(FILTER PCMS_LOCALIZATION_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_localization
  PUBLIC FILE_SET localization TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${PCMS_LOCALIZATION_HEADERS_PUBLIC}
  PRIVATE ${PCMS_LOCALIZATION_HEADERS_PRIVATE})

install(
        TARGETS pcms_localization
        EXPORT pcms_localization-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILE_SET localization DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

install(
        EXPORT pcms_localization-targets
        NAMESPACE pcms::
        DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
