set(
        PCMS_LOCALIZATION_HEADERS
        point_search.h
        queue_visited.hpp
        adj_search.hpp
        mls_support_helpers.h
        localization_factory.h
        point_cloud_localization.h
        mesh_localization.h
)

set(
        PCMS_LOCALIZATION_SOURCES
        point_search.cpp
        adj_search.cpp
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
set_target_properties(pcms_localization PROPERTIES PUBLIC_HEADER "${PCMS_LOCALIZATION_HEADERS}")

install(
        TARGETS pcms_localization
        EXPORT pcms_localization-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/localization
)

install(
        EXPORT pcms_localization-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
