set(
        PCMS_GEOMETRY_HEADERS
        geometry.hpp
        bound_source.hpp
        supports.hpp
)

set(
        PCMS_GEOMETRY_SOURCES
        geometry.cpp
        bound_source.cpp
        # Geometry::Bind binds coordinate maps; the map sits below the geometry
        # (its header stays in pcms/field), so it is compiled here
        ../field/coordinate_map.cpp
)
add_library(pcms_geometry ${PCMS_GEOMETRY_SOURCES})
add_library(pcms::geometry ALIAS pcms_geometry)
target_sources(pcms_geometry PUBLIC
        FILE_SET geometry
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_GEOMETRY_HEADERS})
target_include_directories(
        pcms_geometry
        PUBLIC
        # include path should be pcms/geometry/
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        # include path to search for pcms/config.h
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(
        pcms_geometry PUBLIC
        Kokkos::kokkos
        pcms::utility
        pcms::discretization
        Omega_h::omega_h
)
target_compile_features(pcms_geometry PUBLIC cxx_std_20)

set_target_properties(
        pcms_geometry PROPERTIES OUTPUT_NAME pcmsgeometry EXPORT_NAME
        geometry
)

## export the library
install(
        TARGETS pcms_geometry
        EXPORT pcms_geometry-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILE_SET geometry DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

install(
        EXPORT pcms_geometry-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
