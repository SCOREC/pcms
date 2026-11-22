set(
        PCMS_GEOMETRY_HEADERS
        ../../../include/pcms/geometry/geometry.hpp
        ../../../include/pcms/geometry/bound_points.hpp
        ../../../include/pcms/geometry/quadrature.hpp
        ../../../include/pcms/geometry/supports.hpp
        ../../../include/pcms/geometry/product_geometry.hpp
        ../../../include/pcms/geometry/detail/bound_points.hpp
        ../../../include/pcms/geometry/detail/geometry.hpp
        ../../../include/pcms/geometry/detail/quadrature.hpp
)

set(
        PCMS_GEOMETRY_SOURCES
        geometry.cpp
        bound_points.cpp
        # Geometry::Bind binds coordinate maps; the map sits below the geometry
        # (its header stays in pcms/field), so it is compiled here
        ../field/coordinate_map.cpp
        ../field/value_transform.cpp
        ../field/units.cpp
)
add_library(pcms_geometry ${PCMS_GEOMETRY_SOURCES})
add_library(pcms::geometry ALIAS pcms_geometry)
target_sources(pcms_geometry PUBLIC
        FILE_SET geometry
        TYPE HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
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
)
if(PCMS_ENABLE_UNITS)
  target_link_libraries(pcms_geometry PRIVATE units::units)
endif()
# Omega_h is used by the implementation only (no public header includes it)
if(PCMS_ENABLE_OMEGA_H)
  target_link_libraries(pcms_geometry PRIVATE Omega_h::omega_h)
endif()
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
        DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
