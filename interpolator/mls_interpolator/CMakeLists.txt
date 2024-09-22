cmake_minimum_required(VERSION 3.20)
project(Interpolator)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

set(HEADER_FILES
    adj_search_dega2.hpp
    MLSInterpolation.hpp
    points.hpp
    adj_search.hpp
    MLSCoefficients.hpp
    queue_visited.hpp
)

install(FILES ${HEADER_FILES} DESTINATION include/interpolator)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/interpolatorConfigVersion.cmake"
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}
    COMPATIBILITY SameMajorVersion
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/interpolatorConfigVersion.cmake"
    DESTINATION lib/cmake/interpolator)

add_library(interpolator INTERFACE)
target_include_directories(interpolator INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/interpolator>)

install(TARGETS interpolator
    EXPORT interpolatorTargets
    INCLUDES DESTINATION include/interpolator)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/interpolatorConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/interpolatorConfig.cmake"
    INSTALL_DESTINATION lib/cmake/interpolator
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/interpolatorConfig.cmake"
    DESTINATION lib/cmake/interpolator)

install(EXPORT interpolatorTargets
    FILE interpolatorTargets.cmake
    NAMESPACE interpolator::
    DESTINATION lib/cmake/interpolator)

