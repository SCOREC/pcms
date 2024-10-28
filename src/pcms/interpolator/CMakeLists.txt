cmake_minimum_required(VERSION 3.20)

project(Interpolator)

set(HEADER_FILES
    adj_search_dega2.hpp
    MLSInterpolation.hpp
    points.hpp
    adj_search.hpp
    MLSCoefficients.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
)

install(FILES ${HEADER_FILES} DESTINATION include/pcms/interpolator)

add_library(pcms_interpolator INTERFACE)

target_link_libraries(pcms_interpolator INTERFACE pcms::core)

target_include_directories(pcms_interpolator INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator>)


install(TARGETS pcms_interpolator
    EXPORT interpolatorTargets
    INCLUDES DESTINATION include/pcms/interpolator)


install(EXPORT interpolatorTargets
    FILE interpolatorTargets.cmake
    NAMESPACE pcms::
    DESTINATION lib/cmake/pcms)

