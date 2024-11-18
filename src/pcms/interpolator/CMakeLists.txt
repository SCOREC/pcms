cmake_minimum_required(VERSION 3.20)

project(Interpolator)

set(PCMS_FIELD_TRANSFER_HEADERS
    MLSInterpolation.hpp
    points.hpp
    adj_search.hpp
    MLSCoefficients.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
    MLS_rbf_options.hpp
)

set(PCMS_FIELD_TRANSFER_SOURCES
    MLS_rbf_options.cpp
)
install(FILES ${PCMS_FIELD_TRANSFER_HEADERS} DESTINATION include/pcms/interpolator)


add_library(pcms_interpolator ${PCMS_FIELD_TRANSFER_SOURCES})

target_link_libraries(pcms_interpolator PUBLIC pcms::core)

target_include_directories(pcms_interpolator INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator>)


install(TARGETS pcms_interpolator
    EXPORT interpolatorTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}) 


install(EXPORT interpolatorTargets
    FILE interpolatorTargets.cmake
    NAMESPACE pcms::
    DESTINATION lib/cmake/pcms)

