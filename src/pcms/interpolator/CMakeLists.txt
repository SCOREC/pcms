#cmake_minimum_required(VERSION 3.20)
#project(Interpolator)

find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_TRANSFER_HEADERS
    MLSInterpolation.hpp
    points.hpp
    adj_search.hpp
    MLSCoefficients.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
    MLS_rbf_options.hpp
        interpolation_base.h
)

set(PCMS_FIELD_TRANSFER_SOURCES
    MLS_rbf_options.cpp
        interpolation_base.cpp
)
install(FILES ${PCMS_FIELD_TRANSFER_HEADERS} DESTINATION include/pcms/interpolator)


add_library(pcms_interpolator ${PCMS_FIELD_TRANSFER_SOURCES})
set_target_properties(pcms_interpolator PROPERTIES
        OUTPUT_NAME pcmsinterpolator
        EXPORT_NAME interpolator)
add_library(pcms::interpolator ALIAS pcms_interpolator)


target_link_libraries(pcms_interpolator PUBLIC pcms::core Kokkos::kokkoskernels)

target_include_directories(pcms_interpolator INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


install(TARGETS pcms_interpolator
    EXPORT interpolatorTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}) 


install(EXPORT interpolatorTargets
    FILE pcms_interpolator-targets.cmake
    NAMESPACE pcms::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

