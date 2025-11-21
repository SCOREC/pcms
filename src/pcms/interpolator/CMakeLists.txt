find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_TRANSFER_HEADERS
    pcms_interpolator_aliases.hpp
    adj_search.hpp
    mls_interpolation_impl.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
    mls_interpolation.hpp
    pcms_interpolator_view_utils.hpp
    pcms_interpolator_logger.hpp
    interpolation_base.h
    interpolation_helpers.h
    spline_interpolator.hpp)

set(PCMS_FIELD_TRANSFER_SOURCES
    mls_interpolation.cpp
    interpolation_base.cpp)
install(FILES ${PCMS_FIELD_TRANSFER_HEADERS}
        DESTINATION include/pcms/interpolator)


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


install(
  TARGETS pcms_interpolator
  EXPORT interpolatorTargets
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


install(
    EXPORT interpolatorTargets
    FILE pcms_interpolator-targets.cmake
    NAMESPACE pcms::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
