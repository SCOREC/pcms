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
    interpolation_base.cpp
)

add_library(pcms_interpolator ${PCMS_FIELD_TRANSFER_SOURCES})
set_target_properties(pcms_interpolator PROPERTIES
        OUTPUT_NAME pcmsinterpolator
        EXPORT_NAME interpolator)
target_sources(pcms_interpolator PUBLIC
        FILE_SET transformer
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_FIELD_TRANSFER_HEADERS})
add_library(pcms::interpolator ALIAS pcms_interpolator)
target_compile_features(pcms_interpolator PUBLIC cxx_std_17)

target_link_libraries(pcms_interpolator PUBLIC pcms::core PRIVATE Kokkos::kokkoskernels)

target_include_directories(pcms_interpolator INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


install(
  TARGETS pcms_interpolator
  EXPORT pcms_interpolator-targets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator
  FILE_SET transformer DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


install(
    EXPORT pcms_interpolator-targets
    NAMESPACE pcms::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
