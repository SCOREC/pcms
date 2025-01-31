find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_TRANSFER_HEADERS
    pcms_interpolator_aliases.hpp
    adj_search.hpp
    mls_interpolation_impl.hpp
    queue_visited.hpp
    linear_interpolant.hpp
    multidimarray.hpp
    mls_interpolation.hpp
)

set(PCMS_FIELD_TRANSFER_SOURCES
    mls_interpolation.cpp
)
install(FILES ${PCMS_FIELD_TRANSFER_HEADERS} DESTINATION include/pcms/interpolator)


add_library(pcms_interpolator ${PCMS_FIELD_TRANSFER_SOURCES})
set_target_properties(pcms_interpolator PROPERTIES
        OUTPUT_NAME pcmsinterpolator
        EXPORT_NAME interpolator)
add_library(pcms::interpolator ALIAS pcms_interpolator)


target_link_libraries(pcms_interpolator PUBLIC pcms::core Kokkos::kokkoskernels)


install(TARGETS pcms_interpolator
    EXPORT interpolatorTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/interpolator
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}) 


install(EXPORT interpolatorTargets
    FILE pcms_interpolator-targets.cmake
    NAMESPACE pcms::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

