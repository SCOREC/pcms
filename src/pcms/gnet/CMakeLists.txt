find_package(KokkosKernels REQUIRED)

set(
        PCMS_GNET_HEADERS
        gnet_grid.hpp
        gnet_io.hpp
        gnet_diagnostics.hpp
        gnet_sampling.hpp
        gnet_to_gtc.hpp
)

set(
        PCMS_GNET_SOURCES
        gnet_grid.cpp
        gnet_io.cpp
        gnet_diagnostics.cpp
        gnet_sampling.cpp
        gnet_to_gtc.cpp
)
add_library(pcms_gnet ${PCMS_GNET_SOURCES})
add_library(pcms::gnet ALIAS pcms_gnet)
target_include_directories(
        pcms_gnet
        PUBLIC
        # include path should be pcms/gnet/
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        # include path to search for pcms/config.h
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(
        pcms_gnet
        PUBLIC
        Kokkos::kokkos
        Omega_h::omega_h
        pcms::utility
        pcms::localization
        pcms::field
        PRIVATE
        Kokkos::kokkoskernels
        HDF5::HDF5
)
target_compile_features(pcms_gnet PUBLIC cxx_std_20)

set_target_properties(
        pcms_gnet PROPERTIES OUTPUT_NAME pcmsgnet EXPORT_NAME
        gnet
        # keep the path to HDF5 in the installed library
        INSTALL_RPATH_USE_LINK_PATH TRUE
)

## export the library
set_target_properties(pcms_gnet PROPERTIES PUBLIC_HEADER "${PCMS_GNET_HEADERS}")

install(
        TARGETS pcms_gnet
        EXPORT pcms_gnet-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/gnet
)

install(
        EXPORT pcms_gnet-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
