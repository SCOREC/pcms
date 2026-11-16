set(
        PCMS_UTILITY_HEADERS
        arrays.h
        assert.h
        bounding_box.h
        common.h
        entity_types.h
        eqdsk.h
        mesh_geometry.h
        memory_spaces.h
        mpi_type.h
        omega_h_array_utils.h
        types.h
        array_mask.h
        inclusive_scan.h
        profile.h
        print.h
        uniform_grid.h
)

set(
        PCMS_UTILITY_SOURCES
        assert.cpp
        eqdsk.cpp
        print.cpp
)
add_library(pcms_utility ${PCMS_UTILITY_SOURCES})
add_library(pcms::utility ALIAS pcms_utility)
target_include_directories(
        pcms_utility
        PUBLIC
        # include path should be pcms/utility/
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        # include path to search for pcms/config.h
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(pcms_utility PUBLIC Kokkos::kokkos perfstubs)
if(PCMS_ENABLE_OMEGA_H)
    target_link_libraries(pcms_utility PUBLIC Omega_h::omega_h)
endif()
target_compile_features(pcms_utility PUBLIC cxx_std_20)

set_target_properties(
        pcms_utility PROPERTIES OUTPUT_NAME pcmsutility EXPORT_NAME
        utility
)

## export the library
target_sources(pcms_utility PUBLIC
        FILE_SET utilities
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../
        FILES ${PCMS_UTILITY_HEADERS})
install(
        TARGETS pcms_utility
        EXPORT pcms_utility-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILE_SET utilities DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

install(
        EXPORT pcms_utility-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
