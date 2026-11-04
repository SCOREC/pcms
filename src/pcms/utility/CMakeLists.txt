set(
        PCMS_UTILITY_HEADERS
        arrays.hpp
        coordinate_system.hpp
        coordinate_view.hpp
        coordinate_systems/cartesian.hpp
        coordinate_systems/cylindrical.hpp
        assert.hpp
        bounding_box.hpp
        common.hpp
        entity_types.hpp
        eqdsk.hpp
        mesh_geometry.hpp
        memory_spaces.hpp
        mpi_type.hpp
        omega_h_array_utils.hpp
        types.hpp
        array_mask.hpp
        inclusive_scan.hpp
        profile.hpp
        print.hpp
        uniform_grid.hpp
)

set(
        PCMS_UTILITY_SOURCES
        assert.cpp
        coordinate_system.cpp
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
