set(
        PCMS_UTILITY_HEADERS
        ../../../include/pcms/utility/arrays.hpp
        ../../../include/pcms/utility/coordinate_system.hpp
        ../../../include/pcms/utility/coordinate_view.hpp
        ../../../include/pcms/utility/coordinate_systems/cartesian.hpp
        ../../../include/pcms/utility/coordinate_systems/cylindrical.hpp
        ../../../include/pcms/utility/assert.hpp
        ../../../include/pcms/utility/bounding_box.hpp
        ../../../include/pcms/utility/detail/common.hpp
        ../../../include/pcms/utility/entity_types.hpp
        ../../../include/pcms/utility/eqdsk.hpp
        mesh_geometry.hpp
        ../../../include/pcms/utility/memory_spaces.hpp
        mpi_type.hpp
        omega_h_array_utils.hpp
        ../../../include/pcms/utility/types.hpp
        array_mask.hpp
        ../../../include/pcms/utility/detail/inclusive_scan.hpp
        ../../../include/pcms/utility/detail/profile.hpp
        ../../../include/pcms/utility/print.hpp
        ../../../include/pcms/utility/uniform_grid.hpp
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
        # public headers (include/pcms/...)
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
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
# Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(PCMS_UTILITY_HEADERS_PUBLIC ${PCMS_UTILITY_HEADERS})
list(FILTER PCMS_UTILITY_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(PCMS_UTILITY_HEADERS_PRIVATE ${PCMS_UTILITY_HEADERS})
list(FILTER PCMS_UTILITY_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_utility
  PUBLIC FILE_SET utilities TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${PCMS_UTILITY_HEADERS_PUBLIC}
  PRIVATE ${PCMS_UTILITY_HEADERS_PRIVATE})
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
        DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
