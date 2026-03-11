set(
        PCMS_UTILITY_HEADERS
        arrays.h
        assert.h
        common.h
        memory_spaces.h
        types.h
        array_mask.h
        inclusive_scan.h
        profile.h
        print.h
)

set(
        PCMS_UTILITY_SOURCES
        assert.cpp
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
target_compile_features(pcms_utility PUBLIC cxx_std_17)

set_target_properties(
        pcms_utility PROPERTIES OUTPUT_NAME pcmsutility EXPORT_NAME
        utility
)

if (PCMS_ENABLE_SPDLOG)
    target_link_libraries(pcms_utility PUBLIC spdlog::spdlog)
endif ()

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