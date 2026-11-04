set(
        PCMS_DISCRETIZATION_HEADERS
        ../../../include/pcms/discretization/discretization.hpp
        empty.hpp
)

set(
        PCMS_DISCRETIZATION_SOURCES
        empty.cpp
        factories.cpp
)

add_library(pcms_discretization ${PCMS_DISCRETIZATION_SOURCES})
add_library(pcms::discretization ALIAS pcms_discretization)
# Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(PCMS_DISCRETIZATION_HEADERS_PUBLIC ${PCMS_DISCRETIZATION_HEADERS})
list(FILTER PCMS_DISCRETIZATION_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(PCMS_DISCRETIZATION_HEADERS_PRIVATE ${PCMS_DISCRETIZATION_HEADERS})
list(FILTER PCMS_DISCRETIZATION_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_discretization
  PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${PCMS_DISCRETIZATION_HEADERS_PUBLIC}
  PRIVATE ${PCMS_DISCRETIZATION_HEADERS_PRIVATE})
target_include_directories(
        pcms_discretization
        PUBLIC
        # include path should be pcms/discretization/
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>"
        # include path to search for pcms/config.h
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(pcms_discretization PUBLIC Kokkos::kokkos pcms::utility)

if (PCMS_ENABLE_OMEGA_H)
  target_link_libraries(pcms_discretization PRIVATE Omega_h::omega_h)
endif()

target_compile_features(pcms_discretization PUBLIC cxx_std_20)

set_target_properties(
        pcms_discretization PROPERTIES OUTPUT_NAME pcmsdiscretization EXPORT_NAME
        discretization
)

## export the library
install(
        TARGETS pcms_discretization
        EXPORT pcms_discretization-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILE_SET field_discretization DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
)

install(
        EXPORT pcms_discretization-targets
        NAMESPACE pcms::
        DESTINATION ${PCMS_INSTALL_CMAKEDIR}
)
