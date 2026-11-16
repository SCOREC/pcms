set(
        PCMS_DISCRETIZATION_HEADERS
        discretization.h
        discretization/xgc_reverse_classification.h
        discretization/empty.hpp
        discretization/point_cloud.hpp
        discretization/uniform_grid.hpp
        discretization/xgc.hpp
)
if (PCMS_ENABLE_OMEGA_H)
  list(APPEND PCMS_DISCRETIZATION_HEADERS discretization/omega_h.hpp)
endif()

set(
        PCMS_DISCRETIZATION_SOURCES
        discretization/empty.cpp
        discretization/point_cloud.cpp
        discretization/xgc_reverse_classification.cpp
        discretization/xgc.cpp
)
if (PCMS_ENABLE_OMEGA_H)
  list(APPEND PCMS_DISCRETIZATION_SOURCES discretization/omega_h.cpp)
endif()

add_library(pcms_discretization ${PCMS_DISCRETIZATION_SOURCES})
add_library(pcms::discretization ALIAS pcms_discretization)
target_sources(pcms_discretization PUBLIC
        FILE_SET field_discretization
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_DISCRETIZATION_HEADERS})
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
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms
)
