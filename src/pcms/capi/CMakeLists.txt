add_library(pcms_capi client.cpp kokkos.cpp)
add_library(pcms::capi ALIAS pcms_capi)
target_link_libraries(pcms_capi PUBLIC MPI::MPI_C PRIVATE pcms::core)

target_include_directories(pcms_capi
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes our include path pcms/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/capi>")

set(CAPI_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/client.h ${CMAKE_CURRENT_SOURCE_DIR}/kokkos.h)
set_target_properties(pcms_capi PROPERTIES PUBLIC_HEADER "${CAPI_HEADERS}"
         OUTPUT_NAME pcmscapi
         EXPORT_NAME capi)

add_library(pcms_interpolator_capi interpolator.cpp kokkos.cpp)
add_library(pcms::interpolator::capi ALIAS pcms_interpolator_capi)
target_link_libraries(pcms_interpolator_capi PUBLIC MPI::MPI_C PRIVATE pcms::interpolator)
target_include_directories( pcms_interpolator_capi
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes the module path cpms/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/capi>")
set(INT_CAPI_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.h)
set_target_properties(pcms_interpolator_capi PROPERTIES PUBLIC_HEADERS "${INT_CAPI_HEADERS}"
        OUTPUTNAME pcmsinterpolatorcapi
        EXPORT_NAME interpolatorcapi)

install(
        TARGETS pcms_capi pcms_interpolator_capi
        EXPORT pcms_capi-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/capi/)
install(
        EXPORT pcms_capi-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

