add_library(pcms_capi_core client.cpp kokkos.cpp)
add_library(pcms::capi::core ALIAS pcms_capi_core)

target_link_libraries(
  pcms_capi_core
  PUBLIC MPI::MPI_C
  PRIVATE pcms::core)

target_include_directories(
  pcms_capi_core
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes our
                                                                # include path
                                                                # pcms/capi
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/capi>")

set(CAPI_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/client.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/kokkos.h)
set_target_properties(
  pcms_capi_core PROPERTIES OUTPUT_NAME pcmscapicore
             EXPORT_NAME capi::core)
target_sources(pcms_capi_core PUBLIC
        FILE_SET core
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${CAPI_CORE_HEADERS})
add_library(pcms_capi_transfer mesh.cpp interpolator.cpp kokkos.cpp)
add_library(pcms::capi::transfer ALIAS pcms_capi_transfer)
target_link_libraries(pcms_capi_transfer PUBLIC MPI::MPI_C PRIVATE pcms::transfer)
target_include_directories(pcms_capi_transfer
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes the module path pcms/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pcms/capi>")

set(CAPI_TRANSFER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/mesh.h)

set_target_properties(pcms_capi_transfer
  PROPERTIES OUTPUT_NAME pcmscapitransfer
             EXPORT_NAME capi::transfer)
target_sources(pcms_capi_transfer PUBLIC
        FILE_SET transfer
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
	FILES ${CAPI_TRANSFER_HEADERS})
# high level interface target
add_library(pcms_capi INTERFACE)
add_library(pcms::capi ALIAS pcms_capi)
set_target_properties(pcms_capi PROPERTIES EXPORT_NAME capi)
# link capi libraries to a high level interface library
target_link_libraries(pcms_capi INTERFACE pcms::capi::core)
target_link_libraries(pcms_capi INTERFACE pcms::capi::transfer)
install(
  TARGETS pcms_capi
  EXPORT pcms_capi-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)
install(
  EXPORT pcms_capi-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

# library targets
install(
  TARGETS pcms_capi_core
  EXPORT pcms_capi_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILE_SET core DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)
install(
  EXPORT pcms_capi_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  TARGETS pcms_capi_transfer
  EXPORT pcms_capi_transfer-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILE_SET transfer DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)
install(
  EXPORT pcms_capi_transfer-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
