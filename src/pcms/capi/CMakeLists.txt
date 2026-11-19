# the C API client wraps the coupler with XGC fields
if(PCMS_BUILD_CAPI_CLIENT)
  add_library(pcms_capi_core client.cpp kokkos.cpp)
  add_library(pcms::capi::core ALIAS pcms_capi_core)

  target_link_libraries(
    pcms_capi_core
    PUBLIC MPI::MPI_C
    PRIVATE pcms::coupler Omega_h::omega_h)

  target_include_directories(
    pcms_capi_core
    PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes our
                                                                  # include path
                                                                  # pcms/capi
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

  set(CAPI_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../../include/pcms/capi/client.h
                   ${CMAKE_CURRENT_SOURCE_DIR}/../../../include/pcms/capi/kokkos.h)
  set_target_properties(
    pcms_capi_core PROPERTIES OUTPUT_NAME pcmscapicore
               EXPORT_NAME capi::core)
  # Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(CAPI_CORE_HEADERS_PUBLIC ${CAPI_CORE_HEADERS})
list(FILTER CAPI_CORE_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(CAPI_CORE_HEADERS_PRIVATE ${CAPI_CORE_HEADERS})
list(FILTER CAPI_CORE_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_capi_core
  PUBLIC FILE_SET core TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${CAPI_CORE_HEADERS_PUBLIC}
  PRIVATE ${CAPI_CORE_HEADERS_PRIVATE})
endif()

# the C transfer API interpolates on Omega_h meshes by polynomial
# reconstruction
if(PCMS_BUILD_CAPI_TRANSFER)
add_library(pcms_capi_transfer mesh.cpp interpolator.cpp kokkos.cpp)
add_library(pcms::capi::transfer ALIAS pcms_capi_transfer)
target_link_libraries(pcms_capi_transfer PUBLIC MPI::MPI_C PRIVATE pcms::transfer
                                                         Omega_h::omega_h)
target_include_directories(pcms_capi_transfer
        PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes the module path pcms/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

set(CAPI_TRANSFER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../../include/pcms/capi/interpolator.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/../../../include/pcms/capi/kokkos.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/../../../include/pcms/capi/mesh.h)

set_target_properties(pcms_capi_transfer
  PROPERTIES OUTPUT_NAME pcmscapitransfer
             EXPORT_NAME capi::transfer)
# Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(CAPI_TRANSFER_HEADERS_PUBLIC ${CAPI_TRANSFER_HEADERS})
list(FILTER CAPI_TRANSFER_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(CAPI_TRANSFER_HEADERS_PRIVATE ${CAPI_TRANSFER_HEADERS})
list(FILTER CAPI_TRANSFER_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_capi_transfer
  PUBLIC FILE_SET transfer TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${CAPI_TRANSFER_HEADERS_PUBLIC}
  PRIVATE ${CAPI_TRANSFER_HEADERS_PRIVATE})
endif()

# the headers are C: the header-set verification compiles them as C
set_source_files_properties(${CAPI_CORE_HEADERS} ${CAPI_TRANSFER_HEADERS}
                            PROPERTIES LANGUAGE C)

# high level interface target
add_library(pcms_capi INTERFACE)
add_library(pcms::capi ALIAS pcms_capi)
set_target_properties(pcms_capi PROPERTIES EXPORT_NAME capi)
# link capi libraries to a high level interface library
if(TARGET pcms_capi_core)
  target_link_libraries(pcms_capi INTERFACE pcms::capi::core)
endif()
if(TARGET pcms_capi_transfer)
  target_link_libraries(pcms_capi INTERFACE pcms::capi::transfer)
endif()
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
  DESTINATION ${PCMS_INSTALL_CMAKEDIR})

# library targets
if(TARGET pcms_capi_core)
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
    DESTINATION ${PCMS_INSTALL_CMAKEDIR})
endif()

if(TARGET pcms_capi_transfer)
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
  DESTINATION ${PCMS_INSTALL_CMAKEDIR})
endif()
