add_library(pcms_fortranapi_core pcms.f90 client_wrap.c)
target_compile_definitions(pcms_fortranapi_core PUBLIC HAVE_MPI)
add_library(pcms::fortranapi::core ALIAS pcms_fortranapi_core)
set_target_properties(
  pcms_fortranapi_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran
             OUTPUT_NAME pcmsfortranapi EXPORT_NAME fortranapi::core)
target_include_directories(
  pcms_fortranapi_core
  PUBLIC
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:pcms::fortranapi,Fortran_MODULE_DIRECTORY>>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran>)
target_link_libraries(
  pcms_fortranapi_core
  PRIVATE pcms::capi::core
  PUBLIC Kokkos::kokkos)

add_library(pcms_fortranapi_transfer pcms_interpolator.f90 interpolator_wrap.c pcms_mesh.f90 mesh_wrap.c)
add_library(pcms::fortranapi::transfer ALIAS pcms_fortranapi_transfer)
set_target_properties(pcms_fortranapi_transfer PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran
        OUTPUT_NAME pcmsfortranapi_transfer
        EXPORT_NAME fortranapi::transfer)
target_include_directories(pcms_fortranapi_transfer PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:pcms::fortranapi::transfer,Fortran_MODULE_DIRECTORY>>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran>)
target_link_libraries(pcms_fortranapi_transfer PRIVATE pcms::capi::transfer PUBLIC Kokkos::kokkos)

INSTALL(DIRECTORY $<TARGET_PROPERTY:pcms::fortranapi::core,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(DIRECTORY $<TARGET_PROPERTY:pcms::fortranapi::transfer,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

# high level interface target
add_library(pcms_fortranapi INTERFACE)
add_library(pcms::fortranapi ALIAS pcms_fortranapi)
set_target_properties(pcms_fortranapi PROPERTIES EXPORT_NAME fortranapi)
# link capi libraries to a high level interface library
target_link_libraries(pcms_fortranapi INTERFACE pcms::fortranapi::core)
target_link_libraries(pcms_fortranapi INTERFACE pcms::fortranapi::transfer)


install(
  TARGETS pcms_fortranapi
  EXPORT pcms_fortranapi-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/fortranapi/)

install(
  EXPORT pcms_fortranapi-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  TARGETS pcms_fortranapi_core
  EXPORT pcms_fortranapi_core-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/fortranapi/)

install(
  EXPORT pcms_fortranapi_core-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)

install(
  TARGETS pcms_fortranapi_transfer
  EXPORT pcms_fortranapi_transfer-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/fortranapi/)

install(
  EXPORT pcms_fortranapi_transfer-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
