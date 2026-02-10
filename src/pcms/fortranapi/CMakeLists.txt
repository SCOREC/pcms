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

add_library(pcms_fortranapi_interpolator pcms_interpolator.f90 interpolator_wrap.c pcms_mesh.f90 mesh_wrap.c)
add_library(pcms::fortranapi::interpolator ALIAS pcms_fortranapi_interpolator)
set_target_properties(pcms_fortranapi_interpolator PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran
        OUTPUT_NAME pcmsfortranapiinterpolator
        EXPORT_NAME fortranapi::interpolator)
target_include_directories(pcms_fortranapi_interpolator PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:pcms::fortranapi::interpolator,Fortran_MODULE_DIRECTORY>>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran>)
target_link_libraries(pcms_fortranapi_interpolator PRIVATE pcms::capi::interpolator PUBLIC Kokkos::kokkos)

INSTALL(DIRECTORY $<TARGET_PROPERTY:pcms::fortranapi::core,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(DIRECTORY $<TARGET_PROPERTY:pcms::fortranapi::interpolator,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

# high level interface target
add_library(pcms_fortranapi INTERFACE)
add_library(pcms::fortranapi ALIAS pcms_fortranapi)
set_target_properties(pcms_fortranapi PROPERTIES EXPORT_NAME fortranapi)
# link capi libraries to a high level interface library
target_link_libraries(pcms_fortranapi INTERFACE pcms::fortranapi::core)
target_link_libraries(pcms_fortranapi INTERFACE pcms::fortranapi::interpolator)


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
  TARGETS pcms_fortranapi_interpolator
  EXPORT pcms_fortranapi_interpolator-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/fortranapi/)

install(
  EXPORT pcms_fortranapi_interpolator-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
