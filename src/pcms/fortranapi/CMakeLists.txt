add_library(pcms_fortranapi pcms.f90 client_wrap.c)
target_compile_definitions(pcms_fortranapi PUBLIC HAVE_MPI)
add_library(pcms::fortranapi ALIAS pcms_fortranapi)
set_target_properties(pcms_fortranapi PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran
        OUTPUT_NAME pcmsfortranapi
        EXPORT_NAME fortranapi)
target_include_directories(pcms_fortranapi PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:pcms::fortranapi,Fortran_MODULE_DIRECTORY>>
                                                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran>)
target_link_libraries(pcms_fortranapi PRIVATE pcms::capi PUBLIC Kokkos::kokkos)
INSTALL(DIRECTORY $<TARGET_PROPERTY:pcms::fortranapi,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
        TARGETS pcms_fortranapi
        EXPORT pcms_fortranapi-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/fortranapi/)
install(
        EXPORT pcms_fortranapi-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
