add_library(wdmcpl_fortranapi wdmcpl.f90 client_wrap.c)
target_compile_definitions(wdmcpl_fortranapi PUBLIC HAVE_MPI)
add_library(wdmcpl::fortranapi ALIAS wdmcpl_fortranapi)
set_target_properties(wdmcpl_fortranapi PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran
        OUTPUT_NAME wdmcplfortranapi
        EXPORT_NAME fortranapi)
target_include_directories(wdmcpl_fortranapi PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:wdmcpl::fortranapi,Fortran_MODULE_DIRECTORY>>
                                                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran>)
target_link_libraries(wdmcpl_fortranapi PRIVATE wdmcpl::capi)
INSTALL(DIRECTORY $<TARGET_PROPERTY:wdmcpl::fortranapi,Fortran_MODULE_DIRECTORY>
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
        TARGETS wdmcpl_fortranapi
        EXPORT wdmcpl_fortranapi-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wdmcpl/fortranapi/)
install(
        EXPORT wdmcpl_fortranapi-targets
        NAMESPACE wdmcpl::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)
