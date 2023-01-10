add_library(wdmcpl_capi client.cpp)
add_library(wdmcpl::capi ALIAS wdmcpl_capi)
target_link_libraries(wdmcpl_capi PUBLIC wdmcpl::core)

target_include_directories(wdmcpl_capi
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes our include path wdmcpl/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/wdmcpl/capi>")

set_target_properties(wdmcpl_capi PROPERTIES PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/client.h
         OUTPUT_NAME wdmcplcapi
         EXPORT_NAME capi)
install(
        TARGETS wdmcpl_capi
        EXPORT wdmcpl_capi-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wdmcpl/capi/)
install(
        EXPORT wdmcpl_capi-targets
        NAMESPACE wdmcpl::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)
