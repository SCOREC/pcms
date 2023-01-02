add_library(capi client.cpp)
add_library(wdmcpl::capi ALIAS capi)
target_link_libraries(capi PUBLIC wdmcpl::wdmcpl)

target_include_directories(capi
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>" # this makes our include path wdmcpl/capi
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/capi>")

set_target_properties(capi PROPERTIES
        PUBLIC_HEADER "client.h")
install(
        TARGETS capi
        EXPORT wdmcpl-capi-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/capi
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#configure_package_config_file(
#        "${CMAKE_SOURCE_DIR}/config.cmake.in"
#        "${CMAKE_BINARY_DIR}/wdmcpl-config.cmake"
#        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl
#)
#write_basic_package_version_file(
#        "${PROJECT_BINARY_DIR}/wdmcpl-config-version.cmake"
#        COMPATIBILITY AnyNewerVersion)
#
#install(FILES
#        "${PROJECT_BINARY_DIR}/wdmcpl-config.cmake"
#        "${PROJECT_BINARY_DIR}/wdmcpl-config-version.cmake"
#        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)

install(
        EXPORT wdmcpl-capi-targets
        NAMESPACE wdmcpl::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wdmcpl)
