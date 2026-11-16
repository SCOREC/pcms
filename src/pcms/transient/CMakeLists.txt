set(PCMS_TRANSIENT_HEADERS timestepper.hpp)
set(PCMS_TRANSIENT_SOURCES timestepper.cpp)

add_library(pcms_transient ${PCMS_TRANSIENT_SOURCES})
set_target_properties(pcms_transient PROPERTIES EXPORT_NAME transient)
add_library(pcms::transient ALIAS pcms_transient)

target_sources(
  pcms_transient
  PUBLIC FILE_SET transient
         TYPE HEADERS
         BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
         FILES ${PCMS_TRANSIENT_HEADERS})
target_compile_features(pcms_transient PUBLIC cxx_std_20)
target_link_libraries(pcms_transient PUBLIC pcms::utility)
target_include_directories(
  pcms_transient
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
  TARGETS pcms_transient
  EXPORT pcms_transient-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/transient
  FILE_SET transient DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

install(
  EXPORT pcms_transient-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
