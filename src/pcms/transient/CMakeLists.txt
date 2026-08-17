set(PCMS_TRANSIENT_HEADERS
    participant.hpp
    accelerator.hpp
    timestepper.hpp)

add_library(pcms_transient INTERFACE)
set_target_properties(pcms_transient PROPERTIES EXPORT_NAME transient)
add_library(pcms::transient ALIAS pcms_transient)

target_sources(
  pcms_transient
  INTERFACE FILE_SET transient
            TYPE HEADERS
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
            FILES ${PCMS_TRANSIENT_HEADERS})
target_compile_features(pcms_transient INTERFACE cxx_std_20)
target_link_libraries(pcms_transient INTERFACE pcms::utility)
target_include_directories(
  pcms_transient
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
  TARGETS pcms_transient
  EXPORT pcms_transient-targets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/transient
  FILE_SET transient DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms)

install(
  EXPORT pcms_transient-targets
  NAMESPACE pcms::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
