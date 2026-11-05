# Self-containment check: every installed public header must compile on its
# own (only its own includes, nothing before it). One translation unit is
# generated per header and compiled into an object library; the test
# "public_headers_self_contained" builds that library.

file(GLOB_RECURSE pcms_public_headers RELATIVE "${PROJECT_SOURCE_DIR}/include"
     CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/pcms/*.hpp"
     "${PROJECT_SOURCE_DIR}/include/pcms/*.h")

# headers whose library is not part of this build
set(pcms_skipped_header_regex "")
if(NOT TARGET pcms::coupler)
  list(APPEND pcms_skipped_header_regex "^pcms/coupler/")
endif()
if(NOT TARGET pcms::transient)
  list(APPEND pcms_skipped_header_regex "^pcms/transient/")
endif()
if(NOT PCMS_ENABLE_C)
  list(APPEND pcms_skipped_header_regex "^pcms/capi/")
endif()
if(NOT TARGET pcms::capi::core)
  list(APPEND pcms_skipped_header_regex "^pcms/capi/client\\.h$")
endif()
if(NOT PCMS_ENABLE_OMEGA_H)
  list(APPEND pcms_skipped_header_regex "^pcms/backends/omega_h/"
       "^pcms/field/eqdsk_field\\.hpp$")
endif()
if(NOT PCMS_ENABLE_MFEM)
  list(APPEND pcms_skipped_header_regex "^pcms/backends/mfem/")
endif()
if(NOT PCMS_ENABLE_XGC)
  list(APPEND pcms_skipped_header_regex "^pcms/backends/xgc/xgc\\.hpp$")
endif()

set(pcms_header_check_sources "")
foreach(header IN LISTS pcms_public_headers)
  set(skip FALSE)
  foreach(regex IN LISTS pcms_skipped_header_regex)
    if(header MATCHES "${regex}")
      set(skip TRUE)
    endif()
  endforeach()
  if(skip)
    continue()
  endif()
  string(MAKE_C_IDENTIFIER "${header}" name)
  if(header MATCHES "\\.h$")
    set(source "${CMAKE_CURRENT_BINARY_DIR}/${name}.c")
  else()
    set(source "${CMAKE_CURRENT_BINARY_DIR}/${name}.cpp")
  endif()
  file(CONFIGURE OUTPUT "${source}" CONTENT "#include <${header}>\n")
  list(APPEND pcms_header_check_sources "${source}")
endforeach()

add_library(pcms_public_header_check OBJECT EXCLUDE_FROM_ALL
            ${pcms_header_check_sources})
target_link_libraries(pcms_public_header_check PRIVATE pcms::pcms)
if(TARGET mfem)
  target_link_libraries(pcms_public_header_check PRIVATE mfem)
endif()

add_test(
  NAME public_headers_self_contained
  COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target
          pcms_public_header_check)
