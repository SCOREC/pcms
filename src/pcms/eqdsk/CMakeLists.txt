# EQDSK equilibria as PCMS functions, on the uniform-grid backend (no external
# dependency): part of pcms::field. Public header: include/pcms/eqdsk/.
target_sources(
  pcms_field
  PRIVATE eqdsk.cpp
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
         ${PROJECT_SOURCE_DIR}/include/pcms/eqdsk/eqdsk.hpp)
