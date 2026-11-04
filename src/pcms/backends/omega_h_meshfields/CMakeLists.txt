# Omega_h + MeshFields bridge: the MeshFields field adapter, part of
# pcms::field. Built only when both backends are enabled.
target_sources(pcms_field PRIVATE layout.cpp)
target_sources(
  pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../..
         FILES layout.hpp evaluator_backend.hpp evaluator.hpp field_data.hpp)
target_link_libraries(pcms_field PUBLIC meshfields::meshfields)
