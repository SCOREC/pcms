# Point-cloud backend (no external dependency). No public headers.

# layout and field data: part of pcms::field
target_sources(pcms_field PRIVATE layout.cpp field_data.cpp                                   field_data.hpp)

# MLS evaluation and the geometry
target_sources(pcms_field PRIVATE evaluator.hpp )

target_sources(
  pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/point_cloud/geometry.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/point_cloud/layout.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/point_cloud/mls_evaluator.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/field/evaluator/mls_interpolation.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/localization/mls_support_helpers.hpp)
