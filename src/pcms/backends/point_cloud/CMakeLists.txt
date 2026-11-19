# Point-cloud backend (no external dependency). Public header:
# include/pcms/backends/point_cloud/geometry.hpp.

# layout, field data, the MLS kernel and evaluator, and the geometry: part of
# pcms::field
target_sources(
  pcms_field
  PRIVATE layout.cpp field_data.cpp geometry.cpp mls_interpolation.cpp
          field_data.hpp layout.hpp evaluator.hpp mls_evaluator.hpp
          mls_interpolation.hpp mls_interpolation_impl.hpp
          interpolator_aliases.hpp interpolator_view_utils.hpp
          interpolator_logger.hpp point_cloud_geometry.hpp
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
         ${PROJECT_SOURCE_DIR}/include/pcms/backends/point_cloud/geometry.hpp)
