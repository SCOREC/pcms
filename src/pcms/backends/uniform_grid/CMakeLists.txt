# Uniform-grid backend (no external dependency). No public headers.

# layout, evaluators and spline interpolators: part of pcms::field
target_sources(
  pcms_field
  PRIVATE layout.cpp layout.hpp evaluator.hpp spline_interpolator.hpp
          spline_evaluator.hpp linear_interpolant.hpp
          multidimarray.hpp)

target_sources(
  pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/uniform_grid/geometry.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/localization/uniform_grid_search.hpp)
