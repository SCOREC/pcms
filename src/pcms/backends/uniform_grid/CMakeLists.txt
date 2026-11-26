# Uniform-grid backend (no external dependency). Public header:
# include/pcms/backends/uniform_grid/geometry.hpp.

# geometry, search, layout, evaluators and spline interpolators: part of
# pcms::field
target_sources(
  pcms_field
  PRIVATE geometry.cpp layout.cpp layout.hpp evaluator.hpp
          spline_interpolator.hpp spline_evaluator.hpp uniform_grid_geometry.hpp
          uniform_grid_search.hpp
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
         ${PROJECT_SOURCE_DIR}/include/pcms/backends/uniform_grid/geometry.hpp)
