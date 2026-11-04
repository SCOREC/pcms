# Uniform-grid backend (no external dependency).

# discretization: part of pcms::discretization (header only)
target_sources(
  pcms_discretization
  PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES discretization.hpp)

# layout, evaluators and spline interpolators: part of pcms::field, built
# with PCMS_ENABLE_OMEGA_H as before
if(PCMS_ENABLE_OMEGA_H)
  target_sources(pcms_field PRIVATE layout.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES layout.hpp evaluator.hpp
           spline_interpolator.hpp spline_evaluator.hpp
           linear_interpolant.hpp multidimarray.hpp)
endif()
