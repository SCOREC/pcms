# Uniform-grid backend (no external dependency). No public headers.

# layout, evaluators and spline interpolators: part of pcms::field, built
# with PCMS_ENABLE_OMEGA_H as before
if(PCMS_ENABLE_OMEGA_H)
  target_sources(
    pcms_field
    PRIVATE layout.cpp layout.hpp evaluator.hpp spline_interpolator.hpp
            spline_evaluator.hpp geometry.hpp linear_interpolant.hpp
            multidimarray.hpp)
endif()
