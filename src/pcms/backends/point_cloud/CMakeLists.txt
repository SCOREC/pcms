# Point-cloud backend. Public header:
# include/pcms/backends/point_cloud/geometry.hpp.

# layout, field data and the geometry with its support search: part of
# pcms::field
target_sources(
  pcms_field
  PRIVATE layout.cpp field_data.cpp geometry.cpp field_data.hpp layout.hpp
          mls_evaluator.hpp mls_interpolation.hpp point_cloud_geometry.hpp
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
         ${PROJECT_SOURCE_DIR}/include/pcms/backends/point_cloud/geometry.hpp)

# the polynomial-reconstruction space on point clouds: its MLS kernel and
# evaluator
if(PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION)
  target_sources(
    pcms_field
    PRIVATE ../../field/function_space/polynomial_reconstruction.cpp
            mls_interpolation.cpp evaluator.hpp mls_interpolation_impl.hpp
            interpolator_aliases.hpp interpolator_view_utils.hpp
            interpolator_logger.hpp)
  # the MLS solve uses KokkosBatched
  target_link_libraries(pcms_field PRIVATE Kokkos::kokkoskernels)
endif()
