#ifndef PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP
#define PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP

#include "pcms/backends/point_cloud/mls_evaluator.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/field/point_rows.hpp"

#include <memory>

namespace pcms::point_cloud
{

/// @brief Point evaluation of point-cloud fields by moving least squares.
/// @param geometry geometry whose Bind computes the MLS supports of query
///        points among the layout's DOF holders
/// @param layout point-cloud layout of the evaluated fields
inline std::shared_ptr<const PointEvaluatorFactory> CreateMLSEvaluatorFactory(
  std::shared_ptr<const Geometry> geometry,
  std::shared_ptr<const FieldLayout> layout, MLSOptions options)
{
  using Bound = pcms::detail::MLSBinding;
  const int dim = layout->GetDimension();
  return std::make_shared<pcms::detail::RowsEvaluatorFactory<Bound>>(
    std::move(geometry), std::move(layout), std::set<DifferentialOperator>{},
    [dim, options](const Bound& bound, DifferentialOperator /*op*/) {
      const auto& supports = bound.Supports();
      return pcms::detail::PointRows{
        supports.supports_ptr, supports.supports_idx,
        pcms::detail::mls_coefficients(
          bound.SourceCoordinates(), bound.TargetCoordinates(), supports, dim,
          static_cast<LO>(options.degree), options.basis, options.lambda,
          options.tol, options.decay_factor)};
    },
    true);
}

} // namespace pcms::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP
