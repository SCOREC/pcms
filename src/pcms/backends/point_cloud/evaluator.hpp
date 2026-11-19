#ifndef PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP
#define PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP

#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/backends/point_cloud/mls_evaluator.hpp"
#include "pcms/geometry/geometry.hpp"

#include <memory>

namespace pcms::point_cloud
{

/// @brief Point evaluation of point-cloud fields by moving least squares.
class PointCloudEvaluatorFactory
  : public pcms::detail::TypedPointEvaluatorFactory<Real,
                                                    pcms::detail::MLSBinding>
{
  using Bound = pcms::detail::MLSBinding;

public:
  /// @param geometry geometry whose Bind computes the MLS supports of query
  ///        points among the layout's DOF holders
  /// @param layout point-cloud layout of the evaluated fields
  /// @param options MLS options
  PointCloudEvaluatorFactory(std::shared_ptr<const Geometry> geometry,
                             std::shared_ptr<const FieldLayout> layout,
                             MLSOptions options = {})
    : pcms::detail::TypedPointEvaluatorFactory<Real, Bound>(
        std::move(geometry)),
      layout_(std::move(layout)),
      options_(options)
  {
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<pcms::detail::MLSPointEvaluator<>>(
      std::move(bound), layout_->GetDimension(), options_, layout_);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  MLSOptions options_;
};

} // namespace pcms::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_EVALUATOR_HPP
