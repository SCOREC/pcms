#ifndef PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
#define PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H

#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/evaluator/mls_point_cloud.h"
#include "pcms/geometry/geometry.hpp"

#include <memory>

namespace pcms
{

/// @brief Point evaluation of point-cloud fields by moving least squares.
class PointCloudEvaluatorFactory
  : public TypedPointEvaluatorFactory<Real, MLSBinding>
{
  using Bound = MLSBinding;

public:
  /// @param geometry geometry whose Bind computes the MLS supports of query
  ///        points among the layout's DOF holders
  /// @param layout point-cloud layout of the evaluated fields
  /// @param options MLS options
  PointCloudEvaluatorFactory(std::shared_ptr<const Geometry> geometry,
                             std::shared_ptr<const FieldLayout> layout,
                             MLSOptions options = {})
    : TypedPointEvaluatorFactory<Real, Bound>(std::move(geometry)),
      layout_(std::move(layout)),
      options_(options)
  {
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<MLSPointEvaluator<>>(
      std::move(bound), layout_->GetDimension(), options_, layout_);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
