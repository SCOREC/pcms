#ifndef PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
#define PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H

#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/evaluator/mls_point_cloud.h"
#include "pcms/localization/localization_path_selection.h"
#include "pcms/localization/mesh_localization.h"
#include "pcms/localization/localization_factory.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/omega_h_array_utils.h"

#include <Omega_h_array.hpp>
#include <memory>

namespace pcms
{

// PointCloudEvaluatorFactory implements FieldEvaluatorFactory<Real> for
// reconstructed fields that provide source coordinates through FieldLayout.
//
// Support localization is delegated to a LocalizationFactory, allowing
// different search backends (N² point-cloud or mesh-adjacency BFS) to be
// plugged in without changing this class. CreatePointEvaluator calls
// LocalizationFactory::Build once for the supplied target coordinates and
// returns an MLSPointEvaluator that can be reused at zero additional
// localization cost.
class PointCloudEvaluatorFactory : public FieldEvaluatorFactory<Real>
{
public:
  PointCloudEvaluatorFactory(std::shared_ptr<const FieldLayout> layout,
                             std::shared_ptr<LocalizationFactory> localization,
                             MLSOptions options = {})
    : layout_(std::move(layout)),
      localization_(std::move(localization)),
      options_(options)
  {
  }

  const FieldLayout& GetLayout() const override { return *layout_; }

  CoordinateSystem GetCoordinateSystem() const override
  {
    return layout_->GetDOFHolderCoordinates().GetCoordinateSystem();
  }

  bool HasDOFHolderCoordinates() const override { return true; }

  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return layout_->GetDOFHolderCoordinates();
  }

  bool SupportsNearestBoundary() const override { return false; }

  std::unique_ptr<PointEvaluator<Real>> CreatePointEvaluator(
    const EvaluationRequest& request) const override
  {
    const auto coords = request.coords;
    if (coords.GetCoordinateSystem() != GetCoordinateSystem()) {
      throw pcms_error(
        "PointCloudEvaluatorFactory: coordinate system mismatch");
    }
    if (GetCoordinateSystem() != CoordinateSystem::Cartesian) {
      throw pcms_error(
        "PointCloudEvaluatorFactory: only Cartesian coordinates are "
        "supported for MLS point-cloud evaluation");
    }

    // Extract source coordinates for the MLS solve.
    const auto src_view = layout_->GetDOFHolderCoordinates().GetValues();
    const int dim = layout_->GetDimension();
    Omega_h::Reals source_coords =
      flatten_to_omega_h_reals(src_view, "src_coords");

    // Extract target coordinates for the MLS solve.
    const auto tgt_view = coords.GetValues();
    PCMS_ALWAYS_ASSERT(static_cast<int>(tgt_view.extent(1)) == dim);
    Omega_h::Reals target_coords_oh =
      flatten_to_omega_h_reals(tgt_view, "tgt_coords");

    SupportResults supports;
    auto path =
      detail::SelectLocalizationPath(*layout_, request.GetQueryLayout());
    if (path == detail::LocalizationPath::CentroidToVertexAdjacencySearch) {
      auto* adjacency =
        dynamic_cast<const AdjacencyLocalizationFactory*>(localization_.get());
      PCMS_ALWAYS_ASSERT(adjacency != nullptr);
      supports = adjacency->BuildSameMeshCentroidToVertex();
    } else {
      supports = localization_->Build(coords);
    }

    return std::make_unique<MLSPointEvaluator<>>(
      std::move(source_coords), std::move(target_coords_oh),
      std::move(supports), dim, options_);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::shared_ptr<LocalizationFactory> localization_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
