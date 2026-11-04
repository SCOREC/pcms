#ifndef PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
#define PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H

#include "pcms/field/field_evaluator_factory.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/evaluation_request.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/backends/point_cloud/mls_evaluator.hpp"
#include "pcms/field/localization_path_selection.hpp"
#include "pcms/localization/mesh_localization.hpp"
#include "pcms/localization/localization_factory.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"

#include <Omega_h_array.hpp>
#include <memory>
#include <string>

namespace pcms::backend::point_cloud
{

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
    // MLS distances are based on the unweighted Euclidean measure
    // this is only valid with identity metric
    // this almost certainly should be generalized to take metric for distances
    // for now, we just check and fail so at least the user gets a warning
    if (!HasIdentityMetric(*layout_->GetCoordinateSystem())) {
      throw pcms_error(
        "PointCloudEvaluatorFactory: MLS point-cloud evaluation requires a "
        "coordinate system whose metric is the identity; got '" +
        std::string(layout_->GetCoordinateSystem()->Kind()) + "'");
    }

    // Extract source coordinates for the MLS solve.
    const auto src_view = layout_->GetDOFHolderCoordinates().GetValues();
    const int dim = layout_->GetDimension();
    Omega_h::Reals source_coords =
      pcms::utility::flatten_to_omega_h_reals(src_view, "src_coords");

    // Extract target coordinates for the MLS solve.
    const auto tgt_view = coords.GetValues();
    PCMS_ALWAYS_ASSERT(static_cast<int>(tgt_view.extent(1)) == dim);
    Omega_h::Reals target_coords_oh =
      pcms::utility::flatten_to_omega_h_reals(tgt_view, "tgt_coords");

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

} // namespace pcms::backend::point_cloud

#endif // PCMS_POINT_CLOUD_EVALUATOR_FACTORY_H
