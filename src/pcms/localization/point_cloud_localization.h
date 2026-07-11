#ifndef PCMS_LOCALIZATION_POINT_CLOUD_LOCALIZATION_H
#define PCMS_LOCALIZATION_POINT_CLOUD_LOCALIZATION_H

#include "pcms/localization/localization_factory.h"
#include "pcms/localization/mls_support_helpers.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/omega_h_array_utils.h"

#include <Omega_h_array.hpp>
#include <memory>

namespace pcms
{

// LocalizationFactory implementation using N^2 point-cloud search.
//
// Holds a reference to the PointCloudLayout so source coordinates are not
// copied at construction time; they are extracted from the layout in Build().
class PointCloudLocalizationFactory : public LocalizationFactory
{
public:
  PointCloudLocalizationFactory(std::shared_ptr<const PointCloudLayout> layout,
                                MLSOptions options)
    : layout_(std::move(layout)), options_(options)
  {
  }

  SupportResults Build(
    CoordinateView<DeviceMemorySpace> target_coords) const override
  {
    if (target_coords.GetCoordinateSystem() != CoordinateSystem::Cartesian) {
      throw pcms_error(
        "PointCloudLocalizationFactory: only Cartesian coordinates are "
        "supported");
    }

    const auto src_view = layout_->GetCoordinates();
    const int dim = layout_->GetDimension();
    Omega_h::Reals source_coords =
      flatten_to_omega_h_reals(src_view, "src_coords");

    const auto tgt_view = target_coords.GetValues();
    PCMS_ALWAYS_ASSERT(static_cast<int>(tgt_view.extent(1)) == dim);
    Omega_h::Reals target_coords_oh =
      flatten_to_omega_h_reals(tgt_view, "tgt_coords");

    return BuildPointCloudSupports(source_coords, target_coords_oh, dim,
                                   options_.radius, options_.min_req_supports,
                                   options_.adapt_radius);
  }

private:
  std::shared_ptr<const PointCloudLayout> layout_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_POINT_CLOUD_LOCALIZATION_H
