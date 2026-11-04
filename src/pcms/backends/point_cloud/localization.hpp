#ifndef PCMS_BACKENDS_POINT_CLOUD_LOCALIZATION_HPP
#define PCMS_BACKENDS_POINT_CLOUD_LOCALIZATION_HPP

#include "pcms/localization/localization_factory.hpp"
#include "pcms/localization/mls_support_helpers.hpp"
#include "pcms/backends/point_cloud/layout.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"

#include <Omega_h_array.hpp>
#include <memory>
#include <string>

namespace pcms::backend::point_cloud
{

// LocalizationFactory implementation using N^2 point-cloud search.
//
// Holds a reference to the PointCloudLayout so source coordinates are not
// copied at construction time; they are extracted from the layout in Build().
class PointCloudLocalizationFactory : public pcms::localization::LocalizationFactory
{
public:
  PointCloudLocalizationFactory(std::shared_ptr<const PointCloudLayout> layout,
                                MLSOptions options)
    : layout_(std::move(layout)), options_(options)
  {
  }

  pcms::localization::SupportResults Build(
    CoordinateView<DeviceMemorySpace> target_coords) const override
  {
    // Euclidean search structures are valid wherever
    // the metric is the identity.
    if (!HasIdentityMetric(*target_coords.GetCoordinateSystem())) {
      throw pcms_error(
        "PointCloudLocalizationFactory: requires a coordinate system whose "
        "metric is the identity; "
        "got '" +
        std::string(target_coords.GetCoordinateSystem()->Kind()) + "'");
    }

    const auto src_view = layout_->GetCoordinates();
    const int dim = layout_->GetDimension();
    Omega_h::Reals source_coords =
      pcms::utility::flatten_to_omega_h_reals(src_view, "src_coords");

    const auto tgt_view = target_coords.GetValues();
    PCMS_ALWAYS_ASSERT(static_cast<int>(tgt_view.extent(1)) == dim);
    Omega_h::Reals target_coords_oh =
      pcms::utility::flatten_to_omega_h_reals(tgt_view, "tgt_coords");

    return pcms::localization::BuildPointCloudSupports(source_coords, target_coords_oh, dim,
                                   options_.radius, options_.min_req_supports,
                                   options_.adapt_radius);
  }

private:
  std::shared_ptr<const PointCloudLayout> layout_;
  MLSOptions options_;
};

} // namespace pcms::backend::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_LOCALIZATION_HPP
