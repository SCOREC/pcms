#ifndef PCMS_FIELD_EVALUATOR_POINT_CLOUD_GEOMETRY_HPP
#define PCMS_FIELD_EVALUATOR_POINT_CLOUD_GEOMETRY_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/evaluator/mls_point_cloud.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/geometry/geometry.hpp"
#include "pcms/localization/mls_support_helpers.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/entity_types.h"
#include "pcms/utility/omega_h_array_utils.h"
#include "pcms/utility/profile.h"

#include <Omega_h_array.hpp>
#include <memory>
#include <string>

namespace pcms
{

/// Geometry of a set of sites (the DOF holders of a layout: a point cloud or
/// the entities of one dimension of a mesh): the identity map, located by the
/// MLS supports of each query point. Bind finds the supports with
/// Neighbourhood, using the options the geometry was created with.
///
/// The supports are searched by the geometry the site layout is on
/// (FieldLayout::GetGeometry) when it has one and its DOF holders are the
/// entities of one dimension of it without ghosts, e.g. the adjacency search of
/// an Omega_h mesh for its vertices or element centroids; otherwise, for a bare
/// point cloud, by distance. It holds only the site layout, and the layout
/// never holds a point-cloud geometry built on it.
/// @throws pcms_error on construction if the sites' coordinate system does not
///         have the identity metric
class PointCloudGeometry final : public Geometry
{
public:
  /// Own it by std::shared_ptr, as CreatePointCloudGeometry does: Bind uses
  /// shared_from_this().
  PointCloudGeometry(std::shared_ptr<const FieldLayout> site_layout,
                     NeighbourhoodOptions options)
    : site_layout_(std::move(site_layout)), options_(std::move(options))
  {
    // MLS distances are based on the unweighted Euclidean measure, which is
    // valid only with the identity metric; this almost certainly should be
    // generalized to take the metric for distances
    if (!HasIdentityMetric(*site_layout_->GetCoordinateSystem())) {
      throw pcms_error(
        "PointCloudGeometry: MLS point-cloud evaluation requires a "
        "coordinate system whose metric is the identity; got '" +
        std::string(site_layout_->GetCoordinateSystem()->Kind()) + "'");
    }
  }

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return site_layout_->GetCoordinateSystem();
  }

  [[nodiscard]] const std::shared_ptr<const FieldLayout>& GetSiteLayout()
    const noexcept
  {
    return site_layout_;
  }

  /// A hint naming the cloud or the geometry of its site layout, whose
  /// entities the sites are, is used (and forwarded to that geometry).
  [[nodiscard]] bool TakesHints(const OwnerHints& owners,
                                LO num_points) const noexcept override
  {
    const auto& site_geometry = site_layout_->GetGeometry();
    return owners.Names(*this, num_points) ||
           (site_geometry != nullptr &&
            owners.Names(*site_geometry, num_points));
  }

  /// The supports among the sites: the search of the site layout's geometry
  /// among the entities the sites are, when there is one (see the class
  /// comment), with `owners` forwarded to it (verified there), and the
  /// distance search otherwise. `site_dim` is ignored: the sites are the DOF
  /// holders of the site layout.
  [[nodiscard]] SupportResults Neighbourhood(
    CoordinateView<DeviceMemorySpace> coords, int /*site_dim*/,
    const NeighbourhoodOptions& options,
    const OwnerHints& owners = {}) const override
  {
    detail::CheckNeighbourhoodArguments(*this, coords, options,
                                        "PointCloudGeometry::Neighbourhood");
    const auto& site_geometry = site_layout_->GetGeometry();
    const int site_dim = site_layout_->GetDOFHolderEntityDim();
    // the supports are entity ids, which are site indices only without ghosts
    if (site_geometry != nullptr && site_dim >= 0 &&
        !site_layout_->HasGhosts()) {
      const auto n = static_cast<LO>(coords.GetValues().extent(0));
      // a hint naming the cloud names the entities of its sites
      const OwnerHints forwarded =
        owners.Names(*this, n)
          ? OwnerHints{site_geometry, owners.entity_dim, owners.entities}
          : owners;
      return site_geometry->Neighbourhood(std::move(coords), site_dim, options,
                                          forwarded);
    }
    const auto sites = site_layout_->GetDOFHolderCoordinates().GetValues();
    return BuildPointCloudSupports(
      flatten_to_omega_h_reals(sites, "src_coords"),
      flatten_to_omega_h_reals(coords.GetValues(), "tgt_coords"),
      static_cast<int>(sites.extent(1)), options.radius, options.min_supports,
      options.adapt_radius);
  }

protected:
  std::shared_ptr<BoundSource> BindImpl(
    const BindRequest& request) const override
  {
    PCMS_FUNCTION_TIMER;
    // Extract source coordinates for the MLS solve.
    const auto src_view = site_layout_->GetDOFHolderCoordinates().GetValues();
    Omega_h::Reals source_coords =
      flatten_to_omega_h_reals(src_view, "src_coords");

    const int site_dim = site_layout_->GetDOFHolderEntityDim();
    SupportResults supports =
      Neighbourhood(request.coords, site_dim, options_, request.owners);

    return std::make_shared<MLSBinding>(
      shared_from_this(), request, std::move(source_coords),
      std::move(supports), site_dim, options_);
  }

private:
  std::shared_ptr<const FieldLayout> site_layout_;
  NeighbourhoodOptions options_;
};

/// A point-cloud geometry whose sites are `coords` (copied to the device).
/// @param coords the sites
/// @param options the support search of the geometry
/// @throws pcms_error if the coordinate system does not have the identity
///         metric (MLS distances are Euclidean)
inline std::shared_ptr<const PointCloudGeometry> CreatePointCloudGeometry(
  CoordinateView<HostMemorySpace> coords, NeighbourhoodOptions options = {})
{
  const auto values = coords.GetValues();
  const int dim = static_cast<int>(values.extent(1));
  Kokkos::View<const Real**, Kokkos::HostSpace> host_view(
    values.data_handle(), values.extent(0), values.extent(1));
  auto device_view = Kokkos::View<Real**>("device_view", host_view.extent(0),
                                          host_view.extent(1));
  DeepCopyMismatchLayouts(device_view, host_view);
  auto layout = std::make_shared<PointCloudLayout>(
    dim, device_view, coords.GetCoordinateSystem());
  return std::make_shared<PointCloudGeometry>(std::move(layout),
                                              std::move(options));
}

/// A point-cloud geometry whose sites are the DOF holders of `sites` (e.g. the
/// entities of one dimension of a mesh, through OmegaHEntityLayout, whose
/// geometry's adjacency search then serves).
/// @param sites the layout of the sites
/// @param options the support search of the geometry
/// @throws pcms_error if `sites` is null or its coordinate system does not
///         have the identity metric
inline std::shared_ptr<const PointCloudGeometry> CreatePointCloudGeometry(
  std::shared_ptr<const FieldLayout> sites, NeighbourhoodOptions options = {})
{
  if (sites == nullptr) {
    throw pcms_error("CreatePointCloudGeometry: the sites must not be null");
  }
  return std::make_shared<PointCloudGeometry>(std::move(sites),
                                              std::move(options));
}

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_POINT_CLOUD_GEOMETRY_HPP
