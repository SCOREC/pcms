#ifndef PCMS_BACKENDS_POINT_CLOUD_GEOMETRY_HPP
#define PCMS_BACKENDS_POINT_CLOUD_GEOMETRY_HPP

#include "pcms/field/field_layout.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <memory>

namespace pcms::point_cloud
{

/// A point-cloud geometry whose sites are `coords` (copied to the device).
/// @param coords the sites
/// @param options the support search of the geometry
/// @throws pcms_error if the coordinate system does not have the identity
///         metric (MLS distances are Euclidean)
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  CoordinateView<HostMemorySpace> coords, NeighbourhoodOptions options = {});

/// A point-cloud geometry whose sites are the DOF holders of `sites` (e.g. the
/// entities of one dimension of a mesh, through OmegaHEntityLayout, whose
/// geometry's adjacency search then serves).
/// @param sites the layout of the sites
/// @param options the support search of the geometry
/// @throws pcms_error if `sites` is null or its coordinate system does not
///         have the identity metric
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  std::shared_ptr<const FieldLayout> sites, NeighbourhoodOptions options = {});

} // namespace pcms::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_GEOMETRY_HPP
