#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/entity_layout.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::omega_h
{

std::shared_ptr<const OmegaHGeometry> CreateGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options)
{
  return std::make_shared<OmegaHGeometry>(mesh, std::move(system),
                                          std::move(options));
}

std::shared_ptr<const OmegaHGeometry> CreateGeometry(
  std::shared_ptr<Omega_h::Mesh> mesh,
  std::shared_ptr<Omega_h::Library> library,
  std::shared_ptr<const CoordinateSystem> system, OmegaHGeometryOptions options)
{
  return std::make_shared<OmegaHGeometry>(
    std::move(mesh), std::move(library), std::move(system), std::move(options));
}

std::shared_ptr<const point_cloud::PointCloudGeometry> CreatePointCloudGeometry(
  std::shared_ptr<const OmegaHGeometry> geometry, int entity_dim,
  NeighbourhoodOptions options)
{
  if (geometry == nullptr) {
    throw pcms_error("pcms::omega_h::CreatePointCloudGeometry: the "
                     "geometry must not be null");
  }
  if (entity_dim != 0 && entity_dim != geometry->GetMesh().dim()) {
    throw pcms_error("pcms::omega_h::CreatePointCloudGeometry: "
                     "entity_dim must be 0 (vertices) or the mesh dimension "
                     "(element centroids)");
  }
  return point_cloud::CreateGeometry(
    std::make_shared<OmegaHEntityLayout>(std::move(geometry), entity_dim, 1),
    std::move(options));
}

} // namespace pcms::omega_h
