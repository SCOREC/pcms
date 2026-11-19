#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/entity_layout.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::omega_h
{

std::shared_ptr<const Geometry> CreateGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options)
{
  return std::make_shared<OmegaHGeometry>(mesh, std::move(system),
                                          std::move(options));
}

std::shared_ptr<const Geometry> CreateGeometry(
  std::shared_ptr<Omega_h::Mesh> mesh,
  std::shared_ptr<Omega_h::Library> library,
  std::shared_ptr<const CoordinateSystem> system, OmegaHGeometryOptions options)
{
  return std::make_shared<OmegaHGeometry>(
    std::move(mesh), std::move(library), std::move(system), std::move(options));
}

std::shared_ptr<const Geometry> CreatePointCloudGeometry(
  std::shared_ptr<const Geometry> geometry, int entity_dim,
  NeighbourhoodOptions options)
{
  if (geometry == nullptr) {
    throw pcms_error("pcms::omega_h::CreatePointCloudGeometry: the "
                     "geometry must not be null");
  }
  auto mesh_geometry =
    std::dynamic_pointer_cast<const OmegaHGeometry>(std::move(geometry));
  if (mesh_geometry == nullptr) {
    throw pcms_error("pcms::omega_h::CreatePointCloudGeometry: the "
                     "geometry must be an Omega_h geometry");
  }
  if (entity_dim != 0 && entity_dim != mesh_geometry->GetMesh().dim()) {
    throw pcms_error("pcms::omega_h::CreatePointCloudGeometry: "
                     "entity_dim must be 0 (vertices) or the mesh dimension "
                     "(element centroids)");
  }
  return point_cloud::CreateGeometry(std::make_shared<OmegaHEntityLayout>(
                                       std::move(mesh_geometry), entity_dim, 1),
                                     std::move(options));
}

} // namespace pcms::omega_h
