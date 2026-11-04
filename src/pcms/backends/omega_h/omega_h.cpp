#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/entity_layout.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::backend::omega_h
{

std::shared_ptr<const Geometry> CreateGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options)
{
  return pcms::CreateOmegaHGeometry(mesh, std::move(system),
                                    std::move(options));
}

Omega_h::Mesh& GetMesh(const Geometry& geometry)
{
  const auto* omega_h = dynamic_cast<const OmegaHGeometry*>(&geometry);
  if (omega_h == nullptr) {
    throw pcms_error("pcms::backend::omega_h::GetMesh: the geometry is not an "
                     "Omega_h geometry");
  }
  return omega_h->GetMesh();
}

std::shared_ptr<const Geometry> CreatePointCloudGeometry(
  std::shared_ptr<const Geometry> geometry, int entity_dim,
  NeighbourhoodOptions options)
{
  const auto omega_h =
    std::dynamic_pointer_cast<const OmegaHGeometry>(std::move(geometry));
  if (omega_h == nullptr) {
    throw pcms_error("pcms::backend::omega_h::CreatePointCloudGeometry: the "
                     "geometry must be an Omega_h geometry");
  }
  if (entity_dim != 0 && entity_dim != omega_h->GetMesh().dim()) {
    throw pcms_error("pcms::backend::omega_h::CreatePointCloudGeometry: "
                     "entity_dim must be 0 (vertices) or the mesh dimension "
                     "(element centroids)");
  }
  return pcms::CreatePointCloudGeometry(
    std::make_shared<OmegaHEntityLayout>(omega_h, entity_dim, 1),
    std::move(options));
}

} // namespace pcms::backend::omega_h
