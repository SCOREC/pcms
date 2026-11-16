#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/mesh_geometry.h"

namespace pcms
{

static_assert(std::is_same_v<ClassificationId, Omega_h::ClassId>,
              "ClassificationId must match Omega_h::ClassId for direct views");
static_assert(
  std::is_same_v<ClassificationDimension, Omega_h::I8>,
  "ClassificationDimension must match Omega_h::I8 for direct views");

OmegaHDiscretization::OmegaHDiscretization(Omega_h::Mesh& mesh) : mesh_(mesh)
{
  for (int entity_dim = 0; entity_dim <= max_entity_dim_; ++entity_dim) {
    class_dims_by_dim_[entity_dim] = Omega_h::Read<Omega_h::I8>();
    class_ids_by_dim_[entity_dim] = Omega_h::Read<Omega_h::ClassId>();
  }

  for (int entity_dim = 0; entity_dim <= mesh_.dim(); ++entity_dim) {
    class_dims_by_dim_[entity_dim] = Omega_h::Read<Omega_h::I8>(
      mesh_.get_array<Omega_h::I8>(entity_dim, "class_dim"));
    class_ids_by_dim_[entity_dim] = Omega_h::Read<Omega_h::ClassId>(
      mesh_.get_array<Omega_h::ClassId>(entity_dim, "class_id"));
  }
}

bool OmegaHDiscretization::SameEntities(
  const Discretization& other) const noexcept
{
  auto p = dynamic_cast<const OmegaHDiscretization*>(&other);
  return p != nullptr && &p->mesh_ == &mesh_;
}

int OmegaHDiscretization::GetDimension() const
{
  return mesh_.dim();
}

LO OmegaHDiscretization::GetNumEntities(int entity_dim) const
{
  if (entity_dim < 0 || entity_dim > mesh_.dim())
    return 0;
  return mesh_.nents(entity_dim);
}

Rank1View<const ClassificationDimension, DeviceMemorySpace>
OmegaHDiscretization::GetEntityClassificationDimensions(int entity_dim) const
{
  if (entity_dim < 0 || entity_dim > max_entity_dim_)
    return Rank1View<const ClassificationDimension, DeviceMemorySpace>(nullptr,
                                                                       0);
  const auto& dims = class_dims_by_dim_[entity_dim];
  return Rank1View<const ClassificationDimension, DeviceMemorySpace>(
    dims.data(), dims.size());
}

Rank1View<const ClassificationId, DeviceMemorySpace>
OmegaHDiscretization::GetEntityClassificationIds(int entity_dim) const
{
  if (entity_dim < 0 || entity_dim > max_entity_dim_)
    return Rank1View<const ClassificationId, DeviceMemorySpace>(nullptr, 0);
  const auto& ids = class_ids_by_dim_[entity_dim];
  return Rank1View<const ClassificationId, DeviceMemorySpace>(ids.data(),
                                                              ids.size());
}

} // namespace pcms
