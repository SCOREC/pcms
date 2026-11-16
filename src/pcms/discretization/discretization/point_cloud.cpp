#include "pcms/discretization/discretization/point_cloud.hpp"

namespace pcms
{

PointCloudDiscretization::PointCloudDiscretization(
  int dim, Kokkos::View<const Real**, HostMemorySpace> coords,
  const void* entity_identity)
  : dim_(dim),
    entity_identity_(entity_identity != nullptr
                       ? entity_identity
                       : static_cast<const void*>(coords.data())),
    class_dims_("point_cloud_discretization_class_dims", coords.extent(0)),
    class_ids_("point_cloud_discretization_class_ids", coords.extent(0))
{
  Kokkos::deep_copy(class_dims_, static_cast<ClassificationDimension>(dim_));
  Kokkos::deep_copy(class_ids_, static_cast<ClassificationId>(0));
}

bool PointCloudDiscretization::SameEntities(
  const Discretization& other) const noexcept
{
  auto p = dynamic_cast<const PointCloudDiscretization*>(&other);
  return p != nullptr && p->entity_identity_ == entity_identity_;
}

int PointCloudDiscretization::GetDimension() const
{
  return dim_;
}

LO PointCloudDiscretization::GetNumEntities(int entity_dim) const
{
  return entity_dim == Vertex ? static_cast<LO>(class_dims_.extent(0)) : 0;
}

Rank1View<const ClassificationDimension, DeviceMemorySpace>
PointCloudDiscretization::GetEntityClassificationDimensions(
  int entity_dim) const
{
  return entity_dim == Vertex
           ? make_const_array_view(class_dims_)
           : Rank1View<const ClassificationDimension, DeviceMemorySpace>(
               nullptr, 0);
}

Rank1View<const ClassificationId, DeviceMemorySpace>
PointCloudDiscretization::GetEntityClassificationIds(int entity_dim) const
{
  return entity_dim == Vertex
           ? make_const_array_view(class_ids_)
           : Rank1View<const ClassificationId, DeviceMemorySpace>(nullptr, 0);
}

} // namespace pcms
