#ifndef PCMS_BACKENDS_POINT_CLOUD_DISCRETIZATION_HPP
#define PCMS_BACKENDS_POINT_CLOUD_DISCRETIZATION_HPP

#include "pcms/discretization/discretization.hpp"

#include <Kokkos_Core.hpp>

namespace pcms::backend::point_cloud
{

class PointCloudDiscretization : public Discretization
{
public:
  PointCloudDiscretization(int dim,
                           Kokkos::View<const Real**, HostMemorySpace> coords,
                           const void* entity_identity = nullptr);

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;

private:
  int dim_;
  const void* entity_identity_;
  Kokkos::View<ClassificationDimension*, DeviceMemorySpace> class_dims_;
  Kokkos::View<ClassificationId*, DeviceMemorySpace> class_ids_;
};

} // namespace pcms::backend::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_DISCRETIZATION_HPP
