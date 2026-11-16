#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_POINT_CLOUD_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_POINT_CLOUD_H

#include "pcms/discretization/discretization.h"

#include <Kokkos_Core.hpp>

namespace pcms
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

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_POINT_CLOUD_H
