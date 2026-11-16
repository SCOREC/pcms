#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_XGC_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_XGC_H

#include "pcms/discretization/discretization.h"
#include "pcms/discretization/discretization/xgc_reverse_classification.h"

#include <Kokkos_Core.hpp>

namespace pcms
{

class XGCDiscretization : public Discretization
{
public:
  XGCDiscretization(const ReverseClassificationVertex& reverse_classification,
                    LO num_plane_nodes);

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;

private:
  const ReverseClassificationVertex* reverse_classification_;
  LO num_plane_nodes_;
  Kokkos::View<ClassificationDimension*, DeviceMemorySpace> class_dims_;
  Kokkos::View<ClassificationId*, DeviceMemorySpace> class_ids_;
};

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_XGC_H
