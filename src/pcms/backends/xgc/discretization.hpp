#ifndef PCMS_BACKENDS_XGC_DISCRETIZATION_HPP
#define PCMS_BACKENDS_XGC_DISCRETIZATION_HPP

#include "pcms/discretization/discretization.hpp"
#include "pcms/backends/xgc/reverse_classification.hpp"

#include <Kokkos_Core.hpp>

namespace pcms::backend::xgc
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

} // namespace pcms::backend::xgc

#endif // PCMS_BACKENDS_XGC_DISCRETIZATION_HPP
