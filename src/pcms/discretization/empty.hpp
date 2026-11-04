#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_EMPTY_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_EMPTY_H

#include "pcms/discretization/discretization.hpp"

#include <Kokkos_Core.hpp>

namespace pcms::discretization
{

class EmptyDiscretization : public Discretization
{
public:
  EmptyDiscretization();

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;
};

} // namespace pcms::discretization

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_EMPTY_H
