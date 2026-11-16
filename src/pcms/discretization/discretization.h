#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_H

#include "pcms/utility/arrays.h"
#include "pcms/utility/entity_types.h"

namespace pcms
{

using ClassificationDimension = int8_t;
using ClassificationId = LO;

class Discretization
{
public:
  virtual bool SameEntities(const Discretization& other) const noexcept = 0;

  virtual int GetDimension() const = 0;

  virtual LO GetNumEntities(int entity_dim) const = 0;

  virtual Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const = 0;

  virtual Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const = 0;

  virtual ~Discretization() noexcept = default;
};

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_H
