#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_H

#include "pcms/utility/arrays.hpp"
#include "pcms/utility/entity_types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <memory>

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

/// Discretization of a uniform (structured) grid: its vertices and cells.
std::shared_ptr<const Discretization> CreateUniformGridDiscretization(
  const UniformGrid<2>& grid);
std::shared_ptr<const Discretization> CreateUniformGridDiscretization(
  const UniformGrid<3>& grid);

/// Discretization with no entities, for spaces that only evaluate (e.g.
/// fields provided by an external code).
std::shared_ptr<const Discretization> CreateEmptyDiscretization();

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_H
