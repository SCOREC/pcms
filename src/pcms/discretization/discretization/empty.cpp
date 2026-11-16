#include "pcms/discretization/discretization/empty.hpp"

namespace pcms
{

EmptyDiscretization::EmptyDiscretization() {}

bool EmptyDiscretization::SameEntities(
  const Discretization& other) const noexcept
{
  return dynamic_cast<const EmptyDiscretization*>(&other) != nullptr;
}

int EmptyDiscretization::GetDimension() const
{
  return 0;
}

LO EmptyDiscretization::GetNumEntities(int) const
{
  return 0;
}

Rank1View<const ClassificationDimension, DeviceMemorySpace>
EmptyDiscretization::GetEntityClassificationDimensions(int) const
{
  return {};
}

Rank1View<const ClassificationId, DeviceMemorySpace>
EmptyDiscretization::GetEntityClassificationIds(int) const
{
  return {};
}

} // namespace pcms
