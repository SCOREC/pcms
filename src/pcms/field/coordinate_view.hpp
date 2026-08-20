#ifndef PCMS_FIELD_COORDINATE_VIEW_HPP
#define PCMS_FIELD_COORDINATE_VIEW_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <memory>
#include <string>

namespace pcms
{

template <typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
class CoordinateView
{
public:
  CoordinateView(std::shared_ptr<const CoordinateSystem> coordinate_system,
                 Rank2View<const Real, MemorySpace, LayoutPolicy> coords)
    : system_(ResolveCoordinateSystem(std::move(coordinate_system),
                                      static_cast<int>(coords.extent(1)))),
      coordinates_(coords)
  {
  }

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept
  {
    return system_;
  }

  [[nodiscard]] Rank2View<const Real, MemorySpace, LayoutPolicy> GetValues()
    const noexcept
  {
    return coordinates_;
  }

private:
  std::shared_ptr<const CoordinateSystem> system_;
  Rank2View<const Real, MemorySpace, LayoutPolicy> coordinates_;
};

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_VIEW_HPP
