#ifndef WDM_COUPLING_COORDINATE_H
#define WDM_COUPLING_COORDINATE_H
#include <array>
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
namespace wdmcpl
{
template <typename CoordinateSystem, typename CoordT = Real, auto N = 3>
class Coordinate
{
public:
  using ScalarT = CoordT;
  template <typename... Args>
  constexpr Coordinate(Args... args) : data_{std::forward<Args>(args)...}
  {}
  [[nodiscard]]
  constexpr nonstd::span<const ScalarT, N> Values() const noexcept { return data_; }
  [[nodiscard]]
  constexpr ScalarT operator[](std::size_t i) const { return data_[i]; }

private:
  std::array<CoordT, N> data_;
};
} // namespace wdmcpl

#endif // WDM_COUPLING_COORDINATE_H
