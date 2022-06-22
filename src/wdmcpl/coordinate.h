#ifndef WDM_COUPLING_COORDINATE_H
#define WDM_COUPLING_COORDINATE_H
#include <array>
#include <span>
#include "wdmcpl/types.h"
namespace wdmcpl
{
template <typename CoordinateSystem, typename CoordT, auto N = 3>
class Coordinate
{
public:
  using ScalarT = CoordT;
  Coordinate() = default;
  constexpr Coordinate(const std::array<ScalarT,N> & data) : data_{data} {}
  constexpr Coordinate(std::array<ScalarT,N>&& data) : data_{std::move(data)} {}

  Coordinate(const Coordinate&) = default;
  Coordinate(Coordinate&&) noexcept = default;
  Coordinate& operator=(const Coordinate&) = default;
  Coordinate& operator=(Coordinate&&) noexcept = default;
  ~Coordinate() = default;
  [[nodiscard]]
  constexpr std::span<const ScalarT, N> Values() const noexcept { return data_; }
  [[nodiscard]]
  constexpr ScalarT operator[](std::size_t i) const { return data_[i]; }

private:
  std::array<ScalarT, N> data_;
};
template<typename T=Real>
struct ScalarData {
  T data;
};
} // namespace wdmcpl

#endif // WDM_COUPLING_COORDINATE_H
