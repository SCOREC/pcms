#ifndef WDM_COUPLING_COORDINATE_H
#define WDM_COUPLING_COORDINATE_H
#include <array>
#include <span>
#include "wdmcpl/types.h"
#include <Kokkos_Macros.hpp>
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
  constexpr const std::array<CoordT,N>& Values() const noexcept { return data_; }
  [[nodiscard]]
  constexpr ScalarT operator[](std::size_t i) const { return data_[i]; }

private:
  std::array<CoordT, N> data_;
};

// TODO add explicit conversion to underlying type
template <typename CoordinateSystem, typename ElementType = Real>
struct CoordinateElement
{
  using value_type = std::remove_cvref_t<ElementType>;
  using reference = value_type&;
  using const_reference = const value_type &;
  using coordinate_system = CoordinateSystem;
  CoordinateElement(value_type data) : data_(std::move(data)) {}

  [[nodiscard]]
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr const_reference underlying() const noexcept{ return data_ ; }
  [[nodiscard]]
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr reference underlying() noexcept { return data_ ; }

private:
  value_type data_;
};

} // namespace wdmcpl

#endif // WDM_COUPLING_COORDINATE_H
