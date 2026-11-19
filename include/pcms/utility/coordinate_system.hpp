#ifndef PCMS_UTILITY_COORDINATE_SYSTEM_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEM_HPP

#include <memory>
#include <string_view>

namespace pcms
{

class CoordinateSystem
{
public:
  [[nodiscard]] virtual std::string_view Kind() const noexcept = 0;

  [[nodiscard]] virtual int Dimension() const noexcept = 0;

  [[nodiscard]] virtual bool HasOrthogonalBasis() const noexcept = 0;

  [[nodiscard]] virtual bool HasUnitScaleFactors() const noexcept = 0;

  [[nodiscard]] virtual bool operator==(
    const CoordinateSystem& other) const noexcept = 0;

  virtual ~CoordinateSystem() = default;
};

} // namespace pcms

#include "pcms/utility/detail/coordinate_system.hpp"

#endif // PCMS_UTILITY_COORDINATE_SYSTEM_HPP
