#ifndef PCMS_FIELD_COORDINATE_SYSTEM_HPP
#define PCMS_FIELD_COORDINATE_SYSTEM_HPP

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


[[nodiscard]] inline bool SameCoordinateSystem(
  const std::shared_ptr<const CoordinateSystem>& a,
  const std::shared_ptr<const CoordinateSystem>& b) noexcept
{
  if (a == nullptr || b == nullptr) {
    return a == b;
  }
  return *a == *b;
}

[[nodiscard]] inline bool HasIdentityMetric(
  const CoordinateSystem& system) noexcept
{
  return system.HasOrthogonalBasis() && system.HasUnitScaleFactors();
}

[[nodiscard]] std::shared_ptr<const CoordinateSystem> ResolveCoordinateSystem(
  std::shared_ptr<const CoordinateSystem> system, int dim);

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_SYSTEM_HPP
