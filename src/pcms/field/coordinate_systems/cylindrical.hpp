#ifndef PCMS_FIELD_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
#define PCMS_FIELD_COORDINATE_SYSTEMS_CYLINDRICAL_HPP

#include "pcms/field/coordinate_system.hpp"

namespace pcms
{
namespace csys
{

class CylindricalRThetaZ final : public CoordinateSystem
{
public:
  [[nodiscard]] static std::shared_ptr<const CoordinateSystem> Create()
  {
    return std::shared_ptr<const CoordinateSystem>(new CylindricalRThetaZ());
  }

  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cylindrical_rthetaz";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 3; }
  [[nodiscard]] bool HasOrthogonalBasis() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool HasUnitScaleFactors() const noexcept override
  {
    return false;
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    return dynamic_cast<const CylindricalRThetaZ*>(&other) != nullptr;
  }

private:
  CylindricalRThetaZ() noexcept = default;
};

class CylindricalRZ final : public CoordinateSystem
{
public:
  [[nodiscard]] static std::shared_ptr<const CoordinateSystem> Create()
  {
    return std::shared_ptr<const CoordinateSystem>(new CylindricalRZ());
  }

  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cylindrical_rz";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 2; }
  [[nodiscard]] bool HasOrthogonalBasis() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool HasUnitScaleFactors() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    return dynamic_cast<const CylindricalRZ*>(&other) != nullptr;
  }

private:
  CylindricalRZ() noexcept = default;
};

} // namespace csys
} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
