#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/assert.hpp"
#include <string>

namespace pcms
{

namespace
{

class Cartesian final : public CoordinateSystem
{
public:
  explicit Cartesian(int dim) noexcept : dim_(dim) {}

  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return dim_ == 0 ? "cartesian (dimension deferred)" : "cartesian";
  }
  [[nodiscard]] int Dimension() const noexcept override { return dim_; }
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
    const auto* p = dynamic_cast<const Cartesian*>(&other);
    return p != nullptr && p->dim_ == dim_;
  }

private:
  int dim_;
};

class CylindricalRThetaZ final : public CoordinateSystem
{
public:
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
};

class CylindricalRZ final : public CoordinateSystem
{
public:
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
};

} // namespace

namespace csys
{

std::shared_ptr<const CoordinateSystem> CreateCartesian(int dim)
{
  if (dim < 1 || dim > 3) {
    throw pcms_error("csys::CreateCartesian: dimension must be 1..3, got " +
                     std::to_string(dim));
  }
  return std::make_shared<const Cartesian>(dim);
}

std::shared_ptr<const CoordinateSystem> CreateCylindricalRThetaZ()
{
  return std::make_shared<const CylindricalRThetaZ>();
}

std::shared_ptr<const CoordinateSystem> CreateCylindricalRZ()
{
  return std::make_shared<const CylindricalRZ>();
}

} // namespace csys

namespace detail
{

std::shared_ptr<const CoordinateSystem> DeferredCartesian()
{
  return std::make_shared<const Cartesian>(0);
}

bool IsCartesian(const CoordinateSystem& system) noexcept
{
  return dynamic_cast<const Cartesian*>(&system) != nullptr;
}

std::shared_ptr<const CoordinateSystem> ResolveCoordinateSystem(
  std::shared_ptr<const CoordinateSystem> system, int dim)
{
  if (system == nullptr) {
    throw pcms_error(
      "ResolveCoordinateSystem: coordinate system must not be null");
  }
  // To avoid requiring the user specifying the dimension twice in
  // FunctionSpace constructors, we have a default
  // Cartesian constructor that constructs with dimension 0.
  if (IsCartesian(*system) && system->Dimension() == 0) {
    return csys::CreateCartesian(dim);
  }
  if (system->Dimension() != dim) {
    throw pcms_error("ResolveCoordinateSystem: coordinate system '" +
                     std::string(system->Kind()) + "' has " +
                     std::to_string(system->Dimension()) +
                     " coordinate columns but the data has " +
                     std::to_string(dim));
  }
  return system;
}

} // namespace detail
} // namespace pcms
