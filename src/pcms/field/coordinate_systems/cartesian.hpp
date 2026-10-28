#ifndef PCMS_FIELD_COORDINATE_SYSTEMS_CARTESIAN_HPP
#define PCMS_FIELD_COORDINATE_SYSTEMS_CARTESIAN_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/utility/assert.h"
#include <string>

namespace pcms
{
namespace csys
{

class Cartesian final : public CoordinateSystem
{
public:
  [[nodiscard]] static std::shared_ptr<const CoordinateSystem> Create(int dim)
  {
    if (dim < 1 || dim > 3) {
      throw pcms_error("csys::Cartesian::Create: dimension must be 1..3, got " +
                       std::to_string(dim));
    }
    return std::shared_ptr<const CoordinateSystem>(new Cartesian(dim));
  }

  [[nodiscard]] static std::shared_ptr<const CoordinateSystem> Deferred()
  {
    return std::shared_ptr<const CoordinateSystem>(new Cartesian(0));
  }

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
  explicit Cartesian(int dim) noexcept : dim_(dim) {}

  int dim_;
};

} // namespace csys
} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_SYSTEMS_CARTESIAN_HPP
