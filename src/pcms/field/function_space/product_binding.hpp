#ifndef PCMS_FIELD_FUNCTION_SPACE_PRODUCT_BINDING_HPP
#define PCMS_FIELD_FUNCTION_SPACE_PRODUCT_BINDING_HPP

// Private: the product geometry class and its binding. (Compiled in pcms_field,
// not pcms_geometry: the twist-shift identification evaluates a Function.)

#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/product_geometry.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms
{

class ProductGeometry;

/// BoundPoints of a product geometry: the factors' bound points of the factor
/// parts of the query points.
class ProductBinding final : public BoundPoints
{
public:
  ProductBinding(std::shared_ptr<const Geometry> geometry,
                 const BindRequest& request);

  [[nodiscard]] const std::shared_ptr<const BoundPoints>& FactorA()
    const noexcept
  {
    return a_;
  }
  [[nodiscard]] const std::shared_ptr<const BoundPoints>& FactorB()
    const noexcept
  {
    return b_;
  }

protected:
  void RebindImpl() override;

  /// The block composition of the factors' geometric factors: the
  /// block-diagonal Jacobian and metric, and the product of the volume
  /// factors.
  [[nodiscard]] detail::GeometricFactors ComputeFactors() const override;

private:
  void Locate(bool rebind);

  std::shared_ptr<const BoundPoints> a_;
  std::shared_ptr<const BoundPoints> b_;
};

class ProductGeometry final : public Geometry
{
public:
  ProductGeometry(std::shared_ptr<const Geometry> a,
                  std::shared_ptr<const Geometry> b,
                  std::optional<ProductIdentification> identification);

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return system_;
  }
  [[nodiscard]] bool SupportsNearestBoundary() const noexcept override
  {
    return a_->SupportsNearestBoundary() && b_->SupportsNearestBoundary();
  }
  [[nodiscard]] const std::shared_ptr<const Geometry>& FactorA() const noexcept
  {
    return a_;
  }
  [[nodiscard]] const std::shared_ptr<const Geometry>& FactorB() const noexcept
  {
    return b_;
  }

  [[nodiscard]] const std::optional<ProductIdentification>& Identification()
    const noexcept
  {
    return identification_;
  }

  /// The tensor product of a rule of each factor: the point (i, k) is the
  /// i-th point of `qa` with the k-th of `qb`, at index i * qb points + k,
  /// with the product of the weights.
  [[nodiscard]] Quadrature TensorQuadrature(const Quadrature& qa,
                                            const Quadrature& qb) const;

protected:
  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

  /// The tensor product of the factors' rules.
  [[nodiscard]] Quadrature CreateQuadrature(
    const detail::QuadratureRule& rule) const override;

private:
  std::shared_ptr<const Geometry> a_;
  std::shared_ptr<const Geometry> b_;
  std::shared_ptr<const CoordinateSystem> system_;
  std::optional<ProductIdentification> identification_;
};

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACE_PRODUCT_BINDING_HPP
