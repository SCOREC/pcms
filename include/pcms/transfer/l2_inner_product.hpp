#ifndef PCMS_TRANSFER_L2_INNER_PRODUCT_HPP
#define PCMS_TRANSFER_L2_INNER_PRODUCT_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms
{

/// @brief Integrals of the fields of a space by quadrature. Vectors and
///        covectors are paired by the metric of their coordinate system.
class L2Integral
{
public:
  /// @param space the space of u, whose geometry the quadrature is on
  /// @param other the evaluator source of v (a space, or an operator on one),
  ///        with as many components as `space`; null for `space` itself
  /// @param quadrature a quadrature of the space's geometry in place of the
  ///        space's own
  /// @throws pcms_error if no quadrature is given and the space has none, the
  ///         quadrature is not on the space's geometry, or the evaluators
  ///         cannot be created
  explicit L2Integral(std::shared_ptr<const FunctionSpace> space,
                      std::shared_ptr<const PointEvaluatorFactory> other = {},
                      std::optional<Quadrature> quadrature = {});
  ~L2Integral();

  /// @return the integral of a scalar field
  /// @throws pcms_error for a field of several components
  [[nodiscard]] Real Apply(const Field<Real>& u) const;

  /// @return the integral of u . v
  [[nodiscard]] Real Inner(const Field<Real>& u, const Field<Real>& v) const;

  /// @return sqrt of the integral of |u - v|^2, from the pointwise difference
  [[nodiscard]] Real Distance(const Field<Real>& u, const Field<Real>& v) const;

  /// @return sqrt of the integral of |u|^2
  [[nodiscard]] Real Norm(const Field<Real>& u) const;

private:
  class Form;
  std::unique_ptr<const Form> form_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_L2_INNER_PRODUCT_HPP
