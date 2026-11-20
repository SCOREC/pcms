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

/// The L2 inner product of two fields on the same domain, by quadrature: the
/// quadrature of u's space (its own rule exact for u times a function of the
/// space, unless one is given) is bound once for both. The
/// components of a vector field are paired by the Euclidean product (no
/// metric yet). Setup is cached here; Apply and Distance only evaluate.
class L2InnerProduct
{
public:
  /// @param u_space the space of u, whose geometry the quadrature is on
  /// @param v the evaluator source of v (a space, or an operator on one),
  ///        bound at the quadrature's points; its fields have as many
  ///        components as u's
  /// @param quadrature a quadrature of u_space's geometry to use in place of
  ///        the space's own
  /// @throws pcms_error if no quadrature is given and the space has none, the
  ///         quadrature is not on the space's geometry, or the evaluators
  ///         cannot be created
  L2InnerProduct(std::shared_ptr<const FunctionSpace> u_space,
                 std::shared_ptr<const PointEvaluatorFactory> v,
                 std::optional<Quadrature> quadrature = {});

  /// @return the integral of u . v
  [[nodiscard]] Real Apply(const Field<Real>& u, const Field<Real>& v) const;

  /// @return sqrt of the integral of |u - v|^2, from the pointwise
  ///         difference, never from the expansion |u|^2 - 2 u v + |v|^2
  ///         (which loses all accuracy for nearly equal fields)
  [[nodiscard]] Real Distance(const Field<Real>& u, const Field<Real>& v) const;

private:
  std::shared_ptr<const FunctionSpace> u_space_;
  std::shared_ptr<const PointEvaluatorFactory> v_source_;
  Quadrature quadrature_;
  std::unique_ptr<PointEvaluator<Real>> u_at_points_;
  std::unique_ptr<PointEvaluator<Real>> v_at_points_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> u_values_, v_values_;
};

/// The integral and the L2 norm of a field, by quadrature (the rule as for
/// L2InnerProduct).
class L2Integral
{
public:
  /// @throws pcms_error as L2InnerProduct
  L2Integral(std::shared_ptr<const FunctionSpace> space,
             std::optional<Quadrature> quadrature = {});

  /// @return the integral of a scalar field
  /// @throws pcms_error for a field of several components
  [[nodiscard]] Real Apply(const Field<Real>& u) const;

  /// @return sqrt of the integral of |u|^2
  [[nodiscard]] Real Norm(const Field<Real>& u) const;

private:
  std::shared_ptr<const FunctionSpace> space_;
  Quadrature quadrature_;
  std::unique_ptr<PointEvaluator<Real>> at_points_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> values_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_L2_INNER_PRODUCT_HPP
