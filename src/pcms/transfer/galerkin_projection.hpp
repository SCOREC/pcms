#ifndef PCMS_TRANSFER_GALERKIN_PROJECTION_HPP
#define PCMS_TRANSFER_GALERKIN_PROJECTION_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/transfer/linear_solver.hpp"
#include "pcms/transfer/detail/projection.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include <memory>
#include <optional>

namespace pcms::transfer
{

// Galerkin L2 projection onto a function space V through quadratures:
//   x = M^-1 B_VS u,   M = B(V) on V's natural quadrature,
// with B_VS the form of the source S on the right-hand side's quadrature.
// With a control variate (the interpolation I of the source onto V, or, for a
// V without nodal DOFs, its projection with V's natural quadrature), only the
// residual is projected:
//   x = I u + M^-1 (B_VS u - B_VV I u),
// with B_VV the form of V on the right-hand side's points.
class GalerkinProjection : public TransferOperator<Real>
{
public:
  /// Holds `source` and `target_space`.
  /// @throws pcms_error if the target space's evaluator does not support the
  ///         transpose
  GalerkinProjection(std::shared_ptr<const PointEvaluatorFactory> source,
                     std::shared_ptr<const FunctionSpace> target_space,
                     pcms::detail::ProjectionOptions options);

  /// Supported when the source's evaluator at the points (and the control
  /// variate's interpolation) support the transpose:
  /// P^T v = B_VS^T z + I^T (v - B_VV z), z = M^-1 v. `source` is data on the
  /// target's layout, `target` a functional on the source space's layout.
  [[nodiscard]] bool SupportsTranspose() const noexcept override;

private:
  void ApplyImpl(const Field<Real>& source, Field<Real>& target) const override;
  void ApplyTransposeImpl(const Field<Real>& source, Field<Real>& target,
                          Real beta) const override;

  /// @return `scratch`, a function of the target space in `basis`
  Function<Real>& InBasis(Function<Real>& scratch,
                          const ValueBasis& basis) const;

  std::shared_ptr<const FunctionSpace> target_space_;
  detail::BilinearForm right_hand_side_;
  detail::BilinearForm mass_;
  // the solver of M x = b
  std::unique_ptr<detail::LinearSolver> solver_;
  // the control variate: I, and B_VV on the right-hand side's points
  std::unique_ptr<TransferOperator<Real>> interpolation_;
  std::optional<detail::BilinearForm> control_variate_;
  mutable Function<Real> functional_;
  mutable Function<Real> coefficients_;
  mutable Function<Real> correction_;
};

} // namespace pcms::transfer

#endif // PCMS_TRANSFER_GALERKIN_PROJECTION_HPP
