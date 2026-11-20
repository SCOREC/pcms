#ifndef PCMS_TRANSFER_GALERKIN_PROJECTION_HPP
#define PCMS_TRANSFER_GALERKIN_PROJECTION_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/transfer/linear_solver.hpp"
#include "pcms/transfer/projection.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms::transfer
{

// Galerkin L2 projection onto a function space through quadratures:
//   M x = E^T W (f - g) + M g
// with E the target space evaluated at the quadrature points, W the weights,
// f the source at the points, and g the control variate (the interpolation of
// the source onto the target space, only with options.control_variate; zero
// otherwise). The right-hand side is a LinearForm (E^T W f); M = E^T W E is a
// BilinearForm on the natural quadrature of the target space, solved by a
// LinearSolver (or, lumped, divided by its row sums).
//
// Construction (one-time cost): the quadrature, the evaluators of the source
// and the target at its points, the mass form and the solver's setup.
// Apply() (per-call cost): the source at the points, the transpose, a solve.
class GalerkinProjection : public TransferOperator<Real>
{
public:
  /// Holds `source` and `target_space`.
  GalerkinProjection(std::shared_ptr<const PointEvaluatorFactory> source,
                     std::shared_ptr<const FunctionSpace> target_space,
                     ProjectionOptions options);

  void Apply(const Field<Real>& source, Field<Real>& target) const override;

  /// Supported without a control variate when the source's evaluator at the
  /// quadrature points supports the transpose: P^T = S^T W E M^-1 (M is
  /// symmetric; the same solver), with S the source and E the target at the
  /// points. `source` is data on the target's layout, `target` a functional on
  /// the source space's layout.
  [[nodiscard]] bool SupportsTranspose() const noexcept override;
  void ApplyTranspose(const Field<Real>& source,
                      Field<Real>& target) const override;

private:
  std::shared_ptr<const PointEvaluatorFactory> source_;
  std::shared_ptr<const FunctionSpace> target_space_;
  detail::LinearForm right_hand_side_;
  std::unique_ptr<PointEvaluator<Real>> source_at_points_;
  detail::BilinearForm mass_;
  // exactly one of the two: the solver of M x = b, or, for the lumped mass,
  // its row sums (x = b / row sums)
  std::unique_ptr<detail::LinearSolver> solver_;
  std::unique_ptr<Function<Real>> lumped_;
  // the control variate: the interpolation of the source onto the target
  // space, whose scratch field is overwritten on every Apply
  struct ControlVariate
  {
    ControlVariate(std::shared_ptr<const PointEvaluatorFactory> source,
                   std::shared_ptr<const FunctionSpace> target_space,
                   std::shared_ptr<const BoundPoints> at_holders)
      : interpolator(std::move(source), target_space, std::move(at_holders)),
        values(target_space->CreateFunction<Real>())
    {
    }
    pcms::detail::Interpolator<Real> interpolator;
    mutable Function<Real> values;
  };
  std::optional<ControlVariate> control_variate_;
  mutable Function<Real> functional_;
  mutable Function<Real> coefficients_; // M^-1 b, for the transpose
  mutable Kokkos::View<Real**, DeviceMemorySpace> samples_; // [points][comps]
};

} // namespace pcms::transfer

#endif // PCMS_TRANSFER_GALERKIN_PROJECTION_HPP
