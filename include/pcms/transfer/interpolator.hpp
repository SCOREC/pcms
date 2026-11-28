#ifndef PCMS_TRANSFER_INTERPOLATOR_HPP
#define PCMS_TRANSFER_INTERPOLATOR_HPP

#include "pcms/field/function_space.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/bound_points.hpp"
#include "pcms/transfer/detail/interpolator.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"
#include <memory>

namespace pcms
{

/// Point interpolation of `source` at the DOF holders of `target_space`. The
/// target's DOF-holder coordinates are localized once, here; Apply only
/// evaluates, so call it repeatedly in a coupling loop. The operator holds
/// `source` and `target_space`.
///
///   auto interp = CreateInterpolator(src_space, tgt_space);
///   interp->Apply(src_field, tgt_field);
///
/// @param source evaluated at the target's DOF-holder coordinates: a space or
///        a derivative source (Gradient, Divergence, Curl)
/// @param target_space space whose DOF holders receive the values; its basis
///        must be the basis the source's evaluators write
/// @param policy out-of-bounds policy for the evaluation
/// @throws pcms_error if `source` or `target_space` is null, the source has no
///         geometry, or the target's DOFs are not values at its DOF holders
template <typename T = Real>
[[nodiscard]] std::unique_ptr<TransferOperator<T>> CreateInterpolator(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space,
  OutOfBoundsPolicy policy = {})
{
  return std::make_unique<detail::Interpolator<T>>(
    std::move(source), std::move(target_space), policy);
}

/// As above, reusing a localization shared with other interpolators of the
/// same source geometry at the same target sites:
///
///   auto layout = tgt_space->GetLayout();
///   auto bound = src_space->GetGeometry()->Bind(
///     layout->GetDOFHolderCoordinates(), DOFHolderOwners(*layout));
///   auto values = CreateInterpolator(src_space, tgt_space, bound);
///   auto gradient = CreateInterpolator(Gradient(src_space), grad_space,
///   bound);
///
/// @param bound result of binding the source's geometry to the DOF-holder
///        coordinates of `target_space` or of a space with the same DOF-holder
///        coordinates
/// @throws pcms_error if `bound` is not bound to the target's DOF-holder sites
template <typename T = Real>
[[nodiscard]] std::unique_ptr<TransferOperator<T>> CreateInterpolator(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space,
  std::shared_ptr<const BoundPoints> bound)
{
  return std::make_unique<detail::Interpolator<T>>(
    std::move(source), std::move(target_space), std::move(bound));
}

} // namespace pcms

#endif // PCMS_TRANSFER_INTERPOLATOR_HPP
