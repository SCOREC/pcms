#ifndef PCMS_TRANSFER_INTERPOLATOR_HPP
#define PCMS_TRANSFER_INTERPOLATOR_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include "pcms/transfer/transfer_operator.hpp"

namespace pcms
{

namespace detail
{

/// @brief Whether `bound` was bound to the DOF-holder coordinates of `layout`,
///        in the same order and coordinate system: its located points are
///        those coordinates, mapped by the bind's coordinate map if it had
///        one.
/// @return true if it was
bool SameDOFHolderSites(const BoundSource& bound, const FieldLayout& layout);

} // namespace detail

// Interpolator<T> separates the expensive localization step from the cheap
// repeated evaluation step, making it efficient to use in a coupling loop.
//
// Construct once per source×target function-space pair (localizes target DOF
// coordinates into the source mesh), then call Apply repeatedly for different
// field states at zero additional localization cost. Any Field sharing the
// same target FunctionSpace can be passed to Apply.
//
// Usage:
//   Interpolator<Real> interp(src_space, tgt_space);
//   interp.Apply(src_field, tgt_field);   // cheap; called in coupling loop
//   interp.Apply(src_field_next, tgt_field_next); // reuses cached localization
//
// Interpolators of the same source at the same target share one bound source:
//   auto layout = tgt_space->GetLayout();
//   auto geometry = src_space->GetGeometry();
//   auto bound = geometry->Bind(layout->GetDOFHolderCoordinates(),
//                               DOFHolderOwners(*layout));
//   Interpolator<Real> values(*src_space, *tgt_space, bound);
//   Interpolator<Real> gradient(*detail::Gradient(src_space), *grad_space,
//   bound);
template <typename T>
class Interpolator : public TransferOperator<T>
{
public:
  // Expensive: localizes target DOF coords into source mesh. Called once.
  /// @param source factory evaluated at the target's DOF-holder coordinates
  /// @param target_space space whose DOF holders receive the values
  /// @param policy out-of-bounds policy for the evaluation
  /// @throws pcms_error if the source has no geometry
  Interpolator(const PointEvaluatorFactory& source,
               const FunctionSpace& target_space, OutOfBoundsPolicy policy = {})
    : Interpolator(source, target_space,
                   BindTargetSites(source, target_space, policy))
  {
  }

  /// @brief Creates an interpolator from an existing bound source.
  /// @param source source evaluated at the target's DOF-holder coordinates
  /// @param target_space space whose DOF holders receive the values; its
  ///        basis must be the basis the source's evaluators write
  /// @param bound result of binding the source's geometry to the DOF-holder
  ///        coordinates of `target_space` or of a space with the same
  ///        DOF-holder coordinates
  Interpolator(const PointEvaluatorFactory& source,
               const FunctionSpace& target_space,
               std::shared_ptr<const BoundSource> bound)
    : num_points_(static_cast<LO>(
        target_space.GetLayout()->GetDOFHolderCoordinates().GetValues().extent(
          0))),
      n_comp_(target_space.GetLayout()->GetNumComponents()),
      evaluator_(CreateEvaluator(source, target_space, std::move(bound))),
      target_values_("interp_output", num_points_, n_comp_)
  {
  }

  void Apply(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& sd = source.GetData();
    const auto& td = target.GetData();
    const ValueBasis written = evaluator_->OutputBasis(sd.GetValueBasis());
    if (!detail::SameValueBasis(written, td.GetValueBasis())) {
      throw pcms_error(
        "Interpolator: the basis written for the source's stored basis "
        "differs from the target's declared basis");
    }
    evaluator_->Evaluate(source, ValueView<T, DeviceMemorySpace>(
                                   written, MakeRank2View(target_values_)));
    target.SetDOFHolderDataUnchecked(MakeConstRank2View(target_values_));
    target.Synchronize();
  }

private:
  static std::shared_ptr<const BoundSource> BindTargetSites(
    const PointEvaluatorFactory& source, const FunctionSpace& target_space,
    OutOfBoundsPolicy policy)
  {
    const auto geometry = source.GetGeometry();
    if (geometry == nullptr) {
      throw pcms_error("Interpolator: the source has no geometry");
    }
    const auto layout = target_space.GetLayout();
    return geometry->Bind(layout->GetDOFHolderCoordinates(),
                          DOFHolderOwners(*layout), policy);
  }

  static std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    const PointEvaluatorFactory& source, const FunctionSpace& target_space,
    std::shared_ptr<const BoundSource> bound)
  {
    if (bound == nullptr ||
        !detail::SameDOFHolderSites(*bound, *target_space.GetLayout())) {
      throw pcms_error("Interpolator: the bound source must be bound to a "
                       "space with the target's DOF-holder sites");
    }
    return source.CreatePointEvaluator<T>(std::move(bound));
  }

  LO num_points_;
  int n_comp_;
  std::unique_ptr<PointEvaluator<T>> evaluator_;
  mutable Kokkos::View<T**, DeviceMemorySpace> target_values_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_INTERPOLATOR_HPP
