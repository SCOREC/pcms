#ifndef PCMS_TRANSFER_DETAIL_INTERPOLATOR_HPP
#define PCMS_TRANSFER_DETAIL_INTERPOLATOR_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

namespace pcms::detail
{

/// @brief Whether `bound` was bound to the DOF-holder coordinates of `layout`,
///        in the same order and coordinate system: its located points are
///        those coordinates, mapped by the bind's coordinate map if it had
///        one.
/// @return true if it was
bool SameDOFHolderSites(const BoundPoints& bound, const FieldLayout& layout);

// Interpolator<T> separates the expensive localization step from the cheap
// repeated evaluation step, making it efficient to use in a coupling loop.
// Create it with CreateInterpolator (pcms/transfer/interpolator.hpp).
//
// Construct once per source×target function-space pair (localizes target DOF
// coordinates into the source mesh), then call Apply repeatedly for different
// field states at zero additional localization cost. Any Field sharing the
// same target FunctionSpace can be passed to Apply.
template <typename T>
class Interpolator : public TransferOperator<T>
{
public:
  // Expensive: localizes target DOF coords into source mesh. Called once.
  /// @param source factory evaluated at the target's DOF-holder coordinates
  /// @param target_space space whose DOF holders receive the values
  /// @param policy out-of-bounds policy for the evaluation
  /// @throws pcms_error if the source or the target space is null, or the
  ///         source has no geometry
  Interpolator(std::shared_ptr<const PointEvaluatorFactory> source,
               std::shared_ptr<const FunctionSpace> target_space,
               OutOfBoundsPolicy policy = {})
    : Interpolator(source, target_space,
                   BindTargetSites(source, target_space, policy))
  {
  }

  /// @brief Creates an interpolator from an existing BoundPoints.
  /// @param source source evaluated at the target's DOF-holder coordinates
  /// @param target_space space whose DOF holders receive the values; its
  ///        basis must be the basis the source's evaluators write
  /// @param bound result of binding the source's geometry to the DOF-holder
  ///        coordinates of `target_space` or of a space with the same
  ///        DOF-holder coordinates
  /// @throws pcms_error if the source or the target space is null
  Interpolator(std::shared_ptr<const PointEvaluatorFactory> source,
               std::shared_ptr<const FunctionSpace> target_space,
               std::shared_ptr<const BoundPoints> bound)
    : source_(Require(std::move(source), "source")),
      target_space_(Require(std::move(target_space), "target space")),
      num_points_(static_cast<LO>(target_space_->GetLayout()
                                    ->GetDOFHolderCoordinates()
                                    .GetValues()
                                    .extent(0))),
      n_comp_(target_space_->GetLayout()->GetNumComponents()),
      evaluator_(CreateEvaluator(*source_, *target_space_, std::move(bound))),
      target_values_("interp_output", num_points_, n_comp_)
  {
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return evaluator_->SupportsTranspose();
  }

private:
  void ApplyImpl(const Field<T>& source, Field<T>& target) const override
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

  /// The source evaluator's EvaluateTranspose: `source` holds data at the
  /// target's DOF holders (loads), `target` is a functional of the source
  /// space's layout. The totals of the loads are preserved when the
  /// interpolation reproduces constants.
  void ApplyTransposeImpl(const Field<T>& source,
                          Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    if (&source.GetLayout() != target_space_->GetLayout().get()) {
      throw pcms_error("Interpolator::ApplyTranspose: the source must be a "
                       "field of the target space's layout");
    }
    detail::SetZero(target);
    evaluator_->EvaluateTranspose(ValueView<const T, DeviceMemorySpace>(
                                    target.GetData().GetValueBasis(),
                                    source.GetDOFHolderData().GetValues()),
                                  target);
  }

  template <typename P>
  static P Require(P p, const char* what)
  {
    if (p == nullptr) {
      throw pcms_error(std::string("Interpolator: the ") + what +
                       " must not be null");
    }
    return p;
  }

  static std::shared_ptr<const BoundPoints> BindTargetSites(
    const std::shared_ptr<const PointEvaluatorFactory>& source,
    const std::shared_ptr<const FunctionSpace>& target_space,
    OutOfBoundsPolicy policy)
  {
    const auto geometry = Require(source, "source")->GetGeometry();
    if (geometry == nullptr) {
      throw pcms_error("Interpolator: the source has no geometry");
    }
    const auto layout = Require(target_space, "target space")->GetLayout();
    return geometry->Bind(layout->GetDOFHolderCoordinates(),
                          DOFHolderOwners(*layout), policy);
  }

  static std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    const PointEvaluatorFactory& source, const FunctionSpace& target_space,
    std::shared_ptr<const BoundPoints> bound)
  {
    if (bound == nullptr ||
        !detail::SameDOFHolderSites(*bound, *target_space.GetLayout())) {
      throw pcms_error("Interpolator: the BoundPoints must be bound to a "
                       "space with the target's DOF-holder sites");
    }
    return source.CreatePointEvaluator<T>(std::move(bound));
  }

  std::shared_ptr<const PointEvaluatorFactory> source_;
  std::shared_ptr<const FunctionSpace> target_space_;
  LO num_points_;
  int n_comp_;
  std::unique_ptr<PointEvaluator<T>> evaluator_;
  mutable Kokkos::View<T**, DeviceMemorySpace> target_values_;
};

} // namespace pcms::detail

#endif // PCMS_TRANSFER_DETAIL_INTERPOLATOR_HPP
