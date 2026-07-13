#ifndef PCMS_FIELD_INTERPOLATOR_H
#define PCMS_FIELD_INTERPOLATOR_H

#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/transfer/field_compatibility.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <transfer_operator.hpp>

namespace pcms
{

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
template <typename T>
class Interpolator : public TransferOperator<T>
{
public:
  // Expensive: localizes target DOF coords into source mesh. Called once.
  Interpolator(const FunctionSpace& source_space,
               const FunctionSpace& target_space, OutOfBoundsPolicy policy = {})
    : num_points_(static_cast<LO>(target_space.GetLayout()
                                    ->GetDOFHolderCoordinates()
                                    .GetCoordinates()
                                    .extent(0))),
      n_comp_(target_space.GetLayout()->GetNumComponents()),
      source_layout_(source_space.GetLayout()),
      target_layout_(target_space.GetLayout()),
      evaluator_(source_space.CreatePointEvaluator<T>(
        EvaluationRequest::FromFunctionSpace(target_space, policy)))
  {
  }

  // Cheap: apply to any Field whose layout matches the target space.
  // Localization is not repeated.
  void Apply(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    detail::CheckTransferFieldLayout(source, *source_layout_, "source");
    detail::CheckTransferFieldLayout(target, *target_layout_, "target");
    const LO num_points = num_points_;
    const int n_comp = n_comp_;
    Kokkos::View<T**, DeviceMemorySpace> output("interp_output", num_points,
                                                n_comp);
    auto output_view = MakeRank2View(output);
    evaluator_->Evaluate(source, output_view);
    Kokkos::View<T*, DeviceMemorySpace> flat(
      "interp_flat", static_cast<size_t>(num_points) * n_comp);
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_LAMBDA(LO i) {
        for (int c = 0; c < n_comp; ++c) {
          flat(i * n_comp + c) = output(i, c);
        }
      });
    target.GetData().SetDOFHolderData(make_const_array_view(flat));
  }

private:
  LO num_points_;
  int n_comp_;
  std::shared_ptr<const FieldLayout> source_layout_;
  std::shared_ptr<const FieldLayout> target_layout_;
  std::unique_ptr<PointEvaluator<T>> evaluator_;
};

} // namespace pcms

#endif // PCMS_FIELD_INTERPOLATOR_H
