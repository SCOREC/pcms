#ifndef PCMS_FIELD_INTERPOLATOR_H
#define PCMS_FIELD_INTERPOLATOR_H

#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include "pcms/transfer/transfer_operator.hpp"

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
    : num_points_(static_cast<LO>(
        target_space.GetLayout()->GetDOFHolderCoordinates().GetValues().extent(
          0))),
      n_comp_(target_space.GetLayout()->GetNumComponents()),
      evaluator_(source_space.CreatePointEvaluator<T>(
        EvaluationRequest::FromFunctionSpace(target_space, policy))),
      target_values_("interp_output", num_points_, n_comp_)
  {
  }

  void Apply(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& sd = source.GetData();
    const auto& td = target.GetData();
    if (sd.GetValueType() != td.GetValueType()) {
      throw pcms_error("Interpolator: source and target value types differ");
    }
    if (!SameValueBasis(sd.GetValueBasis(), td.GetValueBasis())) {
      throw pcms_error(
        "Interpolator: the source's stored basis differs from the target's "
        "declared basis");
    }
    Apply(this->MakeTransferKey(), source, MakeRank2View(target_values_));
    target.SetDOFHolderDataUnchecked(MakeConstRank2View(target_values_));
  }

  void Apply(TransferKey, const Field<T>& source,
             Rank2View<T, DeviceMemorySpace> out) const override
  {
    PCMS_FUNCTION_TIMER;
    if (static_cast<LO>(out.extent(0)) != num_points_ ||
        static_cast<int>(out.extent(1)) != n_comp_) {
      throw pcms_error(
        "Interpolator: output buffer extents do not match the target layout");
    }
    evaluator_->Evaluate(source, out);
  }

private:
  LO num_points_;
  int n_comp_;
  std::unique_ptr<PointEvaluator<T>> evaluator_;
  mutable Kokkos::View<T**, DeviceMemorySpace> target_values_;
};

} // namespace pcms

#endif // PCMS_FIELD_INTERPOLATOR_H
