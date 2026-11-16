#ifndef PCMS_OUT_OF_BOUNDS_POLICY_H
#define PCMS_OUT_OF_BOUNDS_POLICY_H

#include "pcms/utility/types.h"

namespace pcms
{

enum class OutOfBoundsMode
{
  ERROR,           // Throw error when points are out of bounds
  FILL,            // Fill out-of-bounds points with a fill value
  NEAREST_BOUNDARY // Map to nearest boundary cell (extrapolate)
};

// OutOfBoundsPolicy wraps OutOfBoundsMode with an optional fill value into a
// single construction-time policy object. Out-of-bounds behavior is a policy
// passed when creating a PointEvaluator rather than mutable state on the field.
//
// NearestBoundary is an optional backend capability. Call
// FieldEvaluatorFactory::SupportsNearestBoundary() before requesting it.
// Backends that do not support it throw a descriptive error if it is requested.
// The policy is fixed when the PointEvaluator is created; it is not mutable
// after that evaluator has been constructed.

struct OutOfBoundsPolicy
{
  OutOfBoundsMode mode = OutOfBoundsMode::ERROR;
  Real fill_value = 0.0; // used only when mode == FILL
};

} // namespace pcms

#endif // PCMS_OUT_OF_BOUNDS_POLICY_H
