#ifndef PCMS_FUNCTION_SPACE_H
#define PCMS_FUNCTION_SPACE_H

#include "coordinate_system.h"
#include "evaluation_request.h"
#include "field.h"
#include "field_data.h"
#include "field_factory.h"
#include "field_layout.h"
#include "field_metadata.h"
#include "out_of_bounds_policy.h"
#include "point_evaluator.h"
#include "pcms/discretization/discretization.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <memory>
#include <variant>

namespace pcms
{

// A FunctionSpace allows you to construct fields and evaluators for those fields
class FunctionSpace : public FieldFactory
{
public:
  virtual std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept
  {
    return GetLayout()->GetDiscretization();
  }

  virtual CoordinateSystem GetCoordinateSystem() const noexcept = 0;

  // Create a point evaluator for the given evaluation request.
  // Compile-time error for unsupported T; runtime error for T or capability
  // unsupported by the concrete backend.
  //
  // EvaluationRequest is a construction-time object: it supplies the query
  // coordinates, out-of-bounds policy, and any optional provenance that may
  // help the backend choose an optimized localization path. The resulting
  // PointEvaluator caches only the resolved state needed for repeated
  // Evaluate(...) calls; it is not required to retain the original request.
  template <typename T>
  [[nodiscard]] std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    const EvaluationRequest& request) const;

protected:
  virtual PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const = 0;
};

template <typename T>
std::unique_ptr<PointEvaluator<T>> FunctionSpace::CreatePointEvaluator(
  const EvaluationRequest& request) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  return std::get<std::unique_ptr<PointEvaluator<T>>>(
    CreatePointEvaluatorImpl(TypeEnumFromType<T>(), request));
}

inline EvaluationRequest EvaluationRequest::FromCoordinates(
  CoordinateView<DeviceMemorySpace> coords, OutOfBoundsPolicy policy)
{
  return EvaluationRequest(coords, nullptr, policy);
}

inline EvaluationRequest EvaluationRequest::FromLayout(
  std::shared_ptr<const FieldLayout> layout, OutOfBoundsPolicy policy)
{
  if (layout == nullptr) {
    throw pcms_error("EvaluationRequest::FromLayout: layout must not be null");
  }
  // Must evaluate GetDOFHolderCoordinates() before std::move(layout)
  auto coords = layout->GetDOFHolderCoordinates();
  return EvaluationRequest(coords, std::move(layout), policy);
}

inline EvaluationRequest EvaluationRequest::FromFunctionSpace(
  const FunctionSpace& function_space, OutOfBoundsPolicy policy)
{
  return FromLayout(function_space.GetLayout(), policy);
}

} // namespace pcms

#endif // PCMS_FUNCTION_SPACE_H
