#ifndef PCMS_EVALUATION_REQUEST_H
#define PCMS_EVALUATION_REQUEST_H

#include "pcms/field/coordinate_view.hpp"
#include "field_layout.h"
#include "out_of_bounds_policy.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"

#include <memory>

namespace pcms
{

class FunctionSpace;

struct EvaluationRequest
{
  // Query coordinates to evaluate at. These are always consumed at evaluator
  // construction time.
  CoordinateView<DeviceMemorySpace> coords;
  // Optional provenance for the query sites. When present, construction-time
  // logic may use the layout and its discretization to select optimized
  // localization paths. Concrete PointEvaluator implementations are not
  // required to retain this layout after construction.
  std::shared_ptr<const FieldLayout> query_layout;
  // Construction-time policy that is baked into the created PointEvaluator.
  OutOfBoundsPolicy policy = {};

  static EvaluationRequest FromCoordinates(
    CoordinateView<DeviceMemorySpace> coords, OutOfBoundsPolicy policy = {});

  static EvaluationRequest FromLayout(std::shared_ptr<const FieldLayout> layout,
                                      OutOfBoundsPolicy policy = {});

  static EvaluationRequest FromFunctionSpace(
    const FunctionSpace& function_space, OutOfBoundsPolicy policy = {});

  [[nodiscard]] const FieldLayout* GetQueryLayout() const noexcept
  {
    return query_layout.get();
  }

  [[nodiscard]] const Discretization* GetQueryDiscretization() const noexcept
  {
    auto* layout = GetQueryLayout();
    if (layout == nullptr) {
      return nullptr;
    }
    auto disc = layout->GetDiscretization();
    return disc.get();
  }

private:
  explicit EvaluationRequest(CoordinateView<DeviceMemorySpace> coords_in,
                             std::shared_ptr<const FieldLayout> query_layout_in,
                             OutOfBoundsPolicy policy_in = {}) noexcept
    : coords(coords_in),
      query_layout(std::move(query_layout_in)),
      policy(policy_in)
  {
  }
};

} // namespace pcms

#endif // PCMS_EVALUATION_REQUEST_H
