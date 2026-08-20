#ifndef PCMS_FUNCTION_SPACE_H
#define PCMS_FUNCTION_SPACE_H

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "evaluation_request.h"
#include "field.h"
#include "field_data.h"
#include "field_factory.h"
#include "field_layout.h"
#include "pcms/field/value_view.hpp"
#include "out_of_bounds_policy.h"
#include "point_evaluator.h"
#include "pcms/discretization/discretization.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <memory>
#include <optional>
#include <variant>

namespace pcms
{

// FunctionSpace is an abstract interface representing an evaluatable field
// space: layout, evaluation rules, and coordinate interpretation. It produces
// Functions (a Field bound to this space).
//
// Concrete implementations (e.g. LagrangeFunctionSpace) provide backends for
// specific discretizations or mesh types, and their From* factories return a
// shared_ptr. A Function co-owns its space via shared_ptr (obtained through
// enable_shared_from_this), which is why concrete spaces are always
// shared-owned; a stack-allocated FunctionSpace is not constructible (the
// backend constructors are private and reachable only via the shared_ptr
// From* factories).
class FunctionSpace : public std::enable_shared_from_this<FunctionSpace>
{
public:
  virtual std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept
  {
    return GetLayout()->GetDiscretization();
  }

  virtual std::shared_ptr<const FieldLayout> GetLayout() const noexcept = 0;

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const
  {
    return GetLayout()->GetCoordinateSystem();
  }

  virtual ~FunctionSpace() noexcept = default;

  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name = "", VarianceSignature variance = values::Scalar,
    std::optional<ComponentScaling> scaling = std::nullopt) const;

  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name, VarianceSignature variance,
    std::shared_ptr<const CoordinateSystem> system,
    std::optional<ComponentScaling> scaling = std::nullopt) const;

  // Expert API: wrap externally constructed field data into a Function for this
  // space. The concrete space validates backend-specific field-data type and
  // storage size compatibility.
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name, std::unique_ptr<FieldData<T>> data) const;

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
  template <typename T>
  static Field<T> WrapField(std::shared_ptr<const FieldLayout> layout,
                            std::unique_ptr<FieldData<T>> data)
  {
    return Field<T>(std::string{}, std::move(layout), std::move(data));
  }

  template <typename T>
  static Function<T> WrapFunction(std::string name,
                                  std::shared_ptr<const FieldLayout> layout,
                                  std::unique_ptr<FieldData<T>> data,
                                  std::shared_ptr<const FunctionSpace> space)
  {
    return Function<T>(std::move(name), std::move(layout), std::move(data),
                       std::move(space));
  }

  virtual FieldVariant CreateFieldImpl(Type storage_type,
                                       ValueBasis basis) const = 0;

  virtual FieldVariant CreateFieldImpl(FieldDataVariant data) const = 0;

  virtual PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const = 0;
};

template <typename T>
Function<T> FunctionSpace::CreateFunction(
  std::string name, VarianceSignature variance,
  std::optional<ComponentScaling> scaling) const
{
  return CreateFunction<T>(std::move(name), variance, nullptr, scaling);
}

template <typename T>
Function<T> FunctionSpace::CreateFunction(
  std::string name, VarianceSignature variance,
  std::shared_ptr<const CoordinateSystem> system,
  std::optional<ComponentScaling> scaling) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  ValueBasis basis = detail::MakeValueBasis(GetCoordinateSystem(), variance,
                                            std::move(system), scaling);
  detail::ValidateValueSemantics(basis, GetLayout()->GetNumComponents());
  Field<T> field = std::get<Field<T>>(
    CreateFieldImpl(TypeEnumFromType<T>(), std::move(basis)));
  return WrapFunction<T>(std::move(name), std::move(field.layout_),
                         std::move(field.data_), shared_from_this());
}

template <typename T>
Function<T> FunctionSpace::CreateFunction(
  std::string name, std::unique_ptr<FieldData<T>> data) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!data) {
    throw pcms_error("FunctionSpace::CreateFunction: data must not be null");
  }
  detail::ValidateValueSemantics(data->GetValueBasis(),
                                 GetLayout()->GetNumComponents());
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
  return WrapFunction<T>(std::move(name), std::move(field.layout_),
                         std::move(field.data_), shared_from_this());
}

template <typename T>
std::unique_ptr<PointEvaluator<T>> FunctionSpace::CreatePointEvaluator(
  const EvaluationRequest& request) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!SameCoordinateSystem(request.coords.GetCoordinateSystem(),
                            GetCoordinateSystem())) {
    throw pcms_error(
      "CreatePointEvaluator: query coordinate system '" +
      std::string(request.coords.GetCoordinateSystem()->Kind()) +
      "' is not the space's coordinate system '" +
      std::string(GetCoordinateSystem()->Kind()) +
      "' (the coordinate system its mesh coordinates were declared in; "
      "coordinate systems compare by "
      "object identity); map the points with a CoordinateMap and "
      "evaluate at the View() of its result");
  }
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
