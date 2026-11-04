#ifndef PCMS_FUNCTION_SPACE_H
#define PCMS_FUNCTION_SPACE_H

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_factory.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <variant>

namespace pcms
{

// A FunctionSpace is an evaluatable field space: the layout of its fields, the
// evaluator source that evaluates them at points, and the coordinate
// interpretation. It produces Functions (a Field bound to this space).
//
// It is the one concrete space; backends differ in the components they hand
// to its constructor, and the Create* factories (CreateLagrangeSpace, ...)
// assemble them. A Function co-owns its space via shared_ptr (obtained through
// enable_shared_from_this), so spaces are always shared-owned:
// CreateFunction throws for a space that is not.
class FunctionSpace final
  : public PointEvaluatorFactory
  , public std::enable_shared_from_this<FunctionSpace>
{
public:
  /// A function space from its components. Own it by std::shared_ptr:
  /// functions co-own their space. Backends build spaces this way; users call
  /// the Create* factories.
  /// @param layout the layout of the space's fields
  /// @param create_field_data makes the data of a field of this space from a
  ///        value type and basis (may throw for spaces whose fields hold no
  ///        values in PCMS)
  /// @param evaluator the evaluator source that evaluates the fields at
  ///        points, or null for a space that does not evaluate at points
  /// @param geometry the geometry of a space without an evaluator source, or
  ///        null (an evaluator source supplies its own)
  explicit FunctionSpace(
    std::shared_ptr<const FieldLayout> layout,
    std::function<FieldDataVariant(Type, ValueBasis)> create_field_data,
    std::shared_ptr<const PointEvaluatorFactory> evaluator,
    std::shared_ptr<const Geometry> geometry = nullptr) noexcept
    : layout_(std::move(layout)),
      create_field_data_(std::move(create_field_data)),
      evaluator_(std::move(evaluator)),
      geometry_(std::move(geometry))
  {
    PCMS_ALWAYS_ASSERT(evaluator_ == nullptr || geometry_ == nullptr);
  }

  /// @return the geometry the space is created on: its evaluator source's, or
  ///         the one given to the constructor; null for a space that does not
  ///         evaluate at points
  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return evaluator_ != nullptr ? evaluator_->GetGeometry() : geometry_;
  }

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout() const noexcept
  {
    return layout_;
  }

  /// @return the geometry's coordinate system, or the layout's for a space
  ///         without a geometry
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept
  {
    const Geometry* geometry =
      evaluator_ != nullptr ? evaluator_->GetGeometry().get() : geometry_.get();
    return geometry != nullptr ? geometry->GetCoordinateSystem()
                               : GetLayout()->GetCoordinateSystem();
  }

  /// @return whether the space's evaluator source supports `op`; for a space
  ///         without one, whether `op` is Value and the space has a geometry
  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override;

  ~FunctionSpace() noexcept override = default;

  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name = "", VarianceSignature variance = values::Scalar,
    std::optional<ComponentScaling> scaling = std::nullopt) const;

  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name, VarianceSignature variance,
    std::shared_ptr<const CoordinateSystem> system,
    std::optional<ComponentScaling> scaling = std::nullopt) const;

  // Expert API: wrap externally constructed field data into a Function for this
  // space. The space validates backend-specific field-data type and storage
  // size compatibility.
  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name, std::unique_ptr<FieldData<T>> data) const;

private:
  /// @brief Creates the evaluator through the evaluator source.
  /// @throws pcms_error if the space has no evaluator source
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op) const override;

  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  std::shared_ptr<const FunctionSpace> Owner() const;

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

  FieldVariant CreateFieldImpl(Type storage_type, ValueBasis basis) const;

  FieldVariant CreateFieldImpl(FieldDataVariant data) const;

  std::shared_ptr<const FieldLayout> layout_;
  std::function<FieldDataVariant(Type, ValueBasis)> create_field_data_;
  std::shared_ptr<const PointEvaluatorFactory> evaluator_;
  std::shared_ptr<const Geometry> geometry_;
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
  auto owner = Owner();
  ValueBasis basis = detail::MakeValueBasis(GetCoordinateSystem(), variance,
                                            std::move(system), scaling);
  detail::ValidateValueSemantics(basis, GetLayout()->GetNumComponents());
  Field<T> field = std::get<Field<T>>(
    CreateFieldImpl(TypeEnumFromType<T>(), std::move(basis)));
  return WrapFunction<T>(std::move(name), std::move(field.layout_),
                         std::move(field.data_), std::move(owner));
}

template <typename T>
Function<T> FunctionSpace::CreateFunction(
  std::string name, std::unique_ptr<FieldData<T>> data) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  auto owner = Owner();
  if (!data) {
    throw pcms_error("FunctionSpace::CreateFunction: data must not be null");
  }
  detail::ValidateValueSemantics(data->GetValueBasis(),
                                 GetLayout()->GetNumComponents());
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
  return WrapFunction<T>(std::move(name), std::move(field.layout_),
                         std::move(field.data_), std::move(owner));
}

} // namespace pcms

#endif // PCMS_FUNCTION_SPACE_H
