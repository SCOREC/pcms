#ifndef PCMS_FIELD_FUNCTION_SPACE_HPP
#define PCMS_FIELD_FUNCTION_SPACE_HPP

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
#include "pcms/geometry/quadrature.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <variant>
#include <string>

namespace pcms
{

// A FunctionSpace is an evaluatable field space: the layout of its fields, the
// evaluator source that evaluates them at points, and the coordinate
// interpretation. It produces Functions (a Field bound to this space).
class FunctionSpace
  : public PointEvaluatorFactory
  , public std::enable_shared_from_this<FunctionSpace>
{
public:
  /// A function space from its components. Own it by std::shared_ptr:
  /// functions co-own their space. Backends build spaces this way; users call
  /// the Create* factories.
  /// @param layout the layout of the space's fields
  /// @param evaluator the evaluator source that evaluates the fields at
  ///        points, or null for a space whose fields are not evaluated at
  ///        points
  /// @param geometry the geometry the space is on: the evaluator source's,
  ///        when there is one
  /// @throws pcms_error if `geometry` is null or is not the evaluator
  ///         source's geometry
  explicit FunctionSpace(std::shared_ptr<const FieldLayout> layout,
                         std::shared_ptr<const PointEvaluatorFactory> evaluator,
                         std::shared_ptr<const Geometry> geometry)
    : layout_(std::move(layout)),
      evaluator_(std::move(evaluator)),
      geometry_(std::move(geometry))
  {
    if (geometry_ == nullptr) {
      throw pcms_error("FunctionSpace: the geometry must not be null");
    }
    if (evaluator_ != nullptr && evaluator_->GetGeometry() != geometry_) {
      throw pcms_error("FunctionSpace: the geometry must be the evaluator "
                       "source's geometry");
    }
  }

  /// @return the geometry the space is created on
  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return geometry_;
  }

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout() const noexcept
  {
    return layout_;
  }

  /// @return the geometry's coordinate system
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept
  {
    return geometry_->GetCoordinateSystem();
  }

  /// @brief The space's natural quadrature on its geometry, which
  ///        should be exact for the product of two functions in this space
  /// @throws pcms_error if the space has none
  [[nodiscard]] const Quadrature& GetQuadrature() const
  {
    if (!quadrature_) {
      quadrature_ = CreateQuadrature();
    }
    return *quadrature_;
  }

  /// @return whether the space's evaluator source supports `op`
  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override;

  /// @return whether the space's evaluator source has nodal DOFs
  [[nodiscard]] bool HasNodalDOFs() const noexcept override
  {
    return evaluator_ != nullptr && evaluator_->HasNodalDOFs();
  }

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

  /// @brief Creates a Function of scalar values.
  /// @param metadata its name and duality
  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(FieldMetadata metadata) const;

  // Expert API: wrap externally constructed field data into a Function for this
  // space. The space validates backend-specific field-data type and storage
  // size compatibility.
  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  template <typename T>
  [[nodiscard]] Function<T> CreateFunction(
    std::string name, std::unique_ptr<FieldData<T>> data) const;

protected:
  /// @brief The space's natural quadrature (GetQuadrature).
  /// @throws pcms_error unless the kind overrides it
  [[nodiscard]] virtual Quadrature CreateQuadrature() const
  {
    throw pcms_error("FunctionSpace::GetQuadrature: this space has no "
                     "quadrature");
  }

  /// @brief Creates the data of a field of this space.
  /// @throws pcms_error for an invalid storage_type
  [[nodiscard]] virtual FieldDataVariant CreateFieldData(
    Type storage_type, ValueBasis basis) const;

  /// @brief Checks wrapped field data based on backend requirements
  /// @throws pcms_error if the data is not a field data of this space
  virtual void ValidateFieldData(const FieldDataVariant& data) const;

private:
  /// @brief Creates the evaluator through the evaluator source.
  /// @throws pcms_error if the space has no evaluator source
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op) const override;

  /// @throws pcms_error if the space is not owned by a std::shared_ptr
  std::shared_ptr<const FunctionSpace> Owner() const;

  template <typename T>
  static Field<T> WrapField(std::shared_ptr<const FieldLayout> layout,
                            std::unique_ptr<FieldData<T>> data)
  {
    return Field<T>(FieldMetadata{}, std::move(layout), std::move(data));
  }

  template <typename T>
  static Function<T> WrapFunction(FieldMetadata metadata,
                                  std::shared_ptr<const FieldLayout> layout,
                                  std::unique_ptr<FieldData<T>> data,
                                  std::shared_ptr<const FunctionSpace> space)
  {
    return Function<T>(std::move(metadata), std::move(layout), std::move(data),
                       std::move(space));
  }

  FieldVariant CreateFieldImpl(Type storage_type, ValueBasis basis) const;

  FieldVariant CreateFieldImpl(FieldDataVariant data) const;

  std::shared_ptr<const FieldLayout> layout_;
  std::shared_ptr<const PointEvaluatorFactory> evaluator_;
  std::shared_ptr<const Geometry> geometry_;
  mutable std::optional<Quadrature> quadrature_;
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
  static_assert(detail::is_supported_field_type_v<T>,
                "T is not a supported field type");
  auto owner = Owner();
  ValueBasis basis = detail::MakeValueBasis(GetCoordinateSystem(), variance,
                                            std::move(system), scaling);
  detail::ValidateValueSemantics(basis, GetLayout()->GetNumComponents());
  Field<T> field = std::get<Field<T>>(
    CreateFieldImpl(detail::TypeEnumFromType<T>(), std::move(basis)));
  return WrapFunction<T>(FieldMetadata{.name = std::move(name)},
                         std::move(field.layout_), std::move(field.data_),
                         std::move(owner));
}

template <typename T>
Function<T> FunctionSpace::CreateFunction(FieldMetadata metadata) const
{
  static_assert(detail::is_supported_field_type_v<T>,
                "T is not a supported field type");
  auto owner = Owner();
  ValueBasis basis = detail::MakeValueBasis(
    GetCoordinateSystem(), values::Scalar, nullptr, std::nullopt);
  detail::ValidateValueSemantics(basis, GetLayout()->GetNumComponents());
  Field<T> field = std::get<Field<T>>(
    CreateFieldImpl(detail::TypeEnumFromType<T>(), std::move(basis)));
  return WrapFunction<T>(std::move(metadata), std::move(field.layout_),
                         std::move(field.data_), std::move(owner));
}

template <typename T>
Function<T> FunctionSpace::CreateFunction(
  std::string name, std::unique_ptr<FieldData<T>> data) const
{
  static_assert(detail::is_supported_field_type_v<T>,
                "T is not a supported field type");
  auto owner = Owner();
  if (!data) {
    throw pcms_error("FunctionSpace::CreateFunction: data must not be null");
  }
  detail::ValidateValueSemantics(data->GetValueBasis(),
                                 GetLayout()->GetNumComponents());
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
  return WrapFunction<T>(FieldMetadata{.name = std::move(name)},
                         std::move(field.layout_), std::move(field.data_),
                         std::move(owner));
}

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACE_HPP
