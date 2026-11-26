#ifndef PCMS_FIELD_FIELD_FACTORY_HPP
#define PCMS_FIELD_FIELD_FACTORY_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <optional>
#include <string>
#include "pcms/field/detail/field_factory.hpp"

namespace pcms
{

// FieldFactory constructs Field<T> bundles ({layout, data}) with no evaluation
// capability and no shared-ownership requirement. Comm-only backends (e.g.
// XGCFieldFactory) derive it directly. FunctionSpace is a separate abstraction
// (not a FieldFactory) that additionally supports evaluation and produces
// Functions; a FieldFactory can never be passed where a FunctionSpace is
// required, which is the compile-time form of "comm-only cannot be evaluated".
class FieldFactory
{
public:
  virtual std::shared_ptr<const FieldLayout> GetLayout() const noexcept = 0;

  virtual ~FieldFactory() noexcept = default;

  /// @brief Creates a Field.
  /// @param metadata its name, units and duality
  /// @param basis its values' basis
  template <typename T>
  [[nodiscard]] Field<T> CreateField(FieldMetadata metadata,
                                     ValueBasis basis) const;

  /// @brief Creates a Field of scalar values, or of independent components
  ///        (ValueKind::Array) for a layout of several components.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(FieldMetadata metadata) const;

  /// @brief Creates a Field named `name`; see the overloads with metadata.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(std::string name, ValueBasis basis) const;

  /// @brief Creates a Field named `name` in the default basis.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(std::string name = "") const;

  // Expert API: wrap externally constructed field data into a Field. The
  // concrete factory validates backend-specific field-data type and storage
  // size compatibility.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(FieldMetadata metadata,
                                     std::unique_ptr<FieldData<T>> data) const;

protected:
  template <typename T>
  static Field<T> WrapField(std::shared_ptr<const FieldLayout> layout,
                            std::unique_ptr<FieldData<T>> data)
  {
    return Field<T>(FieldMetadata{}, std::move(layout), std::move(data));
  }

  virtual FieldVariant CreateFieldImpl(Type storage_type,
                                       ValueBasis basis) const = 0;

  virtual FieldVariant CreateFieldImpl(FieldDataVariant data) const = 0;
};

template <typename T>
Field<T> FieldFactory::CreateField(std::string name, ValueBasis basis) const
{
  return CreateField<T>(FieldMetadata{.name = std::move(name)},
                        std::move(basis));
}

template <typename T>
Field<T> FieldFactory::CreateField(std::string name) const
{
  return CreateField<T>(FieldMetadata{.name = std::move(name)});
}

template <typename T>
Field<T> FieldFactory::CreateField(FieldMetadata metadata) const
{
  return CreateField<T>(
    std::move(metadata),
    detail::DefaultValueBasis(GetLayout()->GetNumComponents()));
}

template <typename T>
Field<T> FieldFactory::CreateField(FieldMetadata metadata,
                                   ValueBasis basis) const
{
  static_assert(detail::is_supported_field_type_v<T>,
                "T is not a supported field type");
  basis = detail::ResolveValueBasis(GetLayout()->GetCoordinateSystem(),
                                    std::move(basis));
  detail::ValidateValueSemantics(basis, GetLayout()->GetNumComponents());
  Field<T> field = std::get<Field<T>>(
    CreateFieldImpl(detail::TypeEnumFromType<T>(), std::move(basis)));
  field.metadata_ = std::move(metadata);
  return field;
}

template <typename T>
Field<T> FieldFactory::CreateField(FieldMetadata metadata,
                                   std::unique_ptr<FieldData<T>> data) const
{
  static_assert(detail::is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!data) {
    throw pcms_error("FieldFactory::CreateField: data must not be null");
  }
  detail::ValidateValueSemantics(data->GetValueBasis(),
                                 GetLayout()->GetNumComponents());
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
  field.metadata_ = std::move(metadata);
  return field;
}

} // namespace pcms

#endif // PCMS_FIELD_FIELD_FACTORY_HPP
