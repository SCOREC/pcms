#ifndef PCMS_FIELD_FACTORY_H
#define PCMS_FIELD_FACTORY_H

#include "field.h"
#include "field_data.h"
#include "field_layout.h"
#include "field_metadata.h"
#include "pcms/utility/types.h"
#include <memory>
#include <string>

namespace pcms
{

namespace detail
{

inline size_t ExpectedFlatFieldDataSize(const FieldLayout& layout)
{
  return static_cast<size_t>(layout.GetNumLocalDofHolder()) *
         static_cast<size_t>(layout.GetNumComponents());
}

} // namespace detail

// Compile-time gate: true only for the five supported field value types.
template <typename T>
inline constexpr bool is_supported_field_type_v =
  std::is_same_v<T, int8_t> || std::is_same_v<T, int32_t> ||
  std::is_same_v<T, int64_t> || std::is_same_v<T, float> ||
  std::is_same_v<T, double>;

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

  // Create a new named field with freshly allocated data. The name is optional
  // (empty by default) and identifies the field to consumers such as a coupler.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(std::string name = "",
                                     FieldMetadata metadata = {}) const;

  // Expert API: wrap externally constructed field data into a Field. The
  // concrete factory validates backend-specific field-data type and storage
  // size compatibility.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(std::string name,
                                     std::unique_ptr<FieldData<T>> data) const;

protected:
  template <typename T>
  static Field<T> WrapField(std::shared_ptr<const FieldLayout> layout,
                            std::unique_ptr<FieldData<T>> data)
  {
    return Field<T>(std::string{}, std::move(layout), std::move(data));
  }

  virtual FieldVariant CreateFieldImpl(Type value_type,
                                       FieldMetadata metadata) const = 0;

  virtual FieldVariant CreateFieldImpl(FieldDataVariant data) const = 0;
};

template <typename T>
Field<T> FieldFactory::CreateField(std::string name,
                                   FieldMetadata metadata) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(TypeEnumFromType<T>(), metadata));
  field.name_ = std::move(name);
  return field;
}

template <typename T>
Field<T> FieldFactory::CreateField(std::string name,
                                   std::unique_ptr<FieldData<T>> data) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!data) {
    throw pcms_error("FieldFactory::CreateField: data must not be null");
  }
  Field<T> field =
    std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
  field.name_ = std::move(name);
  return field;
}

} // namespace pcms

#endif // PCMS_FIELD_FACTORY_H
