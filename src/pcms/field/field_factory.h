#ifndef PCMS_FIELD_FACTORY_H
#define PCMS_FIELD_FACTORY_H

#include "field.h"
#include "field_data.h"
#include "field_layout.h"
#include "field_metadata.h"
#include "pcms/utility/common.h"
#include "pcms/utility/types.h"
#include <memory>

namespace pcms
{

namespace detail
{

inline size_t ExpectedFlatFieldDataSize(const FieldLayout& layout)
{
  return static_cast<size_t>(layout.GetNumOwnedDofHolder()) *
         static_cast<size_t>(layout.GetNumComponents());
}

} // namespace detail

// Compile-time gate: true only for the five supported field value types.
template <typename T>
inline constexpr bool is_supported_field_type_v =
  std::is_same_v<T, int8_t> || std::is_same_v<T, int32_t> ||
  std::is_same_v<T, int64_t> || std::is_same_v<T, float> ||
  std::is_same_v<T, double>;

// FieldFactory is an abstract type that can produce
// Field<T>s over a known layout. This is the entire surface the coupler needs
// (GetLayout + CreateField). A FieldFactory makes discrete fields (values on
// DOF holders); it makes no claim about evaluating them as functions.
//
// Communication-only adapters (MFEM, XGC) are FieldFactories. Evaluatable
// spaces are FunctionSpaces, which derive from FieldFactory
class FieldFactory
{
public:
  virtual std::shared_ptr<const FieldLayout> GetLayout() const noexcept = 0;

  virtual ~FieldFactory() noexcept = default;

  // Create a new field with freshly allocated data for this factory.
  // Compile-time error for unsupported T; runtime error for T unsupported by
  // the concrete backend.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(FieldMetadata metadata = {}) const;

  // Expert API: wrap externally constructed field data into a Field for this
  // factory. The concrete factory validates backend-specific field-data type
  // and storage size compatibility.
  template <typename T>
  [[nodiscard]] Field<T> CreateField(std::unique_ptr<FieldData<T>> data) const;

protected:
  template <typename T>
  static Field<T> WrapField(std::shared_ptr<const FieldLayout> layout,
                            std::unique_ptr<FieldData<T>> data)
  {
    return Field<T>(typename Field<T>::CtorKey{}, std::move(layout),
                    std::move(data));
  }

  virtual FieldVariant CreateFieldImpl(Type value_type,
                                       FieldMetadata metadata) const = 0;

  virtual FieldVariant CreateFieldImpl(FieldDataVariant data) const = 0;
};

template <typename T>
Field<T> FieldFactory::CreateField(FieldMetadata metadata) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  return std::get<Field<T>>(CreateFieldImpl(TypeEnumFromType<T>(), metadata));
}

template <typename T>
Field<T> FieldFactory::CreateField(std::unique_ptr<FieldData<T>> data) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!data) {
    throw pcms_error("FieldFactory::CreateField: data must not be null");
  }
  return std::get<Field<T>>(CreateFieldImpl(FieldDataVariant{std::move(data)}));
}

} // namespace pcms

#endif // PCMS_FIELD_FACTORY_H
