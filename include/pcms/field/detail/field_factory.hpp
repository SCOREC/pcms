#ifndef PCMS_FIELD_DETAIL_FIELD_FACTORY_HPP
#define PCMS_FIELD_DETAIL_FIELD_FACTORY_HPP

#include "pcms/field/field_layout.hpp"
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace pcms::detail
{

// Compile-time gate: true only for the five supported field value types.
template <typename T>
inline constexpr bool is_supported_field_type_v =
  std::is_same_v<T, int8_t> || std::is_same_v<T, int32_t> ||
  std::is_same_v<T, int64_t> || std::is_same_v<T, float> ||
  std::is_same_v<T, double>;

inline size_t ExpectedFlatFieldDataSize(const FieldLayout& layout)
{
  return static_cast<size_t>(layout.GetNumDofHolder()) *
         static_cast<size_t>(layout.GetNumComponents());
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_FACTORY_HPP
