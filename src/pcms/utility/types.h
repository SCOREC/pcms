#ifndef PCMS_COUPLING_TYPES_H
#define PCMS_COUPLING_TYPES_H
#include "assert.h"
#include <cstdint>
#include <type_traits>

namespace pcms
{
enum class Type
{
  Real,
  LO,
  GO,
  Int8,
  Float
};
using Real = double;
using LO = int32_t;
using GO = int64_t;

namespace detail
{
template <typename... T>
struct dependent_always_false : std::false_type
{};

template <typename T>
struct type_identity
{
  using type = T;
};
template <typename T>
using type_identity_t = typename type_identity<T>::type;
} // namespace detail

template <typename T>
constexpr Type TypeEnumFromType()
{
  if constexpr (std::is_same_v<T, double>) {
    return Type::Real;
  } else if constexpr (std::is_same_v<T, int32_t>) {
    return Type::LO;
  } else if constexpr (std::is_same_v<T, int64_t>) {
    return Type::GO;
  } else if constexpr (std::is_same_v<T, int8_t>) {
    return Type::Int8;
  } else if constexpr (std::is_same_v<T, float>) {
    return Type::Float;
  } else {
    static_assert(detail::dependent_always_false<T>::value,
                  "T is not a supported field type");
  }
};

// Dispatches on a runtime Type value by instantiating F with the corresponding
// type tag. F must accept detail::type_identity<T> for each of the five
// supported scalar types and return a consistent type.
template <typename F>
auto apply_to_type(Type t, F&& f)
{
  switch (t) {
    case Type::Int8: return f(detail::type_identity<int8_t>{});
    case Type::LO: return f(detail::type_identity<int32_t>{});
    case Type::GO: return f(detail::type_identity<int64_t>{});
    case Type::Float: return f(detail::type_identity<float>{});
    case Type::Real: return f(detail::type_identity<double>{});
  }
  throw pcms_error("apply_to_type: unhandled Type value");
}

} // namespace pcms

#endif // PCMS_COUPLING_TYPES_H
