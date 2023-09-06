#ifndef PCMS_COUPLING_TYPES_H
#define PCMS_COUPLING_TYPES_H
#include <redev_types.h>
namespace pcms
{
enum class Type
{
  Real,
  LO,
  GO
};
using Real = double;
using LO = int32_t;
using GO = int64_t;
template <typename T>
constexpr Type TypeEnumFromType(T)
{
  if constexpr (std::is_same_v<T, Real>) {
    return Type::Real;
  } else if constexpr (std::is_same_v<T, LO>) {
    return Type::LO;
  } else if constexpr (std::is_same_v<T, GO>) {
    return Type::GO;
  }
};

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

} // namespace pcms

#endif // PCMS_COUPLING_TYPES_H
