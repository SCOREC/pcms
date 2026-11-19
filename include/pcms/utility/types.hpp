#ifndef PCMS_UTILITY_TYPES_HPP
#define PCMS_UTILITY_TYPES_HPP
#include "pcms/utility/assert.hpp"
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

} // namespace pcms

#include "pcms/utility/detail/types.hpp"

#endif // PCMS_UTILITY_TYPES_HPP
