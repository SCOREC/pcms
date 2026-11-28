#ifndef PCMS_FIELD_FIELD_METADATA_HPP
#define PCMS_FIELD_FIELD_METADATA_HPP

#include <cstdint>
#include <optional>
#include <string>

namespace pcms
{

/// @brief Primal: the values or coefficients of a function, the input and
///        output of Apply. Dual: functionals E^T W f, the input and
///        output of ApplyTranspose.
enum class Duality : std::uint8_t
{
  Primal,
  Dual,
};

/// @brief Name and duality of a field; an unset duality is unknown and never
///        checked.
struct FieldMetadata
{
  std::string name;
  std::optional<Duality> duality;
};

} // namespace pcms

#include "pcms/field/detail/field_metadata.hpp"

#endif // PCMS_FIELD_FIELD_METADATA_HPP
