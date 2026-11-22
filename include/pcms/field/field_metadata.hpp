#ifndef PCMS_FIELD_FIELD_METADATA_HPP
#define PCMS_FIELD_FIELD_METADATA_HPP

#include <cstdint>
#include <optional>
#include <string>

namespace pcms
{

/// @brief How a quantity scales with the size of its region: Apply takes
///        intensive values, ApplyTranspose extensive ones.
enum class Extent : std::uint8_t
{
  Intensive,
  Extensive,
};

/// @brief Name and extent of a field; an unset extent is unknown and never
///        checked.
struct FieldMetadata
{
  std::string name;
  std::optional<Extent> extent;
};

} // namespace pcms

#include "pcms/field/detail/field_metadata.hpp"

#endif // PCMS_FIELD_FIELD_METADATA_HPP
