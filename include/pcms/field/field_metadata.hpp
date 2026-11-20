#ifndef PCMS_FIELD_FIELD_METADATA_HPP
#define PCMS_FIELD_FIELD_METADATA_HPP

#include "pcms/field/units.hpp"
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

/// @brief Name, units and extent of a field; unset units and extent are
///        unknown and never checked.
struct FieldMetadata
{
  std::string name;
  std::optional<Units> units;
  std::optional<Extent> extent;
};

} // namespace pcms

#include "pcms/field/detail/field_metadata.hpp"

#endif // PCMS_FIELD_FIELD_METADATA_HPP
