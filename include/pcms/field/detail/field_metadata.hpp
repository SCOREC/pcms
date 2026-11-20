// Part of pcms/field/field_metadata.hpp, which includes it after the
// declarations it needs.
#include "pcms/field/field_metadata.hpp"
#ifndef PCMS_FIELD_DETAIL_FIELD_METADATA_HPP
#define PCMS_FIELD_DETAIL_FIELD_METADATA_HPP

#include <string>

namespace pcms::detail
{

/// @brief Checks that the units of `source` and `target` agree.
/// @param operation names the transfer in the message
/// @throws pcms_error if both have units and they differ; never without
///         PCMS_ENABLE_UNITS
void CheckSameUnits(const FieldMetadata& source, const FieldMetadata& target,
                    const std::string& operation);

/// @brief Checks the extents against the direction of the transfer.
/// @throws pcms_error if a side is Extensive for Apply or Intensive for
///         ApplyTranspose
void CheckExtent(const FieldMetadata& source, const FieldMetadata& target,
                 bool transpose);

/// @brief Sets the unknown units and extent of `target` from `source`.
void LearnMetadata(const FieldMetadata& source, FieldMetadata& target);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_METADATA_HPP
