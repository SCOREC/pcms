// Part of pcms/field/field_metadata.hpp, which includes it after the
// declarations it needs.
#include "pcms/field/field_metadata.hpp"
#ifndef PCMS_FIELD_DETAIL_FIELD_METADATA_HPP
#define PCMS_FIELD_DETAIL_FIELD_METADATA_HPP

namespace pcms::detail
{

/// @brief Checks the dualities against the direction of the transfer.
/// @throws pcms_error if a side is Dual for Apply or Primal for
///         ApplyTranspose
void CheckDuality(const FieldMetadata& source, const FieldMetadata& target,
                  bool transpose);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_METADATA_HPP
