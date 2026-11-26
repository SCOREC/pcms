// Part of pcms/field/units.hpp, which includes it after the declaration of
// Units.
#include "pcms/field/units.hpp"
#ifndef PCMS_FIELD_DETAIL_UNITS_HPP
#define PCMS_FIELD_DETAIL_UNITS_HPP

namespace pcms::detail
{

/// @return whether `a` and `b` are equal, scale and offset included; always
///         true without PCMS_ENABLE_UNITS
[[nodiscard]] bool SameUnits(const Units& a, const Units& b);

/// @return whether `a` converts to `b` by a scale and an offset; always true
///         without PCMS_ENABLE_UNITS
[[nodiscard]] bool SameDimensions(const Units& a, const Units& b);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_UNITS_HPP
