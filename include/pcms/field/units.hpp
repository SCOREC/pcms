#ifndef PCMS_FIELD_UNITS_HPP
#define PCMS_FIELD_UNITS_HPP

#include <memory>
#include <string>
#include <string_view>

namespace pcms
{

/// @brief The convention of a unit string: SI notation, CF/UDUNITS, or IMAS
///        (SI plus "-" for dimensionless).
enum class UnitSyntax
{
  Default,
  CF,
  IMAS,
};

class Units;

} // namespace pcms

#include "pcms/field/detail/units.hpp"

namespace pcms
{

/// @brief Physical units. Without PCMS_ENABLE_UNITS a Units holds its text
///        and is never compared.
class Units
{
public:
  /// @brief Parses a unit string.
  /// @param text the unit, e.g. "T" or "m2 s-1"
  /// @param syntax the convention of `text`
  /// @throws pcms_error if `text` is not a unit in `syntax` (IMAS "mixed" and
  ///         "as_parent" included); without PCMS_ENABLE_UNITS, if it is empty
  [[nodiscard]] static Units Parse(std::string_view text,
                                   UnitSyntax syntax = UnitSyntax::Default);

  /// @return a canonical form that Parse reads back; without
  ///         PCMS_ENABLE_UNITS, the trimmed text given to Parse
  [[nodiscard]] std::string ToString() const;

private:
  friend bool detail::SameUnits(const Units&, const Units&);
  friend bool detail::SameDimensions(const Units&, const Units&);

  struct Impl;
  explicit Units(std::shared_ptr<const Impl> impl);

  std::shared_ptr<const Impl> impl_;
};

} // namespace pcms

#endif // PCMS_FIELD_UNITS_HPP
