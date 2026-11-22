#ifndef PCMS_FIELD_IDENTIFICATION_HPP
#define PCMS_FIELD_IDENTIFICATION_HPP

#include "pcms/field/field.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <vector>

namespace pcms
{

/// @brief The map of a twist-shift identification: per length crossed along
///        the identified axis, coordinate `shift_axis` moves by S at the
///        coordinates `input_axes`.
/// @param shift S, a scalar Function on a geometry of dimension
///        input_axes.size() in another coordinate system
/// @throws pcms_error if an axis is negative, `shift` is null or its
///         geometry's dimension is not the number of input axes; Bind throws
///         if an axis is not one of the coordinate system
[[nodiscard]] std::shared_ptr<const IdentificationMap> CreateTwistShiftMap(
  int shift_axis, std::shared_ptr<const Function<Real>> shift,
  std::vector<int> input_axes);

} // namespace pcms

#endif // PCMS_FIELD_IDENTIFICATION_HPP
