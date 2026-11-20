#include "pcms/field/field_metadata.hpp"
#include "pcms/configuration.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::detail
{

void CheckSameUnits(const FieldMetadata& source, const FieldMetadata& target,
                    const std::string& operation)
{
#ifdef PCMS_ENABLE_UNITS
  if (!source.units || !target.units ||
      SameUnits(*source.units, *target.units)) {
    return;
  }
  const std::string what = operation + " from '" + source.name + "' [" +
                           source.units->ToString() + "] to '" + target.name +
                           "' [" + target.units->ToString() + "]: ";
  if (!SameDimensions(*source.units, *target.units)) {
    throw pcms_error(what + "the units have different dimensions");
  }
  throw pcms_error(what + "the units differ in scale or offset; PCMS does not "
                          "convert units, so convert one of the fields");
#else
  (void)source;
  (void)target;
  (void)operation;
#endif
}

void CheckExtent(const FieldMetadata& source, const FieldMetadata& target,
                 bool transpose)
{
  const Extent refused = transpose ? Extent::Intensive : Extent::Extensive;
  for (const auto* side : {&source, &target}) {
    if (side->extent == refused) {
      throw pcms_error(
        std::string(transpose ? "transposed transfer" : "transfer") +
        " from '" + source.name + "' to '" + target.name + "': '" + side->name +
        "' is " +
        (transpose ? "intensive; ApplyTranspose takes extensive values"
                   : "extensive; Apply takes intensive values, use the "
                     "transpose to sum parts"));
    }
  }
}

void LearnMetadata(const FieldMetadata& source, FieldMetadata& target)
{
  if (!target.units) {
    target.units = source.units;
  }
  if (!target.extent) {
    target.extent = source.extent;
  }
}

} // namespace pcms::detail
