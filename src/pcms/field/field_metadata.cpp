#include "pcms/field/field_metadata.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::detail
{

void CheckDuality(const FieldMetadata& source, const FieldMetadata& target,
                  bool transpose)
{
  const Duality refused = transpose ? Duality::Primal : Duality::Dual;
  for (const auto* side : {&source, &target}) {
    if (side->duality == refused) {
      throw pcms_error(
        std::string(transpose ? "transposed transfer" : "transfer") +
        " from '" + source.name + "' to '" + target.name + "': '" + side->name +
        "' is " +
        (transpose ? "primal; ApplyTranspose takes dual values"
                   : "dual; Apply takes primal values, use the "
                     "transpose to sum parts"));
    }
  }
}

} // namespace pcms::detail
