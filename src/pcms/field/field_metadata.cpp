#include "pcms/field/field_metadata.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::detail
{

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

} // namespace pcms::detail
