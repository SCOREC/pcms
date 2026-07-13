#ifndef PCMS_TRANSFER_FIELD_COMPATIBILITY_H
#define PCMS_TRANSFER_FIELD_COMPATIBILITY_H

#include "pcms/discretization/discretization.h"
#include "pcms/field/field.h"
#include "pcms/field/field_layout.h"
#include "pcms/utility/common.h"
#include <string>

namespace pcms
{
namespace detail
{

// Verify that a field handed to a TransferOperator::Apply belongs to the same
// space the operator was constructed from.
// Throws pcms_error on mismatch.
template <typename T>
inline void CheckTransferFieldLayout(const Field<T>& field,
                                     const FieldLayout& expected,
                                     const char* role)
{
  const FieldLayout& actual = field.GetLayout();
  if (&actual == &expected) {
    return;
  }
  auto actual_disc = actual.GetDiscretization();
  auto expected_disc = expected.GetDiscretization();
  if (actual_disc && expected_disc &&
      actual_disc->SameEntities(*expected_disc)) {
    return;
  }
  throw pcms_error(std::string("TransferOperator: ") + role +
                   " field does not match the function space the operator was "
                   "constructed from");
}

} // namespace detail
} // namespace pcms

#endif // PCMS_TRANSFER_FIELD_COMPATIBILITY_H
