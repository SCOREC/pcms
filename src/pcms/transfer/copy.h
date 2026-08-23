#ifndef PCMS_TRANSFER_FIELD2_H_
#define PCMS_TRANSFER_FIELD2_H_
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/function_space.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include "pcms/transfer/transfer_operator.hpp"

namespace pcms
{

namespace detail
{

template <typename T>
void CheckCopyCompatible(const Field<T>& source, const Field<T>& target)
{
  if (&source.GetLayout() != &target.GetLayout()) {
    throw pcms_error("Copy: source and target layouts differ");
  }
  const auto& sd = source.GetData();
  const auto& td = target.GetData();
  if (sd.GetValueType() != td.GetValueType() ||
      !SameValueBasis(sd.GetValueBasis(), td.GetValueBasis())) {
    throw pcms_error("Copy: source and target value type/basis differ");
  }
}

} // namespace detail

template <typename T>
class Copy : public TransferOperator<T>
{
public:
  Copy(const FunctionSpace& source_space, const FunctionSpace& target_space)
  {
    auto source_layout = source_space.GetLayout();
    auto target_layout = target_space.GetLayout();
    if (source_layout.get() != target_layout.get()) {
      throw pcms_error("Copy: source and target function spaces have "
                       "different layouts");
    }
  }

  void Apply(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    detail::CheckCopyCompatible(source, target);
    target.SetDOFHolderDataHost(source.GetDOFHolderDataHost());
  }
};

} // namespace pcms

#endif // PCMS_TRANSFER_FIELD2_H_
