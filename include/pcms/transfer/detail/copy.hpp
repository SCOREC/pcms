#ifndef PCMS_TRANSFER_DETAIL_COPY_HPP
#define PCMS_TRANSFER_DETAIL_COPY_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_ops.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include <memory>

namespace pcms::detail
{

// Copies DOF values between fields of spaces with the same layout. Create it
// with CreateCopy (pcms/transfer/copy.hpp).
template <typename T>
class Copy : public TransferOperator<T>
{
public:
  /// @throws pcms_error if a space is null or the spaces have different
  ///         layouts
  Copy(std::shared_ptr<const FunctionSpace> source_space,
       std::shared_ptr<const FunctionSpace> target_space)
    : source_space_(std::move(source_space)),
      target_space_(std::move(target_space))
  {
    if (source_space_ == nullptr || target_space_ == nullptr) {
      throw pcms_error("Copy: the function spaces must not be null");
    }
    if (source_space_->GetLayout().get() != target_space_->GetLayout().get()) {
      throw pcms_error("Copy: source and target function spaces have "
                       "different layouts");
    }
  }

private:
  void ApplyImpl(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    pcms::Copy(source, target);
  }

  std::shared_ptr<const FunctionSpace> source_space_;
  std::shared_ptr<const FunctionSpace> target_space_;
};

} // namespace pcms::detail

#endif // PCMS_TRANSFER_DETAIL_COPY_HPP
