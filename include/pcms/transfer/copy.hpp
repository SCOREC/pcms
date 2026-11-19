#ifndef PCMS_TRANSFER_COPY_HPP
#define PCMS_TRANSFER_COPY_HPP

#include "pcms/field/function_space.hpp"
#include "pcms/transfer/detail/copy.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"
#include <memory>

namespace pcms
{

/// Copies DOF values between spaces with the same DOF holders (the same
/// layout). The operator holds both spaces.
/// @throws pcms_error if a space is null or the spaces have different layouts
template <typename T = Real>
[[nodiscard]] std::unique_ptr<TransferOperator<T>> CreateCopy(
  std::shared_ptr<const FunctionSpace> source_space,
  std::shared_ptr<const FunctionSpace> target_space)
{
  return std::make_unique<detail::Copy<T>>(std::move(source_space),
                                           std::move(target_space));
}

} // namespace pcms

#endif // PCMS_TRANSFER_COPY_HPP
