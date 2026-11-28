// Part of pcms/transfer/transfer_operator.hpp, which includes it after
// declaring TransferOperator and before defining it.
#include "pcms/transfer/transfer_operator.hpp"
#ifndef PCMS_TRANSFER_DETAIL_TRANSFER_OPERATOR_HPP
#define PCMS_TRANSFER_DETAIL_TRANSFER_OPERATOR_HPP

namespace pcms::detail
{

/// @brief target = A^T source + beta target, with the checks of
///        TransferOperator::ApplyTranspose; beta == 0 does not read the
///        target
template <typename T>
void ApplyTranspose(const TransferOperator<T>& op, const Field<T>& source,
                    Field<T>& target, T beta);

} // namespace pcms::detail

#endif // PCMS_TRANSFER_DETAIL_TRANSFER_OPERATOR_HPP
