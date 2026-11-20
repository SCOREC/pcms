#ifndef PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
#define PCMS_TRANSFER_TRANSFER_OPERATOR_HPP

#include "pcms/utility/assert.hpp"

namespace pcms
{

template <typename>
class Field;

template <typename T>
class TransferOperator
{
public:
  using value_type = T;

  virtual void Apply(const Field<T>& source, Field<T>& target) const = 0;

  /// @return whether ApplyTranspose is implemented
  [[nodiscard]] virtual bool SupportsTranspose() const noexcept
  {
    return false;
  }

  /// target = A^T source, where A is Apply: the adjoint of the linear map on
  /// the DOF vectors, <A x, y> = <x, A^T y> (coefficient dot products). Its
  /// source is data on Apply's target layout (loads, sensitivities); its
  /// target is a functional on Apply's source layout, overwritten (unlike
  /// PointEvaluator::EvaluateTranspose, which accumulates).
  /// @throws pcms_error unless SupportsTranspose()
  virtual void ApplyTranspose(const Field<T>& /*source*/,
                              Field<T>& /*target*/) const
  {
    throw pcms_error("TransferOperator::ApplyTranspose: not supported by this "
                     "operator");
  }

  virtual ~TransferOperator() noexcept = default;
};

} // namespace pcms

#endif // PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
