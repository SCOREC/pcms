#ifndef PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
#define PCMS_TRANSFER_TRANSFER_OPERATOR_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_metadata.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms
{

template <typename T>
class TransferOperator
{
public:
  using value_type = T;

  /// @brief Writes `source` transferred to `target`.
  /// @throws pcms_error if a field is Dual, the value bases (units
  ///         included) differ, or `source` and `target` share their data
  void Apply(const Field<T>& source, Field<T>& target) const
  {
    if (&source.GetData() == &target.GetData()) {
      throw pcms_error("TransferOperator::Apply: source and target share "
                       "their data");
    }
    detail::CheckDuality(source.GetMetadata(), target.GetMetadata(), false);
    ApplyImpl(source, target);
  }

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
  /// @throws pcms_error unless SupportsTranspose(), or if a field is
  ///         Primal or `source` and `target` share their data
  void ApplyTranspose(const Field<T>& source, Field<T>& target) const
  {
    if (&source.GetData() == &target.GetData()) {
      throw pcms_error("TransferOperator::ApplyTranspose: source and target "
                       "share their data");
    }
    detail::CheckDuality(source.GetMetadata(), target.GetMetadata(), true);
    ApplyTransposeImpl(source, target);
  }

  virtual ~TransferOperator() noexcept = default;

private:
  virtual void ApplyImpl(const Field<T>& source, Field<T>& target) const = 0;

  virtual void ApplyTransposeImpl(const Field<T>& /*source*/,
                                  Field<T>& /*target*/) const
  {
    throw pcms_error("TransferOperator::ApplyTranspose: not supported by this "
                     "operator");
  }
};

} // namespace pcms

#endif // PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
