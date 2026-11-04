#ifndef PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
#define PCMS_TRANSFER_TRANSFER_OPERATOR_HPP

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

  virtual ~TransferOperator() noexcept = default;
};

} // namespace pcms

#endif // PCMS_TRANSFER_TRANSFER_OPERATOR_HPP
