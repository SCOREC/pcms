#ifndef PCMS_FIELD_FIELD_OPERATOR_HPP
#define PCMS_FIELD_FIELD_OPERATOR_HPP

namespace pcms
{

template <typename>
class Field;

/// An operation on a single field, applied in place. Contrast with
/// TransferOperator, which maps one field onto another.
template <typename T>
class FieldOperator
{
public:
  virtual void Apply(Field<T>& field) const = 0;
  virtual ~FieldOperator() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_FIELD_OPERATOR_HPP
