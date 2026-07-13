#ifndef PCMS_TRANSFER_TRANSFER_H
#define PCMS_TRANSFER_TRANSFER_H

namespace pcms
{

template <typename>
class Field;
class FunctionSpace;

// A TransferOperator maps a source field to a target field -- both discrete
// fields on function spaces. That space-to-space contract is what distinguishes
// a transfer from a PointEvaluator (field-on-a-space -> values at arbitrary
// coordinates), which is the layer to reach for when the target is not a field
// on a space. Every operator is therefore built from a (source, target) space
// pair, recorded here so callers can check a field's space against the one the
// operator was built for before applying it.
template <typename T>
class TransferOperator
{
public:
  virtual void Apply(const Field<T>& source, Field<T>& target) const = 0;
  virtual ~TransferOperator() noexcept = default;

  [[nodiscard]] const FunctionSpace& SourceSpace() const noexcept
  {
    return *source_space_;
  }
  [[nodiscard]] const FunctionSpace& TargetSpace() const noexcept
  {
    return *target_space_;
  }

protected:
  TransferOperator(const FunctionSpace& source_space,
                   const FunctionSpace& target_space) noexcept
    : source_space_(&source_space), target_space_(&target_space)
  {
  }

private:
  const FunctionSpace* source_space_;
  const FunctionSpace* target_space_;
};

} // namespace pcms

#endif
