#ifndef PCMS_TRANSFER_TRANSFER_H
#define PCMS_TRANSFER_TRANSFER_H

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"

namespace pcms
{

template <typename>
class Field;

template <typename>
class TransferOperator;

// passkey for handling internal optimized Apply path
class TransferKey
{
  explicit TransferKey() = default;
  template <typename>
  friend class TransferOperator;
};

template <typename T>
class TransferOperator
{
public:
  using value_type = T;

  virtual void Apply(const Field<T>& source, Field<T>& target) const = 0;

  // internal-use Apply that writes into a caller-supplied buffer instead of a
  // Field. This is needed for wrapping transfer operators, e.g. for
  // coordinate transformations.
  //
  // With no target Field to consult, the output buffer carries its own value
  // basis: the caller tags the ValueView with the basis it expects the
  // operator to write, and the operator throws pcms_error unless that matches
  // the basis it actually writes. Operators that only move components write
  // the source's stored basis; an operator that re-expresses them writes the
  // basis it transforms into. All rank-0 bases compare equal, so scalar
  // callers may tag with ValueBasis{}.
  virtual void Apply(TransferKey, const Field<T>& source,
                     ValueView<T, DeviceMemorySpace> out) const = 0;

  virtual ~TransferOperator() noexcept = default;

protected:
  [[nodiscard]] static TransferKey MakeTransferKey() noexcept
  {
    return TransferKey{};
  }

  /// Throws pcms_error unless the output view claims the basis this operator
  /// writes.
  static void CheckApplyWriteTag(const char* context,
                                 const ValueView<T, DeviceMemorySpace>& out,
                                 const ValueBasis& written)
  {
    detail::CheckWrittenValueBasis(context, out.GetBasis(), written);
  }
};

} // namespace pcms

#endif
