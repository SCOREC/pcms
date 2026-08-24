#ifndef PCMS_TRANSFER_TRANSFER_H
#define PCMS_TRANSFER_TRANSFER_H

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

  // internal-use Apply that can work on Rank2Views. This is needed
  // for wrapping transfer operators. e.g., for coordinate transformations
  virtual void Apply(TransferKey, const Field<T>& source,
                     Rank2View<T, DeviceMemorySpace> out) const = 0;

  virtual ~TransferOperator() noexcept = default;

protected:
  [[nodiscard]] static TransferKey MakeTransferKey() noexcept
  {
    return TransferKey{};
  }
};

} // namespace pcms

#endif
