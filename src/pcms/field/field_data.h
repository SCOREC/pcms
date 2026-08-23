#ifndef PCMS_FIELD_DATA_H
#define PCMS_FIELD_DATA_H

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <memory>
#include <variant>

namespace pcms
{

template <typename T>
class FieldData
{
public:
  using value_type = T;

  virtual FieldValueType GetValueType() const = 0;
  virtual const ValueBasis& GetValueBasis() const = 0;

  virtual Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const = 0;
  virtual void SetDOFHolderDataHost(
    Rank2View<const T, HostMemorySpace> values) = 0;

  virtual Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const = 0;
  virtual void SetDOFHolderData(
    Rank2View<const T, DeviceMemorySpace> values) = 0;

  virtual ~FieldData() noexcept = default;
};

using FieldDataVariant = std::variant<
  std::unique_ptr<FieldData<int8_t>>, std::unique_ptr<FieldData<int32_t>>,
  std::unique_ptr<FieldData<int64_t>>, std::unique_ptr<FieldData<float>>,
  std::unique_ptr<FieldData<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_DATA_H
