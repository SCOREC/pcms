#ifndef PCMS_COUPLING_FIELD_H
#define PCMS_COUPLING_FIELD_H

#include "field_data.h"
#include "field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <memory>
#include <variant>

namespace pcms
{

class FieldFactory;

// Field<T> is FieldLayout (topology / coupling identity) plus owned FieldData<T>
template <typename T>
class Field
{
public:
  Field(Field&&) = default;
  Field& operator=(Field&&) = default;
  Field(const Field&) = delete;
  Field& operator=(const Field&) = delete;

  FieldData<T>& GetData() noexcept { return *data_; }
  const FieldData<T>& GetData() const noexcept { return *data_; }

  const FieldLayout& GetLayout() const { return *layout_; }

  Rank1View<const T, HostMemorySpace> GetDOFHolderDataHost() const
  {
    return data_->GetDOFHolderDataHost();
  }

  void SetDOFHolderDataHost(Rank1View<const T, HostMemorySpace> v)
  {
    data_->SetDOFHolderDataHost(v);
  }

  Rank1View<const T, DeviceMemorySpace> GetDOFHolderData() const
  {
    return data_->GetDOFHolderData();
  }

  void SetDOFHolderData(Rank1View<const T, DeviceMemorySpace> v)
  {
    data_->SetDOFHolderData(v);
  }

private:
  class CtorKey
  {
    CtorKey() = default;
    friend class FieldFactory;
  };

  Field(CtorKey, std::shared_ptr<const FieldLayout> layout,
        std::unique_ptr<FieldData<T>> data)
    : layout_(std::move(layout)), data_(std::move(data))
  {
  }

  friend class FieldFactory;

  std::shared_ptr<const FieldLayout> layout_;
  std::unique_ptr<FieldData<T>> data_;
};

using FieldVariant = std::variant<Field<int8_t>, Field<int32_t>, Field<int64_t>,
                                  Field<float>, Field<double>>;

} // namespace pcms

#endif // PCMS_COUPLING_FIELD_H
