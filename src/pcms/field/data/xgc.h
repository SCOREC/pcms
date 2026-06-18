#ifndef PCMS_XGC_FIELD_DATA_H
#define PCMS_XGC_FIELD_DATA_H

#include "pcms/field/layout/xgc.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/assert.h"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

template <typename T>
class XGCFieldData : public FieldData<T>
{
public:
  // Externally-managed storage: the caller owns the underlying data buffer.
  // The view must remain valid for the lifetime of this object.
  XGCFieldData(std::shared_ptr<const XGCFieldLayout> layout,
               FieldMetadata metadata, Rank1View<T, HostMemorySpace> data)
    : layout_(std::move(layout)), metadata_(metadata), data_(data)
  {
    PCMS_ALWAYS_ASSERT(layout_ != nullptr);
    PCMS_ALWAYS_ASSERT(static_cast<LO>(data_.size()) ==
                       layout_->GetFullDataSize());
  }

  // Self-allocating constructor: XGCFieldFactory::CreateFieldImpl uses this
  // to produce a field with internally-managed storage.
  XGCFieldData(std::shared_ptr<const XGCFieldLayout> layout,
               FieldMetadata metadata)
    : layout_(std::move(layout)),
      metadata_(metadata),
      owned_data_("xgc_field_data",
                  static_cast<size_t>(layout_->GetFullDataSize())),
      data_(owned_data_.data(), owned_data_.extent(0))
  {
    PCMS_ALWAYS_ASSERT(layout_ != nullptr);
  }

  const FieldMetadata& GetMetadata() const override { return metadata_; }

  Rank1View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    return Rank1View<const T, HostMemorySpace>(data_.data_handle(),
                                               data_.size());
  }

  void SetDOFHolderDataHost(Rank1View<const T, HostMemorySpace> values) override
  {
    if (values.size() != data_.size()) {
      throw pcms_error("XGCFieldData::SetDOFHolderDataHost: size mismatch");
    }
    for (size_t i = 0; i < values.size(); ++i) {
      data_(i) = values[i];
    }
  }

  Rank1View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    Kokkos::View<T*, HostMemorySpace> host_view("GetDOFHolderData_host_view",
                                                data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
      host_view(i) = data_(i);
    }
    Kokkos::deep_copy(device_data_, host_view);
    return make_const_array_view(device_data_);
  }

  void SetDOFHolderData(Rank1View<const T, DeviceMemorySpace> values) override
  {
    CopyDeviceRank1ViewToDeviceView(device_data_, values);
    CopyRank1ViewToHost(data_, values);
  }

private:
  std::shared_ptr<const XGCFieldLayout> layout_;
  FieldMetadata metadata_;
  // owned_data_ is non-empty only when the self-allocating constructor is used.
  Kokkos::View<T*, HostMemorySpace> owned_data_;
  Rank1View<T, HostMemorySpace> data_;
  mutable Kokkos::View<T*, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_DATA_H
