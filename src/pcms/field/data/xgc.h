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

  // Self-allocating constructor: XGCFunctionSpace::CreateFieldImpl uses this
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

  Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    const auto nc = layout_->GetNumComponents();
    return Rank2View<const T, HostMemorySpace>(
      data_.data_handle(), static_cast<LO>(data_.size()) / nc, nc);
  }

  void SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace> values) override
  {
    if (values.size() != data_.size()) {
      throw pcms_error("XGCFieldData::SetDOFHolderDataHost: size mismatch");
    }
    const auto num_dof = values.extent(0);
    const auto num_comp = values.extent(1);
    for (size_t i = 0; i < num_dof; ++i) {
      for (size_t c = 0; c < num_comp; ++c) {
        data_(i * num_comp + c) = values(i, c);
      }
    }
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    Kokkos::View<T*, HostMemorySpace> host_view("GetDOFHolderData_host_view",
                                                data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
      host_view(i) = data_(i);
    }
    Kokkos::deep_copy(device_data_, host_view);
    const auto nc = layout_->GetNumComponents();
    return Rank2View<const T, DeviceMemorySpace>(
      device_data_.data(), static_cast<LO>(device_data_.extent(0)) / nc, nc);
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
    // Fill the node-major device scratch from the (possibly non-node-major)
    // 2D view, then mirror it back into the canonical host storage.
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
    CopyRank1ViewToHost(data_, make_const_array_view(device_data_));
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
