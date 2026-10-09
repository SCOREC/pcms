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
  XGCFieldData(std::shared_ptr<const XGCFieldLayout> layout, ValueBasis basis,
               Rank1View<T, HostMemorySpace> data)
    : layout_(std::move(layout)), basis_(std::move(basis)), data_(data)
  {
    PCMS_ALWAYS_ASSERT(layout_ != nullptr);
    PCMS_ALWAYS_ASSERT(static_cast<LO>(data_.size()) ==
                       layout_->GetFullDataSize());
  }

  // Self-allocating constructor: XGCFunctionSpace::CreateFieldImpl uses this
  // to produce a field with internally-managed storage.
  XGCFieldData(std::shared_ptr<const XGCFieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      owned_data_("xgc_field_data",
                  static_cast<size_t>(layout_->GetFullDataSize())),
      data_(owned_data_.data(), owned_data_.extent(0))
  {
    PCMS_ALWAYS_ASSERT(layout_ != nullptr);
  }

  FieldValueType GetValueType() const override
  {
    return ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

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
    // Explicit node-major serialization into the adapter buffer.
    for (size_t i = 0; i < num_dof; ++i) {
      for (size_t c = 0; c < num_comp; ++c) {
        data_(i * num_comp + c) = values(i, c);
      }
    }
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    EnsureDeviceStaging();
    CopyHostRank2ViewToDeviceView(device_data_, GetDOFHolderDataHost());
    return MakeConstRank2View(device_data_);
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
    EnsureDeviceStaging();
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
    // Explicit node-major serialization into the adapter buffer (the mirror
    // keeps the device layout, so it is indexed, not aliased).
    auto host_mirror = Kokkos::create_mirror_view(device_data_);
    Kokkos::deep_copy(host_mirror, device_data_);
    const auto nc = host_mirror.extent(1);
    for (size_t i = 0; i < host_mirror.extent(0); ++i) {
      for (size_t c = 0; c < nc; ++c) {
        data_(i * nc + c) = host_mirror(i, c);
      }
    }
  }

private:
  void EnsureDeviceStaging() const
  {
    const auto nc = static_cast<size_t>(layout_->GetNumComponents());
    const auto n = data_.size() / nc;
    if (device_data_.extent(0) != n || device_data_.extent(1) != nc) {
      device_data_ = Kokkos::View<T**, DeviceMemorySpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing,
                           "xgc_field_data_device"),
        n, nc);
    }
  }

  std::shared_ptr<const XGCFieldLayout> layout_;
  ValueBasis basis_;
  // owned_data_ is non-empty only when the self-allocating constructor is used.
  Kokkos::View<T*, HostMemorySpace> owned_data_;
  Rank1View<T, HostMemorySpace> data_;
  // Shaped device staging, allocated on first device access.
  mutable Kokkos::View<T**, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_DATA_H
