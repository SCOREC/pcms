#ifndef PCMS_XGC_FIELD_DATA_H
#define PCMS_XGC_FIELD_DATA_H

#include "pcms/field/layout/xgc.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/mpi_type.h"
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <memory>

namespace pcms
{

template <typename T>
class XGCFieldData
  : public FieldData<T>
  , public DOFAccess<T>
{
public:
  // Externally-managed storage: the caller owns the underlying data buffer.
  // The view must remain valid for the lifetime of this object.
  //
  // XGC keeps a full copy of a plane's data on every rank of the plane, and
  // only the plane's root rank (rank 0 of plane_comm) exchanges it. When
  // plane_comm is given, SetOwnedValues broadcasts the root's values to the
  // other ranks of the plane, so every rank of plane_comm must call it.
  XGCFieldData(std::shared_ptr<const XGCFieldLayout> layout, ValueBasis basis,
               Rank1View<T, HostMemorySpace> data,
               MPI_Comm plane_comm = MPI_COMM_NULL)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      data_(data),
      plane_comm_(plane_comm)
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

  const DOFAccess<T>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<T>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const T, DeviceMemorySpace> GetOwnedValues() const override
  {
    EnsureDeviceStaging();
    const auto nc = layout_->GetNumComponents();
    CopyHostRank2ViewToDeviceView(
      device_data_,
      Rank2View<const T, HostMemorySpace>(
        data_.data_handle(), static_cast<LO>(data_.size()) / nc, nc));
    return MakeConstRank2View(device_data_);
  }

  void SetOwnedValues(Rank2View<const T, DeviceMemorySpace> values) override
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
    if (plane_comm_ != MPI_COMM_NULL) {
      MPI_Bcast(data_.data_handle(), static_cast<int>(data_.size()),
                GetMPIType(T{}), 0, plane_comm_);
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
  MPI_Comm plane_comm_ = MPI_COMM_NULL;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_DATA_H
