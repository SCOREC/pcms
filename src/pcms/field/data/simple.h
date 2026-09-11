#ifndef PCMS_SIMPLE_FIELD_DATA_H
#define PCMS_SIMPLE_FIELD_DATA_H

#include "../field_data.h"
#include "../field_layout.h"
#include "../field_metadata.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <Omega_h_mesh.hpp>
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/omega_h_array_utils.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>
#include <type_traits>

namespace pcms
{

// SimpleFieldData<T> is a generic concrete FieldData<T> backed by a flat
// Kokkos::View<T*, HostMemorySpace>. It works for any backend whose DOF data
// is a flat coefficient array (OmegaH, UniformGrid, PointCloud, etc.).
//
// Ownership of the layout is shared — the layout is typically held by the
// factory that created this field data object.
template <typename T>
class SimpleFieldData : public FieldData<T>
{
public:
  SimpleFieldData(std::shared_ptr<const FieldLayout> layout,
                  FieldMetadata metadata)
    : layout_(std::move(layout)),
      metadata_(metadata),
      host_data_("simple_field_data",
                 static_cast<size_t>(layout_->LocalSize())),
      device_data_("simple_field_data_device",
                   static_cast<size_t>(layout_->LocalSize()))
  {
  }

  const FieldMetadata& GetMetadata() const override { return metadata_; }

  Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    Kokkos::deep_copy(host_data_, device_data_);
    return Rank2View<const T, HostMemorySpace>(host_data_.data(),
                                               layout_->GetNumLocalDofHolder(),
                                               layout_->GetNumComponents());
  }

  Rank2View<const T, HostMemorySpace> GetOwnedDOFHolderDataHost() const override
  {
    return GatherOwnedHostData(*layout_, device_data_, host_data_,
                               owned_host_data_);
  }

  Rank2View<const T, DeviceMemorySpace> GetOwnedDOFHolderData() const override
  {
    return GatherOwnedDeviceData(*layout_, device_data_, owned_device_data_);
  }

  void SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->LocalSize()));
    CopyHostRank2ViewToDeviceView(device_data_, values);
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    return Rank2View<const T, DeviceMemorySpace>(
      device_data_.data(), layout_->GetNumLocalDofHolder(),
      layout_->GetNumComponents());
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->LocalSize()));
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
  }

  void SynchronizeGhosts() override
  {
    auto* oh = dynamic_cast<const OmegaHDiscretization*>(
      layout_->GetDiscretization().get());
    if (!oh) {
      return;
    }
    const int dim = layout_->GetDOFHolderEntityDim();
    if (dim < 0) {
      return;
    }
    SynchronizeOmegaHBlock<T>(oh->GetMesh(), dim, layout_->GetNumComponents(),
                              device_data_);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  FieldMetadata metadata_;
  mutable Kokkos::View<T*, HostMemorySpace> host_data_;
  mutable Kokkos::View<T*, HostMemorySpace> owned_host_data_;
  mutable Kokkos::View<T*, DeviceMemorySpace> owned_device_data_;
  Kokkos::View<T*, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_SIMPLE_FIELD_DATA_H
