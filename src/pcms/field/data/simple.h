#ifndef PCMS_SIMPLE_FIELD_DATA_H
#define PCMS_SIMPLE_FIELD_DATA_H

#include "../field_data.h"
#include "../field_layout.h"
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

// SimpleFieldData<T> is a generic concrete FieldData<T> backed by shaped
// [dof_holder][component] Kokkos Views.
//
// Ownership of the layout is shared — the layout is typically held by the
// factory that created this field data object.
template <typename T>
class SimpleFieldData : public FieldData<T>
{
public:
  SimpleFieldData(std::shared_ptr<const FieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      host_data_("simple_field_data",
                 static_cast<size_t>(layout_->GetNumLocalDofHolder()),
                 static_cast<size_t>(layout_->GetNumComponents())),
      device_data_("simple_field_data_device",
                   static_cast<size_t>(layout_->GetNumLocalDofHolder()),
                   static_cast<size_t>(layout_->GetNumComponents()))
  {
  }

  FieldValueType GetValueType() const override
  {
    return ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

  Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    DeepCopyMismatchLayouts(host_data_, device_data_);
    return MakeConstRank2View(host_data_);
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
    CopyHostRank2ViewToDeviceView(device_data_, values);
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    return MakeConstRank2View(device_data_);
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
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
    const int nc = layout_->GetNumComponents();
    const LO num_rows = static_cast<LO>(device_data_.extent(0));

    // device_data_ may be component-major on device; flatten to dof-major
    // order for Omega_h's sync_array, then scatter the synchronized values
    // back.
    Kokkos::View<T*, DeviceMemorySpace> flat(
      "ghost_sync_flat",
      static_cast<size_t>(num_rows) * static_cast<size_t>(nc));
    CopyDeviceRank2ViewToDeviceView(flat, MakeConstRank2View(device_data_));

    SynchronizeOmegaHBlock<T>(oh->GetMesh(), dim, nc, flat);

    CopyDeviceViewToRank2View<T>(device_data_, flat);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  ValueBasis basis_;
  mutable Kokkos::View<T**, HostMemorySpace> host_data_;
  mutable Kokkos::View<T**, HostMemorySpace> owned_host_data_;
  mutable Kokkos::View<T**, DeviceMemorySpace> owned_device_data_;
  Kokkos::View<T**, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_SIMPLE_FIELD_DATA_H
