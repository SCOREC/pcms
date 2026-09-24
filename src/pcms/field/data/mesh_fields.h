#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H

#include "pcms/field/layout/mesh_fields.h"
#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/omega_h_array_utils.h"

#include <Kokkos_Core.hpp>
#include <memory>
#include <type_traits>

namespace pcms
{

template <typename T>
class MeshFieldsFieldData : public FieldData<T>
{
public:
  MeshFieldsFieldData(std::shared_ptr<const MeshFieldsAdapterLayout> layout,
                      FieldMetadata metadata)
    : layout_(std::move(layout)),
      metadata_(metadata),
      mesh_field_(MakeMeshFieldBackend<T>(*layout_)),
      host_data_("meshfields_field_data",
                 static_cast<size_t>(layout_->LocalSize())),
      device_data_("meshfields_field_data_device",
                   static_cast<size_t>(layout_->LocalSize()))
  {
    if (!mesh_field_) {
      throw pcms_error(
        "MeshFieldsFieldData does not support this layout/order");
    }
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
    SyncBackend(GetDOFHolderData());
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    // The Rank2View will wrap the dof-major data with layout left when device
    // memory is enabled. This may cause issues in multi component cases. See
    // issue #342
    return Rank2View<const T, DeviceMemorySpace>(
      device_data_.data(), layout_->GetNumLocalDofHolder(),
      layout_->GetNumComponents());
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->LocalSize()));
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
    SyncBackend(GetDOFHolderData());
  }

  void SynchronizeGhosts() override
  {
    const int nc = layout_->GetNumComponents();
    auto& mesh = layout_->GetMesh();
    const auto nodes_per_dim = layout_->GetNodesPerDim();

    size_t row_offset = 0;
    for (int dim = 0; dim <= mesh.dim(); ++dim) {
      if (!nodes_per_dim[dim]) {
        continue;
      }
      const LO num_rows = static_cast<LO>(mesh.nents(dim)) * nodes_per_dim[dim];
      const LO flat_len = num_rows * nc;
      const LO flat_off = static_cast<LO>(row_offset * static_cast<size_t>(nc));

      auto block = Kokkos::subview(
        device_data_, Kokkos::make_pair(flat_off, flat_off + flat_len));
      SynchronizeOmegaHBlock<T>(mesh, dim, nc, block);

      row_offset += static_cast<size_t>(num_rows);
    }

    SyncBackend(GetDOFHolderData());
  }

  std::shared_ptr<MeshFieldBackend<T>> GetMeshFieldBackend() const
  {
    return mesh_field_;
  }

private:
  void SyncBackend(Rank2View<const T, DeviceMemorySpace> data)
  {
    auto nodes_per_dim = layout_->GetNodesPerDim();
    auto num_components = layout_->GetNumComponents();
    auto& mesh = layout_->GetMesh();
    // data is [dof_holder][component], contiguous node-major, so each mesh
    // dimension owns a contiguous block of rows; SetData consumes a flat
    // node-major span over that block.
    size_t row_offset = 0;
    for (int i = 0; i <= mesh.dim(); ++i) {
      if (nodes_per_dim[i]) {
        size_t num_rows = static_cast<size_t>(mesh.nents(i)) *
                          static_cast<size_t>(nodes_per_dim[i]);
        size_t len = num_rows * static_cast<size_t>(num_components);
        Rank1View<const T, DeviceMemorySpace> subspan{
          data.data_handle() + row_offset * static_cast<size_t>(num_components),
          len};
        mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
        row_offset += num_rows;
      }
    }
  }

  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  FieldMetadata metadata_;
  std::shared_ptr<MeshFieldBackend<T>> mesh_field_;
  mutable Kokkos::View<T*, HostMemorySpace> host_data_;
  mutable Kokkos::View<T*, HostMemorySpace> owned_host_data_;
  mutable Kokkos::View<T*, DeviceMemorySpace> owned_device_data_;
  Kokkos::View<T*, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
