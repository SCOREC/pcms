#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_FIELD_DATA_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_FIELD_DATA_H

#include "pcms/field/layout/mesh_fields.h"
#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"

#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

// MeshFieldsReducedQuinticFieldData<T> is a near-duplicate of
// MeshFieldsFieldData<T> that creates a reduced quintic backend instead of
// a standard Lagrange backend.  The two classes are kept separate to
// preserve the original API of MeshFieldsFieldData.
template <typename T>
class MeshFieldsReducedQuinticFieldData : public FieldData<T>
{
public:
  MeshFieldsReducedQuinticFieldData(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout,
    FieldMetadata metadata)
    : layout_(std::move(layout)),
      metadata_(metadata),
      mesh_field_(MakeMeshFieldReducedQuinticBackend<T>(*layout_)),
      host_data_("meshfields_reduced_quintic_field_data",
                 static_cast<size_t>(layout_->OwnedSize())),
      device_data_("meshfields_reduced_quintic_field_data_device",
                   static_cast<size_t>(layout_->OwnedSize()))
  {
    if (!mesh_field_) {
      throw pcms_error(
        "MeshFieldsReducedQuinticFieldData does not support this layout/order");
    }
  }

  const FieldMetadata& GetMetadata() const override { return metadata_; }

  Rank1View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    Kokkos::deep_copy(host_data_, device_data_);
    return make_const_array_view(host_data_);
  }

  void SetDOFHolderDataHost(Rank1View<const T, HostMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->OwnedSize()));
    CopyHostRank1ViewToDeviceView(device_data_, values);
    SyncBackend(make_const_array_view(device_data_));
  }

  Rank1View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    return make_const_array_view(device_data_);
  }

  void SetDOFHolderData(Rank1View<const T, DeviceMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->OwnedSize()));
    CopyDeviceRank1ViewToDeviceView(device_data_, values);
    SyncBackend(make_const_array_view(device_data_));
  }

  std::shared_ptr<MeshFieldBackend<T>> GetMeshFieldBackend() const
  {
    return mesh_field_;
  }

private:
  void SyncBackend(Rank1View<const T, DeviceMemorySpace> flat)
  {
    auto nodes_per_dim = layout_->GetNodesPerDim();
    auto num_components = layout_->GetNumComponents();
    auto& mesh = layout_->GetMesh();
    size_t offset = 0;
    for (int i = 0; i <= mesh.dim(); ++i) {
      if (nodes_per_dim[i]) {
        size_t len = static_cast<size_t>(mesh.nents(i)) *
                     static_cast<size_t>(nodes_per_dim[i]) *
                     static_cast<size_t>(num_components);
        Rank1View<const T, DeviceMemorySpace> subspan{
          flat.data_handle() + offset, len};
        mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
        offset += len;
      }
    }
  }

  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  FieldMetadata metadata_;
  std::shared_ptr<MeshFieldBackend<T>> mesh_field_;
  mutable Kokkos::View<T*, HostMemorySpace> host_data_;
  Kokkos::View<T*, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_FIELD_DATA_H