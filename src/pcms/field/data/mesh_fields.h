#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H

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

namespace
{
// Functor that flattens component-major 2D data (LayoutLeft) into dof-major
// order.
template <typename T>
struct FlattenFunctor
{
  Rank2View<const T, DeviceMemorySpace> data_;
  Kokkos::View<T*, DeviceMemorySpace> flat_;
  LO row_off_;
  size_t nc_;
  FlattenFunctor(Rank2View<const T, DeviceMemorySpace> d,
                 Kokkos::View<T*, DeviceMemorySpace> f, LO ro, size_t nc)
    : data_(d), flat_(f), row_off_(ro), nc_(nc)
  {
  }
  KOKKOS_INLINE_FUNCTION void operator()(LO local) const
  {
    LO global_dof = row_off_ + local;
    for (size_t c = 0; c < nc_; ++c) {
      flat_(local * nc_ + c) = data_(global_dof, c);
    }
  }
};
} // namespace

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
                 static_cast<size_t>(layout_->GetNumOwnedDofHolder()),
                 static_cast<size_t>(layout_->GetNumComponents())),
      device_data_("meshfields_field_data_device",
                   static_cast<size_t>(layout_->GetNumOwnedDofHolder()),
                   static_cast<size_t>(layout_->GetNumComponents()))
  {
    if (!mesh_field_) {
      throw pcms_error(
        "MeshFieldsFieldData does not support this layout/order");
    }
  }

  const FieldMetadata& GetMetadata() const override { return metadata_; }

  Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const override
  {
    DeepCopyMismatchLayouts(host_data_, device_data_);
    return MakeConstRank2View(host_data_);
  }

  void SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->OwnedSize()));
    CopyHostRank2ViewToDeviceView(device_data_, values);
    SyncBackend(GetDOFHolderData());
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const override
  {
    return MakeConstRank2View(device_data_);
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> values) override
  {
    PCMS_ALWAYS_ASSERT(values.size() ==
                       static_cast<size_t>(layout_->OwnedSize()));
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
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
    // device_data_ is rank-2 LayoutLeft (component-major), but SetData
    // expects a flat dof-major span.
    size_t row_offset = 0;
    for (int i = 0; i <= mesh.dim(); ++i) {
      if (nodes_per_dim[i]) {
        size_t num_rows = static_cast<size_t>(mesh.nents(i)) *
                          static_cast<size_t>(nodes_per_dim[i]);
        size_t len = num_rows * static_cast<size_t>(num_components);
        Kokkos::View<T*, DeviceMemorySpace> flat("sync_flat", len);
        Kokkos::parallel_for(
          "SyncBackendReorder",
          Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
            0, static_cast<LO>(num_rows)),
          FlattenFunctor<T>(data, flat, static_cast<LO>(row_offset),
                            num_components));
        Rank1View<const T, DeviceMemorySpace> subspan(flat.data(), len);
        mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
        row_offset += num_rows;
      }
    }
  }

  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  FieldMetadata metadata_;
  std::shared_ptr<MeshFieldBackend<T>> mesh_field_;
  mutable Kokkos::View<T**, HostMemorySpace> host_data_;
  Kokkos::View<T**, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
