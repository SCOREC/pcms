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

// Functor that scatters a flat dof-major buffer back into component-major 2D
// data (LayoutLeft). Inverse of FlattenFunctor; indexes the full destination
// view directly (rather than via Kokkos::subview) because a partial row-range
// subview of a LayoutLeft view is not representable as LayoutLeft/LayoutRight
// (Kokkos falls back to LayoutStride), which Rank2View does not detect.
template <typename T>
struct UnflattenFunctor
{
  Kokkos::View<const T*, DeviceMemorySpace> flat_;
  Rank2View<T, DeviceMemorySpace> data_;
  LO row_off_;
  size_t nc_;
  UnflattenFunctor(Kokkos::View<const T*, DeviceMemorySpace> f,
                   Rank2View<T, DeviceMemorySpace> d, LO ro, size_t nc)
    : flat_(f), data_(d), row_off_(ro), nc_(nc)
  {
  }
  KOKKOS_INLINE_FUNCTION void operator()(LO local) const
  {
    LO global_dof = row_off_ + local;
    for (size_t c = 0; c < nc_; ++c) {
      data_(global_dof, c) = flat_(local * nc_ + c);
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
                 static_cast<size_t>(layout_->GetNumLocalDofHolder()),
                 static_cast<size_t>(layout_->GetNumComponents())),
      device_data_("meshfields_field_data_device",
                   static_cast<size_t>(layout_->GetNumLocalDofHolder()),
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
    return MakeConstRank2View(device_data_);
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
      const LO row_off = static_cast<LO>(row_offset);

      // device_data_ may be component-major on device; flatten this block to
      // dof-major order for Omega_h's sync_array, then scatter the
      // synchronized values back. We index the full device_data_ view
      // directly (not a Kokkos::subview of it) because a partial row-range
      // subview of a LayoutLeft view is not representable as
      // LayoutLeft/LayoutRight (Kokkos falls back to LayoutStride), which
      // MakeConstRank2View/MakeRank2View do not detect.
      Kokkos::View<T*, DeviceMemorySpace> flat(
        "ghost_sync_flat", static_cast<size_t>(num_rows) * static_cast<size_t>(nc));
      Kokkos::parallel_for(
        "SynchronizeGhostsFlatten",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_rows),
        FlattenFunctor<T>(MakeConstRank2View(device_data_), flat, row_off,
                          static_cast<size_t>(nc)));

      SynchronizeOmegaHBlock<T>(mesh, dim, nc, flat);

      Kokkos::parallel_for(
        "SynchronizeGhostsUnflatten",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_rows),
        UnflattenFunctor<T>(flat, MakeRank2View(device_data_), row_off,
                           static_cast<size_t>(nc)));

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
  mutable Kokkos::View<T**, HostMemorySpace> owned_host_data_;
  mutable Kokkos::View<T**, DeviceMemorySpace> owned_device_data_;
  Kokkos::View<T**, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
