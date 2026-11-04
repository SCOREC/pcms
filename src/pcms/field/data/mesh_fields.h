#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H

#include "pcms/field/layout/mesh_fields.h"
#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/field_data.h"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/omega_h_array_utils.h"

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
class MeshFieldsFieldData
  : public FieldData<T>
  , public DOFAccess<T>
{
public:
  // Owned values are the stored data. The MeshFields backend holds the values
  // of every DOF holder, owned and ghost, which its evaluator reads: it
  // shares the owned values and gets the ghost values from SynchronizeGhosts.
  MeshFieldsFieldData(std::shared_ptr<const MeshFieldsAdapterLayout> layout,
                      ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      mesh_field_(MakeMeshFieldBackend<T>(*layout_)),
      owned_data_("meshfields_field_data_device",
                  static_cast<size_t>(layout_->GetNumDofHolder()),
                  static_cast<size_t>(layout_->GetNumComponents()))
  {
    if (!mesh_field_) {
      throw pcms_error(
        "MeshFieldsFieldData does not support this layout/order");
    }
  }

  FieldValueType GetValueType() const override
  {
    return ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

  const DOFAccess<T>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<T>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const T, DeviceMemorySpace> GetValues() const override
  {
    return MakeConstRank2View(owned_data_);
  }

  void SetValues(Rank2View<const T, DeviceMemorySpace> values) override
  {
    CopyDeviceRank2ViewToDeviceView(owned_data_, values);
    if (!layout_->HasGhosts()) {
      SyncBackend(MakeConstRank2View(owned_data_));
      return;
    }
    // Write the owned rows into the backend's values with ghosts, keeping its
    // ghost values.
    auto local = ReadBackend();
    ScatterRows(MakeConstRank2View(owned_data_), layout_->GetOwnedIndices(),
                MakeRank2View(local));
    SyncBackend(MakeConstRank2View(local));
  }

  // Values of every DOF holder, owned and ghost, in mesh order, copied
  // from the backend. Valid until the next call.
  Rank2View<const T, DeviceMemorySpace> GetValuesWithGhosts() const
  {
    values_with_ghosts_ = ReadBackend();
    return MakeConstRank2View(values_with_ghosts_);
  }

  void SynchronizeGhosts() override
  {
    const int nc = layout_->GetNumComponents();
    auto& mesh = layout_->GetMesh();
    const auto nodes_per_dim = layout_->GetNodesPerDim();
    for (int dim = 0; dim <= mesh.dim(); ++dim) {
      if (!nodes_per_dim[dim]) {
        continue;
      }
      // The backend's per-entity data is already dof-major, the order
      // Omega_h's sync_array works in.
      const int width = nodes_per_dim[dim] * nc;
      const size_t len =
        static_cast<size_t>(mesh.nents(dim)) * static_cast<size_t>(width);
      Kokkos::View<T*, DeviceMemorySpace> flat("ghost_sync_flat", len);
      mesh_field_->GetData(Rank1View<T, DeviceMemorySpace>(flat.data(), len),
                           nodes_per_dim[dim], nc, dim);
      SynchronizeOmegaHBlock<T>(mesh, dim, width, flat);
      mesh_field_->SetData(
        Rank1View<const T, DeviceMemorySpace>(flat.data(), len),
        nodes_per_dim[dim], nc, dim);
    }
  }

  std::shared_ptr<MeshFieldBackend<T>> GetMeshFieldBackend() const
  {
    return mesh_field_;
  }

private:
  // The backend's values of every DOF holder with ghosts, in mesh order.
  Kokkos::View<T**, DeviceMemorySpace> ReadBackend() const
  {
    const auto nodes_per_dim = layout_->GetNodesPerDim();
    const auto nc = static_cast<size_t>(layout_->GetNumComponents());
    auto& mesh = layout_->GetMesh();
    Kokkos::View<T**, DeviceMemorySpace> local(
      "meshfields_local_values",
      static_cast<size_t>(layout_->GetNumDofHolderWithGhosts()), nc);
    size_t row_offset = 0;
    for (int i = 0; i <= mesh.dim(); ++i) {
      if (!nodes_per_dim[i]) {
        continue;
      }
      const size_t num_rows = static_cast<size_t>(mesh.nents(i)) *
                              static_cast<size_t>(nodes_per_dim[i]);
      const size_t len = num_rows * nc;
      Kokkos::View<T*, DeviceMemorySpace> flat("read_flat", len);
      mesh_field_->GetData(Rank1View<T, DeviceMemorySpace>(flat.data(), len),
                           nodes_per_dim[i], nc, i);
      Kokkos::parallel_for(
        "ReadBackendReorder",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
          0, static_cast<LO>(num_rows)),
        UnflattenFunctor<T>(flat, MakeRank2View(local),
                            static_cast<LO>(row_offset), nc));
      row_offset += num_rows;
    }
    return local;
  }

  // Serialization boundary: meshfields' SetData consumes flat node-major
  // spans, so the shaped data is explicitly repacked into a flat node-major
  // staging buffer here — the one place this backend handles flat memory.
  void SyncBackend(Rank2View<const T, DeviceMemorySpace> data)
  {
    auto nodes_per_dim = layout_->GetNodesPerDim();
    auto num_components = layout_->GetNumComponents();
    auto& mesh = layout_->GetMesh();
    // The shaped data may be LayoutLeft (component-major), but SetData
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
  ValueBasis basis_;
  std::shared_ptr<MeshFieldBackend<T>> mesh_field_;
  Kokkos::View<T**, DeviceMemorySpace> owned_data_;
  // Backs GetValuesWithGhosts.
  mutable Kokkos::View<T**, DeviceMemorySpace> values_with_ghosts_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_FIELD_DATA_H
