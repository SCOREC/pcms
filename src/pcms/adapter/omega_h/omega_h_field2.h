#ifndef PCMS_OMEGA_H_FIELD2_H
#define PCMS_OMEGA_H_FIELD2_H

#include <Kokkos_Core.hpp>
#include <MeshField.hpp>
#include <Omega_h_for.hpp>
#include <memory>

#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "pcms/types.h"
#include "pcms/assert.h"
#include "pcms/profile.h"
#include "pcms/field.h"
#include "pcms/coordinate_system.h"
#include "pcms/point_search.h"

namespace pcms
{
template <typename T>
class MeshFieldBackend
{
public:
  virtual ~MeshFieldBackend() = default;
  virtual Kokkos::View<T* [1]> evaluate(Kokkos::View<Real**> localCoords,
                                        Kokkos::View<LO*> offsets) const = 0;
  virtual void SetData(Rank1View<const T, HostMemorySpace> data,
                       size_t num_nodes, size_t num_components, int dim) = 0;
  virtual void GetData(Rank1View<T, HostMemorySpace> data, size_t num_nodes,
                       size_t num_components, int dim) const = 0;
};

template <typename T, int Dim, int Order>
class MeshFieldBackendImpl : public MeshFieldBackend<T>
{
public:
  MeshFieldBackendImpl(Omega_h::Mesh& mesh)
    : mesh_(mesh),
      mesh_field_(mesh),
      shape_field_(mesh_field_.template CreateLagrangeField<T, Order, 1>())
  {
  }

  Kokkos::View<T* [1]> evaluate(Kokkos::View<Real**> localCoords,
                                Kokkos::View<LO*> offsets) const override {
    auto self = const_cast<MeshFieldBackendImpl<T, Dim, Order>*>(this);
    return self->mesh_field_.triangleLocalPointEval(localCoords, offsets,
                                                    shape_field_);
  }

  void SetData(Rank1View<const T, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Omega_h::parallel_for(
      mesh_.nents(dim), OMEGA_H_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            shape_field_(ent, n, c, topo) =
              data[ent * stride + n * num_components + c];
          }
        }
      });
  }

  void GetData(Rank1View<T, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) const override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Omega_h::parallel_for(
      mesh_.nents(dim), OMEGA_H_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            data[ent * stride + n * num_components + c] =
              shape_field_(ent, n, c, topo);
          }
        }
      });
  }

private:
  Omega_h::Mesh& mesh_;
  MeshField::OmegahMeshField<DefaultExecutionSpace, Dim> mesh_field_;
  using ShapeField =
    decltype(mesh_field_.template CreateLagrangeField<T, Order, 1>());
  ShapeField shape_field_;
};

struct OmegaHField2LocalizationHint
{
  OmegaHField2LocalizationHint(
    Omega_h::Mesh& mesh, Kokkos::View<GridPointSearch::Result*> search_results)
    : offsets_("", mesh.nelems() + 1),
      coordinates_("", search_results.size(), mesh.dim() + 1),
      indices_("", search_results.size())
  {
    Kokkos::View<LO*> elem_counts("", mesh.nelems());

    for (size_t i = 0; i < search_results.size(); ++i) {
      auto [dim, elem_idx, coord] = search_results(i);
      elem_counts[elem_idx] += 1;
    }

    LO total;
    Kokkos::parallel_scan(
      mesh.nelems(),
      KOKKOS_LAMBDA(LO i, LO & partial, bool is_final) {
        if (is_final) {
          offsets_(i) = partial;
        }
        partial += elem_counts(i);
      },
      total);
    offsets_(mesh.nelems()) = total;

    for (size_t i = 0; i < search_results.size(); ++i) {
      auto [dim, elem_idx, coord] = search_results(i);
      // currently don't handle case where point is on a boundary
      PCMS_ALWAYS_ASSERT(static_cast<int>(dim) == mesh.dim());
      // element should be inside the domain (positive)
      PCMS_ALWAYS_ASSERT(elem_idx >= 0 && elem_idx < mesh.nelems());
      elem_counts(elem_idx) -= 1;
      LO index = offsets_(elem_idx) + elem_counts(elem_idx);
      for (int j = 0; j < (mesh.dim() + 1); ++j) {
        coordinates_(index, j) = coord[j];
      }
      // coordinates_(index, mesh.dim()) = coord[0];
      indices_(index) = i;
    }
  }

  // offsets is the number of points in each element
  Kokkos::View<LO*> offsets_;
  // coordinates are the parametric coordinates of each point
  Kokkos::View<Real**> coordinates_;
  // indices are the index of the original point
  Kokkos::View<LO*> indices_;
};

template <typename T>
class OmegaHField2 : public FieldT<T>
{
public:
  OmegaHField2(const OmegaHFieldLayout& layout)
    : layout_(layout),
      mesh_(layout.GetMesh()),
      search_(mesh_, 10, 10),
      dof_holder_data_("", static_cast<size_t>(layout.OwnedSize()))
  {
    auto nodes_per_dim = layout.GetNodesPerDim();
    if (nodes_per_dim[2] == 0 && nodes_per_dim[3] == 0) {
      if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 0) {
        switch (mesh_.dim()) {
          case 1:
            mesh_field_ =
              std::make_unique<MeshFieldBackendImpl<T, 1, 1>>(mesh_);
            break;
          case 2:
            mesh_field_ =
              std::make_unique<MeshFieldBackendImpl<T, 2, 1>>(mesh_);
            break;
          default: break; // backend is null
        }
      } else if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 1) {
        switch (mesh_.dim()) {
          case 2:
            mesh_field_ =
              std::make_unique<MeshFieldBackendImpl<T, 2, 2>>(mesh_);
            break;
          case 3:
            mesh_field_ =
              std::make_unique<MeshFieldBackendImpl<T, 3, 2>>(mesh_);
            break;
          default: break; // backend is null
        }
      }
    }
  }

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) const override
  {
    PCMS_FUNCTION_TIMER;
    // TODO decide if we want to implicitly perform the coordinate
    // transformations when possible
    if (coordinate_view.GetCoordinateSystem() !=
        layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
      // TODO when moved to PCMS throw PCMS exception
      throw std::runtime_error("Coordinate system mismatch");
    }

    auto coordinates = coordinate_view.GetCoordinates();
    Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
    Kokkos::parallel_for(
      coordinates.size() / 2, KOKKOS_LAMBDA(LO i) {
        coords(i, 0) = coordinates(i, 0);
        coords(i, 1) = coordinates(i, 1);
      });
    auto results = search_(coords);
    auto hint = std::make_shared<OmegaHField2LocalizationHint>(mesh_, results);

    return LocalizationHint{hint};
  }

  void Evaluate(LocalizationHint location,
                FieldDataView<Real, HostMemorySpace> results) const override
  {
    PCMS_FUNCTION_TIMER;
    // TODO decide if we want to implicitly perform the coordinate
    // transformations when possible
    if (results.GetCoordinateSystem() !=
        layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
      // TODO when moved to PCMS throw PCMS exception
      throw std::runtime_error("Coordinate system mismatch");
    }

    OmegaHField2LocalizationHint hint =
      *reinterpret_cast<OmegaHField2LocalizationHint*>(location.data.get());

    auto eval_results = mesh_field_->evaluate(hint.coordinates_, hint.offsets_);

    Rank1View<Real, HostMemorySpace> values = results.GetValues();

    Kokkos::parallel_for(
      eval_results.size(),
      KOKKOS_LAMBDA(LO i) { values[hint.indices_(i)] = eval_results(i, 0); });
  }

  void EvaluateGradient(
    FieldDataView<Real, HostMemorySpace> /* unused */) override
  {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Not implemented");
  }

  const FieldLayout& GetLayout() const override { return layout_; }

  bool CanEvaluateGradient() override
  {
    // TODO compute the gradient field using element shape functions
    return false;
  }

  int Serialize(
    Rank1View<Real, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const override
  {
    PCMS_FUNCTION_TIMER;
    // host copy of filtered field data array
    const auto array_h = GetDOFHolderData();
    if (buffer.size() > 0) {
      auto owned = layout_.GetOwned();
      for (size_t i = 0; i < array_h.size(); i++) {
        if (owned[i])
          buffer[permutation[i]] = array_h[i];
      }
    }
    return array_h.size();
  }

  void Deserialize(
    Rank1View<const Real, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) override
  {
    PCMS_FUNCTION_TIMER;
    Omega_h::HostWrite<Real> sorted_buffer(permutation.size());
    auto owned = layout_.GetOwned();
    for (LO i = 0; i < sorted_buffer.size(); ++i) {
      if (owned[i])
        sorted_buffer[i] = buffer[permutation[i]];
    }
    const auto sorted_buffer_d = Omega_h::Read<Real>(sorted_buffer);
    SetDOFHolderData(pcms::make_const_array_view(sorted_buffer_d));
  }

  Rank1View<const Real, HostMemorySpace> GetDOFHolderData() const override
  {
    PCMS_FUNCTION_TIMER;
    auto nodes_per_dim = layout_.GetNodesPerDim();
    auto num_components = layout_.GetNumComponents();
    size_t offset = 0;
    for (int i = 0; i <= mesh_.dim(); ++i) {
      if (nodes_per_dim[i]) {
        size_t len = static_cast<size_t>(mesh_.nents(i)) *
                     static_cast<size_t>(nodes_per_dim[i]) *
                     static_cast<size_t>(num_components);
        Rank1View<Real, HostMemorySpace> subspan{
          std::data(dof_holder_data_) + offset, len};
        mesh_field_->GetData(subspan, nodes_per_dim[i], num_components, i);
        offset += len;
      }
    }

    return make_const_array_view(dof_holder_data_);
  };

  void SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data) override
  {
    PCMS_FUNCTION_TIMER;

    auto nodes_per_dim = layout_.GetNodesPerDim();
    auto num_components = layout_.GetNumComponents();
    PCMS_ALWAYS_ASSERT(static_cast<LO>(data.size()) ==
                       layout_.GetNumOwnedDofHolder() * num_components);
    size_t offset = 0;
    for (int i = 0; i <= mesh_.dim(); ++i) {
      if (nodes_per_dim[i]) {
        size_t len = static_cast<size_t>(mesh_.nents(i)) *
                     static_cast<size_t>(nodes_per_dim[i]) *
                     static_cast<size_t>(num_components);
        Rank1View<const Real, HostMemorySpace> subspan{
          data.data_handle() + offset, len};
        mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
        offset += len;
      }
    }
  }

  ~OmegaHField2() noexcept = default;

private:
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
  std::unique_ptr<MeshFieldBackend<T>> mesh_field_;
  GridPointSearch search_;
  Kokkos::View<T*> dof_holder_data_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
