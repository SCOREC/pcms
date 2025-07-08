#include "omega_h_field2.h"
#include "omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <memory>

namespace pcms
{
template <int Dim, int Order>
class MeshFieldBackendImpl : public MeshFieldBackend
{
public:
  MeshFieldBackendImpl(Omega_h::Mesh& mesh)
    : mesh_(mesh),
      mesh_field_(mesh),
      shape_field_(mesh_field_.template CreateLagrangeField<Real, Order>())
  {
  }

  Kokkos::View<Real* [1]> evaluate(Kokkos::View<Real**> localCoords,
                                   Kokkos::View<LO*> offsets) const override {
    auto self = const_cast<MeshFieldBackendImpl<Dim, Order>*>(this);
    return self->mesh_field_.triangleLocalPointEval(localCoords, offsets,
                                                    shape_field_);
  }

  void SetData(Rank1View<const Real, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Omega_h::parallel_for(
      mesh_.nents(dim), OMEGA_H_LAMBDA(size_t ent) {
        for (int n = 0; n < num_nodes; ++n) {
          for (int c = 0; c < num_components; ++c) {
            shape_field_(ent, n, c, topo) =
              data[ent * stride + n * num_components + c];
          }
        }
      });
  }

  void GetData(Rank1View<Real, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) const override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Omega_h::parallel_for(
      mesh_.nents(dim), OMEGA_H_LAMBDA(size_t ent) {
        for (int n = 0; n < num_nodes; ++n) {
          for (int c = 0; c < num_components; ++c) {
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
    decltype(mesh_field_.template CreateLagrangeField<Real, Order>());
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

    for (LO i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        elem_counts[elem_idx] += 1;
    }

    LO total;
    Kokkos::parallel_scan(
      mesh.nelems(),
      KOKKOS_LAMBDA(LO i, LO & partial, bool is_final) {
        if (is_final)
          offsets_(i) = partial;
        partial += elem_counts(i);
      }, total);
    offsets_(mesh.nelems()) = total;

    for (LO i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        elem_counts(elem_idx) -= 1;
        LO index = offsets_(elem_idx) + elem_counts(elem_idx);
        for (int j = 1; j <= mesh.dim(); ++j)
          coordinates_(index, j - 1) = coord[j];
        coordinates_(index, mesh.dim()) = coord[0];
        indices_(index) = i;
    }
  }

  Kokkos::View<LO*> offsets_;
  Kokkos::View<Real**> coordinates_;
  Kokkos::View<LO*> indices_;
};

/*
 * Field
 */
OmegaHField2::OmegaHField2(std::string name, CoordinateSystem coordinate_system,
                           const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh)
  : name_(name),
    coordinate_system_(coordinate_system),
    layout_(layout),
    mesh_(mesh),
    search_(mesh, 10, 10),
    dof_holder_data_("", layout.GetNumOwnedDofHolder() * layout.GetNumComponents())
{
  auto nodes_per_dim = layout.GetNodesPerDim();
  if (nodes_per_dim[2] == 0 && nodes_per_dim[3] == 0) {
    if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 0) {
        switch (mesh.dim()) {
          case 1:
            mesh_field_ = std::make_unique<MeshFieldBackendImpl<1, 1>>(mesh);
            break;
          case 2:
            mesh_field_ = std::make_unique<MeshFieldBackendImpl<2, 1>>(mesh);
            break;
        }
    } else if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 1) {
        switch (mesh.dim()) {
          case 2:
            mesh_field_ = std::make_unique<MeshFieldBackendImpl<2, 2>>(mesh);
            break;
          case 3:
            mesh_field_ = std::make_unique<MeshFieldBackendImpl<3, 2>>(mesh);
            break;
        }
    }
  }
}

const std::string& OmegaHField2::GetName() const
{
  return name_;
}

CoordinateSystem OmegaHField2::GetCoordinateSystem() const
{
  return coordinate_system_;
}

FieldDataView<const Real, HostMemorySpace> OmegaHField2::GetDOFHolderData()
  const
{
  PCMS_FUNCTION_TIMER;
  auto nodes_per_dim = layout_.GetNodesPerDim();
  auto num_components = layout_.GetNumComponents();
  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim[i]) {
      size_t len = mesh_.nents(i) * nodes_per_dim[i] * num_components;
      Rank1View<Real, HostMemorySpace> subspan{std::data(dof_holder_data_) + offset, len};
      mesh_field_->GetData(subspan, nodes_per_dim[i], num_components, i);
      offset += len;
    }
  }

  Rank1View<const Real, pcms::HostMemorySpace> array_view{
    std::data(dof_holder_data_), std::size(dof_holder_data_)};
  FieldDataView<const Real, HostMemorySpace> data_view{array_view,
                                                       GetCoordinateSystem()};
  return data_view;
};

CoordinateView<HostMemorySpace> OmegaHField2::GetDOFHolderCoordinates() const
{
  PCMS_FUNCTION_TIMER;
  return CoordinateView<HostMemorySpace>{GetCoordinateSystem(),
                                         layout_.GetDOFHolderCoordinates()};
}

void OmegaHField2::SetDOFHolderData(FieldDataView<const Real, HostMemorySpace> data) {
  PCMS_FUNCTION_TIMER;
  if (data.GetCoordinateSystem() != coordinate_system_) {
    throw std::runtime_error("Coordinate system mismatch");
  }

  Rank1View<const Real, HostMemorySpace> values = data.GetValues();
  auto nodes_per_dim = layout_.GetNodesPerDim();
  auto num_components = layout_.GetNumComponents();
  PCMS_ALWAYS_ASSERT(static_cast<LO>(values.size()) ==
                     layout_.GetNumOwnedDofHolder() * num_components);
  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim[i]) {
      size_t len = mesh_.nents(i) * nodes_per_dim[i] * num_components;
      Rank1View<const Real, HostMemorySpace> subspan{values.data_handle() + offset, len};
      mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
      offset += len;
    }
  }
}

LocalizationHint OmegaHField2::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (coordinate_view.GetCoordinateSystem() != coordinate_system_) {
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

void OmegaHField2::Evaluate(LocalizationHint location,
                            FieldDataView<double, HostMemorySpace> results) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (results.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  OmegaHField2LocalizationHint hint =
    *((OmegaHField2LocalizationHint*)location.data.get());

  auto eval_results = mesh_field_->evaluate(hint.coordinates_, hint.offsets_);

  Rank1View<double, HostMemorySpace> values = results.GetValues();

  Kokkos::parallel_for(
    eval_results.size(),
    KOKKOS_LAMBDA(LO i) { values[hint.indices_(i)] = eval_results(i, 0); });
}

void OmegaHField2::EvaluateGradient(FieldDataView<double, HostMemorySpace> results)
{
  // TODO when moved to PCMS throw PCMS exception
  throw std::runtime_error("Not implemented");
}

const FieldLayout& OmegaHField2::GetLayout() const
{
  return layout_;
}

bool OmegaHField2::CanEvaluateGradient()
{
  // TODO compute the gradient field using element shape functions
  return false;
}

int OmegaHField2::Serialize(
  Rank1View<double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;
  // host copy of filtered field data array
  const auto array_h = GetDOFHolderData().GetValues();
  if (buffer.size() > 0) {
    for (LO i = 0; i < array_h.size(); i++) {
      buffer[i] = array_h[permutation[i]];
    }
  }
  return array_h.size();
}

void OmegaHField2::Deserialize(
  Rank1View<const double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
  Omega_h::HostWrite<Real> sorted_buffer(buffer.size());
  for (size_t i = 0; i < buffer.size(); ++i) {
    sorted_buffer[permutation[i]] = buffer[i];
  }
  const auto sorted_buffer_d = Omega_h::Read<Real>(sorted_buffer);
  Rank1View<const Real, HostMemorySpace> sorted_buffer_view{std::data(sorted_buffer_d), std::size(sorted_buffer_d)};
  FieldDataView<const Real, HostMemorySpace> field_view{sorted_buffer_view, GetCoordinateSystem()};
  SetDOFHolderData(field_view);
}
} // namespace pcms
