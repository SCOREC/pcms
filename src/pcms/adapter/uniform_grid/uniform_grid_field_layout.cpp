#include "uniform_grid_field_layout.h"
#include "uniform_grid_field.h"
#include "pcms/utility/profile.h"
#include <memory>

namespace pcms
{

template <unsigned Dim>
UniformGridFieldLayout<Dim>::UniformGridFieldLayout(
  UniformGrid<Dim>& grid, int num_components,
  CoordinateSystem coordinate_system)
  : grid_(grid),
    num_components_(num_components),
    coordinate_system_(coordinate_system),
    gids_("gids", GetNumVertices()),
    dof_holder_coords_("dof_holder_coords", GetNumVertices(), Dim),
    owned_("owned", GetNumVertices())
{
  PCMS_FUNCTION_TIMER;

  LO num_vertices = GetNumVertices();

  // Initialize global IDs and ownership
  for (LO i = 0; i < num_vertices; ++i) {
    gids_[i] = static_cast<GO>(i);
    owned_[i] = true;
  }

  // Initialize DOF holder coordinates at grid vertices
  Real vertex_spacing[Dim];
  for (unsigned d = 0; d < Dim; ++d) {
    vertex_spacing[d] = grid_.edge_length[d] / grid_.divisions[d];
  }

  if constexpr (Dim == 2) {
    LO vertex_idx = 0;
    for (LO j = 0; j <= grid_.divisions[1]; ++j) {
      for (LO i = 0; i <= grid_.divisions[0]; ++i) {
        dof_holder_coords_(vertex_idx, 0) =
          grid_.bot_left[0] + i * vertex_spacing[0];
        dof_holder_coords_(vertex_idx, 1) =
          grid_.bot_left[1] + j * vertex_spacing[1];
        ++vertex_idx;
      }
    }
  } else if constexpr (Dim == 3) {
    LO vertex_idx = 0;
    for (LO k = 0; k <= grid_.divisions[2]; ++k) {
      for (LO j = 0; j <= grid_.divisions[1]; ++j) {
        for (LO i = 0; i <= grid_.divisions[0]; ++i) {
          dof_holder_coords_(vertex_idx, 0) =
            grid_.bot_left[0] + i * vertex_spacing[0];
          dof_holder_coords_(vertex_idx, 1) =
            grid_.bot_left[1] + j * vertex_spacing[1];
          dof_holder_coords_(vertex_idx, 2) =
            grid_.bot_left[2] + k * vertex_spacing[2];
          ++vertex_idx;
        }
      }
    }
  }
}

template <unsigned Dim>
std::unique_ptr<FieldT<Real>> UniformGridFieldLayout<Dim>::CreateField() const
{
  return std::make_unique<UniformGridField<Dim>>(*this);
}

template <unsigned Dim>
int UniformGridFieldLayout<Dim>::GetNumComponents() const
{
  return num_components_;
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::GetNumOwnedDofHolder() const
{
  return GetNumVertices();
}

template <unsigned Dim>
GO UniformGridFieldLayout<Dim>::GetNumGlobalDofHolder() const
{
  return GetNumVertices();
}

template <unsigned Dim>
Rank1View<const bool, HostMemorySpace> UniformGridFieldLayout<Dim>::GetOwned()
  const
{
  return make_const_array_view(owned_);
}

template <unsigned Dim>
GlobalIDView<HostMemorySpace> UniformGridFieldLayout<Dim>::GetGids() const
{
  return GlobalIDView<HostMemorySpace>(gids_.data(), gids_.size());
}

template <unsigned Dim>
CoordinateView<HostMemorySpace>
UniformGridFieldLayout<Dim>::GetDOFHolderCoordinates() const
{
  Rank2View<const Real, HostMemorySpace> coords_view(
    dof_holder_coords_.data(), dof_holder_coords_.extent(0), Dim);
  return CoordinateView<HostMemorySpace>{coordinate_system_, coords_view};
}

template <unsigned Dim>
bool UniformGridFieldLayout<Dim>::IsDistributed()
{
  return false;
}

template <unsigned Dim>
UniformGrid<Dim>& UniformGridFieldLayout<Dim>::GetGrid() const
{
  return grid_;
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::GetNumCells() const
{
  return grid_.GetNumCells();
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::GetNumVertices() const
{
  LO num_vertices = 1;
  for (unsigned d = 0; d < Dim; ++d) {
    num_vertices *= (grid_.divisions[d] + 1);
  }
  return num_vertices;
}

template <unsigned Dim>
EntOffsetsArray UniformGridFieldLayout<Dim>::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  offsets[0] = 0;
  offsets[1] = grid_.GetNumCells();
  offsets[2] = grid_.GetNumCells();
  offsets[3] = grid_.GetNumCells();
  offsets[4] = grid_.GetNumCells();
  return offsets;
}

template <unsigned Dim>
ReversePartitionMap2 UniformGridFieldLayout<Dim>::GetReversePartitionMap(
  const redev::Partition& partition) const
{
  throw std::runtime_error("Unimplemented");
}

// Explicit template instantiations
template class UniformGridFieldLayout<2>;
template class UniformGridFieldLayout<3>;

} // namespace pcms
