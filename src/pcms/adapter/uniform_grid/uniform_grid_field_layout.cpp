#include "uniform_grid_field_layout.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include <memory>

namespace pcms
{

template <unsigned Dim>
UniformGridFieldLayout<Dim>::UniformGridFieldLayout(
  UniformGrid<Dim> grid, int num_components,
  CoordinateSystem coordinate_system, int order)
  : grid_(std::move(grid)),
    num_components_(num_components),
    coordinate_system_(coordinate_system),
    order_(order),
    gids_("gids", GetNumDofHolders()),
    dof_holder_coords_("dof_holder_coords", GetNumDofHolders(), Dim),
    owned_("owned", GetNumDofHolders())
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(order_ == 0 || order_ == 1);

  LO num_dofs = GetNumDofHolders();

  // Initialize global IDs and ownership
  for (LO i = 0; i < num_dofs; ++i) {
    gids_[i] = static_cast<GO>(i);
    owned_[i] = true;
  }

  // Initialize DOF holder coordinates at grid vertices or cell centers.
  Real vertex_spacing[Dim];
  for (unsigned d = 0; d < Dim; ++d) {
    vertex_spacing[d] = grid_.edge_length[d] / grid_.divisions[d];
  }

  if constexpr (Dim == 2) {
    LO dof_idx = 0;
    if (order_ == 1) {
      for (LO j = 0; j <= grid_.divisions[1]; ++j) {
        for (LO i = 0; i <= grid_.divisions[0]; ++i) {
          dof_holder_coords_(dof_idx, 0) = grid_.bot_left[0] + i * vertex_spacing[0];
          dof_holder_coords_(dof_idx, 1) = grid_.bot_left[1] + j * vertex_spacing[1];
          ++dof_idx;
        }
      }
    } else {
      for (LO j = 0; j < grid_.divisions[1]; ++j) {
        for (LO i = 0; i < grid_.divisions[0]; ++i) {
          dof_holder_coords_(dof_idx, 0) = grid_.bot_left[0] + (i + 0.5) * vertex_spacing[0];
          dof_holder_coords_(dof_idx, 1) = grid_.bot_left[1] + (j + 0.5) * vertex_spacing[1];
          ++dof_idx;
        }
      }
    }
  } else if constexpr (Dim == 3) {
    LO dof_idx = 0;
    if (order_ == 1) {
      for (LO k = 0; k <= grid_.divisions[2]; ++k) {
        for (LO j = 0; j <= grid_.divisions[1]; ++j) {
          for (LO i = 0; i <= grid_.divisions[0]; ++i) {
            dof_holder_coords_(dof_idx, 0) = grid_.bot_left[0] + i * vertex_spacing[0];
            dof_holder_coords_(dof_idx, 1) = grid_.bot_left[1] + j * vertex_spacing[1];
            dof_holder_coords_(dof_idx, 2) = grid_.bot_left[2] + k * vertex_spacing[2];
            ++dof_idx;
          }
        }
      }
    } else {
      for (LO k = 0; k < grid_.divisions[2]; ++k) {
        for (LO j = 0; j < grid_.divisions[1]; ++j) {
          for (LO i = 0; i < grid_.divisions[0]; ++i) {
            dof_holder_coords_(dof_idx, 0) = grid_.bot_left[0] + (i + 0.5) * vertex_spacing[0];
            dof_holder_coords_(dof_idx, 1) = grid_.bot_left[1] + (j + 0.5) * vertex_spacing[1];
            dof_holder_coords_(dof_idx, 2) = grid_.bot_left[2] + (k + 0.5) * vertex_spacing[2];
            ++dof_idx;
          }
        }
      }
    }
  }
}

template <unsigned Dim>
int UniformGridFieldLayout<Dim>::GetNumComponents() const
{
  return num_components_;
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::GetNumOwnedDofHolder() const
{
  return GetNumDofHolders();
}

template <unsigned Dim>
GO UniformGridFieldLayout<Dim>::GetNumGlobalDofHolder() const
{
  return GetNumDofHolders();
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
const UniformGrid<Dim>& UniformGridFieldLayout<Dim>::GetGrid() const
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
LO UniformGridFieldLayout<Dim>::GetNumDofHolders() const
{
  return order_ == 0 ? GetNumCells() : GetNumVertices();
}

template <unsigned Dim>
EntOffsetsArray UniformGridFieldLayout<Dim>::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  LO n = GetNumDofHolders();
  int entity_dim = (order_ == 0) ? static_cast<int>(Dim) : 0;
  for (int i = entity_dim + 1; i < ent_offsets_len; ++i) {
    offsets[i] = n;
  }
  return offsets;
}

template <unsigned Dim>
ReversePartitionMap2 UniformGridFieldLayout<Dim>::GetReversePartitionMap(
  const redev::Partition& partition) const
{
  throw std::runtime_error("Unimplemented");
}

template <unsigned Dim>
int UniformGridFieldLayout<Dim>::GetOrder() const
{
  return order_;
}

// Explicit template instantiations
template class UniformGridFieldLayout<2>;
template class UniformGridFieldLayout<3>;

} // namespace pcms
