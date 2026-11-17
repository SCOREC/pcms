#include "pcms/field/layout/uniform_grid.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include <memory>

namespace pcms
{

namespace
{

// Functor for initializing 2D grid coordinates on device
template <unsigned Dim>
struct InitGridCoordsFunctor2D
{
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  UniformGrid<Dim> grid_;
  int order_;
  Real vertex_spacing_0_;
  Real vertex_spacing_1_;
  LO nx_;

  InitGridCoordsFunctor2D(Kokkos::View<Real**, DeviceMemorySpace> coords,
                          const UniformGrid<Dim>& grid, int order, Real vs0,
                          Real vs1, LO nx)
    : coords_(coords),
      grid_(grid),
      order_(order),
      vertex_spacing_0_(vs0),
      vertex_spacing_1_(vs1),
      nx_(nx)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO dof_idx) const
  {
    if (order_ == 1) {
      LO i = dof_idx % nx_;
      LO j = dof_idx / nx_;
      coords_(dof_idx, 0) = grid_.bot_left[0] + i * vertex_spacing_0_;
      coords_(dof_idx, 1) = grid_.bot_left[1] + j * vertex_spacing_1_;
    } else {
      LO i = dof_idx % grid_.divisions[0];
      LO j = dof_idx / grid_.divisions[0];
      coords_(dof_idx, 0) = grid_.bot_left[0] + (i + 0.5) * vertex_spacing_0_;
      coords_(dof_idx, 1) = grid_.bot_left[1] + (j + 0.5) * vertex_spacing_1_;
    }
  }
};

// Functor for initializing 3D grid coordinates on device
template <unsigned Dim>
struct InitGridCoordsFunctor3D
{
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  UniformGrid<Dim> grid_;
  int order_;
  Real vertex_spacing_0_;
  Real vertex_spacing_1_;
  Real vertex_spacing_2_;
  LO nx_;
  LO ny_;

  InitGridCoordsFunctor3D(Kokkos::View<Real**, DeviceMemorySpace> coords,
                          const UniformGrid<Dim>& grid, int order, Real vs0,
                          Real vs1, Real vs2, LO nx, LO ny)
    : coords_(coords),
      grid_(grid),
      order_(order),
      vertex_spacing_0_(vs0),
      vertex_spacing_1_(vs1),
      vertex_spacing_2_(vs2),
      nx_(nx),
      ny_(ny)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO dof_idx) const
  {
    if (order_ == 1) {
      LO i = dof_idx % nx_;
      LO j = (dof_idx / nx_) % ny_;
      LO k = dof_idx / (nx_ * ny_);
      coords_(dof_idx, 0) = grid_.bot_left[0] + i * vertex_spacing_0_;
      coords_(dof_idx, 1) = grid_.bot_left[1] + j * vertex_spacing_1_;
      coords_(dof_idx, 2) = grid_.bot_left[2] + k * vertex_spacing_2_;
    } else {
      LO i = dof_idx % grid_.divisions[0];
      LO j = (dof_idx / grid_.divisions[0]) % grid_.divisions[1];
      LO k = dof_idx / (grid_.divisions[0] * grid_.divisions[1]);
      coords_(dof_idx, 0) = grid_.bot_left[0] + (i + 0.5) * vertex_spacing_0_;
      coords_(dof_idx, 1) = grid_.bot_left[1] + (j + 0.5) * vertex_spacing_1_;
      coords_(dof_idx, 2) = grid_.bot_left[2] + (k + 0.5) * vertex_spacing_2_;
    }
  }
};

struct InitilizeGidsAndOwnedFunctor
{
  Kokkos::View<GO*> gids_;
  Kokkos::View<bool*> owned_;

  InitilizeGidsAndOwnedFunctor(Kokkos::View<GO*> gids,
                               Kokkos::View<bool*> owned)
    : gids_(gids), owned_(owned)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    gids_[i] = static_cast<GO>(i);
    owned_[i] = true;
  }
};

} // anonymous namespace

template <unsigned Dim>
UniformGridFieldLayout<Dim>::UniformGridFieldLayout(
  UniformGrid<Dim> grid, int num_components,
  std::shared_ptr<const CoordinateSystem> coordinate_system, int order)
  : grid_(std::move(grid)),
    num_components_(num_components),
    order_(order),
    gids_("gids", GetNumDofHolders()),
    gids_host_("gids_host", GetNumDofHolders()),
    dof_holder_coords_("dof_holder_coords", GetNumDofHolders(), Dim),
    owned_("owned", GetNumDofHolders()),
    owned_host_("owned_host", GetNumDofHolders())
{
  PCMS_FUNCTION_TIMER;
  SetCoordinateSystem(ResolveCoordinateSystem(std::move(coordinate_system),
                                              static_cast<int>(Dim)));
  PCMS_ALWAYS_ASSERT(order_ == 0 || order_ == 1);

  LO num_dofs = GetNumDofHolders();

  Kokkos::parallel_for(
    "InitUniformGridGidsAndOwned",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_dofs),
    InitilizeGidsAndOwnedFunctor(gids_, owned_));

  Kokkos::deep_copy(gids_host_, gids_);
  Kokkos::deep_copy(owned_host_, owned_);

  // Initialize DOF holder coordinates directly on device using parallel
  // dispatch
  Real vertex_spacing[Dim];
  for (unsigned d = 0; d < Dim; ++d) {
    vertex_spacing[d] = grid_.edge_length[d] / grid_.divisions[d];
  }

  if constexpr (Dim == 2) {
    LO nx = (order_ == 1) ? (grid_.divisions[0] + 1) : grid_.divisions[0];
    InitGridCoordsFunctor2D<Dim> functor(dof_holder_coords_, grid_, order_,
                                         vertex_spacing[0], vertex_spacing[1],
                                         nx);
    Kokkos::parallel_for(
      "InitUniformGrid2DCoords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_dofs),
      functor);
  } else if constexpr (Dim == 3) {
    LO nx = (order_ == 1) ? (grid_.divisions[0] + 1) : grid_.divisions[0];
    LO ny = (order_ == 1) ? (grid_.divisions[1] + 1) : grid_.divisions[1];
    InitGridCoordsFunctor3D<Dim> functor(dof_holder_coords_, grid_, order_,
                                         vertex_spacing[0], vertex_spacing[1],
                                         vertex_spacing[2], nx, ny);
    Kokkos::parallel_for(
      "InitUniformGrid3DCoords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_dofs),
      functor);
  }

  int entity_dim = (order_ == 0) ? static_cast<int>(Dim) : 0;
  LO n = GetNumDofHolders();
  classification_dims_ =
    Kokkos::View<LO*, DeviceMemorySpace>("classification_dims", n);
  classification_ids_ =
    Kokkos::View<LO*, DeviceMemorySpace>("classification_ids", n);
  Kokkos::deep_copy(classification_dims_host_, static_cast<LO>(entity_dim));
  Kokkos::deep_copy(classification_ids_host_, LO{0});

  classification_dims_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_dims", n);
  classification_ids_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_ids", n);
  Kokkos::deep_copy(classification_dims_host_, classification_dims_);
  Kokkos::deep_copy(classification_ids_host_, classification_ids_);
  discretization_ = std::make_shared<UniformGridDiscretization<Dim>>(grid_);
}

template <unsigned Dim>
std::shared_ptr<const Discretization>
UniformGridFieldLayout<Dim>::GetDiscretization() const noexcept
{
  return discretization_;
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
Rank1View<const bool, HostMemorySpace>
UniformGridFieldLayout<Dim>::GetOwnedHost() const
{
  return make_const_array_view(owned_host_);
}

template <unsigned Dim>
GlobalIDView<HostMemorySpace> UniformGridFieldLayout<Dim>::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size());
}

template <unsigned Dim>
CoordinateView<DeviceMemorySpace>
UniformGridFieldLayout<Dim>::GetDOFHolderCoordinates() const
{
  auto coords_view = MakeConstRank2View(dof_holder_coords_);
  return CoordinateView<DeviceMemorySpace>{GetCoordinateSystem(), coords_view};
}

template <unsigned Dim>
bool UniformGridFieldLayout<Dim>::IsDistributed() const
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
int UniformGridFieldLayout<Dim>::GetDimension() const
{
  return static_cast<int>(Dim);
}

template <unsigned Dim>
Rank1View<const LO, HostMemorySpace>
UniformGridFieldLayout<Dim>::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

template <unsigned Dim>
Rank1View<const LO, HostMemorySpace>
UniformGridFieldLayout<Dim>::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
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
int UniformGridFieldLayout<Dim>::GetOrder() const
{
  return order_;
}

// Explicit template instantiations
template class UniformGridFieldLayout<2>;
template class UniformGridFieldLayout<3>;

} // namespace pcms
