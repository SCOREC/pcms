#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/assert.hpp"
#include <memory>

namespace pcms::uniform_grid
{

namespace
{

// DOF-holder coordinates: order 1 the vertices, order 0 the cell centers
template <unsigned Dim>
struct InitGridCoordsFunctor
{
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  UniformGrid<Dim> grid_;
  int order_;
  unsigned wrapped_;

  KOKKOS_INLINE_FUNCTION
  void operator()(LO dof_idx) const
  {
    LO rest = dof_idx;
    for (unsigned d = 0; d < Dim; ++d) {
      const LO n =
        grid_.divisions[d] + (order_ == 1 && !((wrapped_ >> d) & 1U) ? 1 : 0);
      const Real spacing = grid_.edge_length[d] / grid_.divisions[d];
      coords_(dof_idx, d) =
        grid_.bot_left[d] + (rest % n + (order_ == 1 ? 0.0 : 0.5)) * spacing;
      rest /= n;
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

namespace
{
template <unsigned Dim>
const UniformGrid<Dim>& GridOf(
  const std::shared_ptr<const UniformGridGeometry<Dim>>& geometry)
{
  if (geometry == nullptr) {
    throw pcms_error("UniformGridFieldLayout: the geometry must not be null");
  }
  return geometry->GetGrid();
}
} // anonymous namespace

template <unsigned Dim>
UniformGridFieldLayout<Dim>::UniformGridFieldLayout(
  std::shared_ptr<const UniformGridGeometry<Dim>> geometry, int num_components,
  int order, std::string name)
  : FieldLayout(std::move(name)),
    grid_(GridOf(geometry)),
    num_components_(num_components),
    order_(order),
    wrapped_(geometry->WrappedAxes()),
    gids_("gids", CountDofHolders()),
    gids_host_("gids_host", CountDofHolders()),
    dof_holder_coords_("dof_holder_coords", CountDofHolders(), Dim),
    owned_("owned", CountDofHolders()),
    owned_host_("owned_host", CountDofHolders())
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(order_ == 0 || order_ == 1);

  LO num_dofs = CountDofHolders();

  Kokkos::parallel_for(
    "InitUniformGridGidsAndOwned",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_dofs),
    InitilizeGidsAndOwnedFunctor(gids_, owned_));

  Kokkos::deep_copy(gids_host_, gids_);
  Kokkos::deep_copy(owned_host_, owned_);

  Kokkos::parallel_for(
    "InitUniformGridCoords",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_dofs),
    InitGridCoordsFunctor<Dim>{dof_holder_coords_, grid_, order_, wrapped_});

  int entity_dim = (order_ == 0) ? static_cast<int>(Dim) : 0;
  LO n = CountDofHolders();
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
  SetGeometry(std::move(geometry));
}

template <unsigned Dim>
int UniformGridFieldLayout<Dim>::GetNumComponents() const
{
  return num_components_;
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::GetNumDofHolder() const
{
  return CountDofHolders();
}

template <unsigned Dim>
GO UniformGridFieldLayout<Dim>::GetNumGlobalDofHolder() const
{
  return CountDofHolders();
}

template <unsigned Dim>
Rank1View<const bool, HostMemorySpace>
UniformGridFieldLayout<Dim>::GetOwnedHost() const
{
  return pcms::detail::make_const_array_view(owned_host_);
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
MPI_Comm UniformGridFieldLayout<Dim>::GetCommunicator() const
{
  return MPI_COMM_SELF;
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
    num_vertices *= grid_.divisions[d] + (((wrapped_ >> d) & 1U) ? 0 : 1);
  }
  return num_vertices;
}

template <unsigned Dim>
LO UniformGridFieldLayout<Dim>::CountDofHolders() const
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
  return pcms::detail::make_const_array_view(classification_dims_host_);
}

template <unsigned Dim>
Rank1View<const LO, HostMemorySpace>
UniformGridFieldLayout<Dim>::GetDOFHolderClassificationIdsHost() const
{
  return pcms::detail::make_const_array_view(classification_ids_host_);
}

template <unsigned Dim>
EntOffsetsArray UniformGridFieldLayout<Dim>::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  LO n = CountDofHolders();
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
template class UniformGridFieldLayout<1>;
template class UniformGridFieldLayout<2>;
template class UniformGridFieldLayout<3>;

} // namespace pcms::uniform_grid
