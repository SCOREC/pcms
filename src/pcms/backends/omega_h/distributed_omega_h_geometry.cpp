#include "pcms/backends/omega_h/distributed_omega_h_geometry.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::omega_h
{

namespace
{

// The points copied into the fixed-dimension layout the search takes.
template <int Dim>
Kokkos::View<Real* [Dim], DeviceMemorySpace> ToSearchCoords(
  Kokkos::View<const Real**, DeviceMemorySpace> points) {
  Kokkos::View<Real* [Dim], DeviceMemorySpace> out("search coords",
                                                   points.extent(0));
  Kokkos::parallel_for(
    "pcms::ToSearchCoords",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0,
                                                            points.extent(0)),
    KOKKOS_LAMBDA(LO i) {
      for (int d = 0; d < Dim; ++d) {
        out(i, d) = points(i, d);
      }
    });
  return out;
}

// The communicator of `local`'s mesh.
MPI_Comm CommOf(const std::shared_ptr<const OmegaHGeometry>& local)
{
  if (local == nullptr) {
    throw pcms_error(
      "DistributedOmegaHGeometry: the local geometry must not be null");
  }
  return local->GetMesh().comm()->get_impl();
}

} // namespace

DistributedOmegaHGeometry::DistributedOmegaHGeometry(
  std::shared_ptr<const OmegaHGeometry> local, LO grid_divisions)
  : DistributedGeometry(local, CommOf(local)), mesh_(local->GetMesh())
{
  const auto& local_geometry =
    static_cast<const OmegaHGeometry&>(*GetLocalGeometry());
  const Real tolerance = local_geometry.GetSearchOptions().tolerance;
  const auto& local_search = local_geometry.Search();
  if (mesh_.dim() == 2) {
    search_ = std::make_unique<DistributedPointSearch<2>>(
      mesh_, std::get<GridPointSearch2D>(local_search), tolerance,
      grid_divisions);
  } else {
    search_ = std::make_unique<DistributedPointSearch<3>>(
      mesh_, std::get<GridPointSearch3D>(local_search), tolerance,
      grid_divisions);
  }
}

std::vector<PointOwner> DistributedOmegaHGeometry::LocateOwners(
  Kokkos::View<const Real**, DeviceMemorySpace> points,
  const OwnerHints& owners) const
{
  const auto guesses =
    owners.Names(*GetLocalGeometry(), static_cast<LO>(points.extent(0)))
      ? detail::ElementGuesses(mesh_, owners, "DistributedOmegaHGeometry")
      : Kokkos::View<const LO*, DeviceMemorySpace>{};
  return std::visit(
    [&](const auto& search) {
      using Search = std::decay_t<decltype(*search)>;
      constexpr int Dim = Search::LocalSearch::Dimension;
      const auto located = search->Locate(ToSearchCoords<Dim>(points), guesses);
      std::vector<PointOwner> out(located.size());
      for (size_t i = 0; i < located.size(); ++i) {
        if (located[i].status == Search::ResultStatus::LOCATED) {
          out[i] = {located[i].owner_rank, located[i].element_id};
        }
      }
      return out;
    },
    search_);
}

std::shared_ptr<const DistributedGeometry> CreateDistributedGeometry(
  std::shared_ptr<const Geometry> local, LO grid_divisions)
{
  auto omega_h_local =
    std::dynamic_pointer_cast<const OmegaHGeometry>(std::move(local));
  if (omega_h_local == nullptr) {
    throw pcms_error("CreateDistributedGeometry: the local geometry is not "
                     "an Omega_h geometry");
  }
  return std::make_shared<DistributedOmegaHGeometry>(std::move(omega_h_local),
                                                     grid_divisions);
}

std::shared_ptr<const DistributedGeometry> CreateDistributedGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options, LO grid_divisions)
{
  return CreateDistributedGeometry(
    CreateGeometry(mesh, std::move(system), std::move(options)),
    grid_divisions);
}

} // namespace pcms::omega_h
