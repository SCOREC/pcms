#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/field/data/simple.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>

using pcms::GO;
using pcms::LO;
using pcms::Real;

namespace
{
void require(bool cond, const char* msg)
{
  if (!cond) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] distributed field check failed: %s\n", rank,
                 msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

// Copy a device Rank2View (mdspan-like) to a host Kokkos view, staging through
// a same-layout device scratch so host/device layout mismatches are handled.
template <typename CoordView>
Kokkos::View<Real**, pcms::HostMemorySpace> copy_coords_to_host(
  const CoordView& coords, int n, int dim)
{
  Kokkos::View<Real**, pcms::HostMemorySpace> host("coords_host", n, dim);
  auto device = Kokkos::create_mirror_view(pcms::DeviceMemorySpace(), host);
  Kokkos::parallel_for(
    n, KOKKOS_LAMBDA(int i) {
      for (int d = 0; d < dim; ++d)
        device(i, d) = coords(i, d);
    });
  Kokkos::deep_copy(host, device);
  return host;
}
} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int result = 0;
  {
    Kokkos::ScopeGuard kokkos{argc, argv};
    Omega_h::Library lib(&argc, &argv);
    auto world = lib.world();
    const int rank = world->rank();
    const int nproc = world->size();
    require(nproc >= 2, "test requires at least 2 MPI ranks");

    // Distributed 2D simplex box with a ghost layer so each rank has owned and
    // ghost (non-owned) vertices.
    Omega_h::Mesh mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0,
                                            0.0, /*nx=*/8, /*ny=*/8, /*nz=*/0,
                                            /*symmetric=*/false);
    mesh.set_parting(OMEGA_H_GHOSTED, 1, false);

    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, /*order=*/1, /*num_components=*/1,
      pcms::CoordinateSystem::Cartesian, "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    const pcms::FieldLayout& layout = *space->GetLayout();

    const LO n_owned = layout.GetNumOwnedDofHolder();
    const LO n_local = layout.GetNumLocalDofHolder();

    require(layout.IsDistributed(), "OmegaH layout should be distributed");
    // Count owned vertices from the owned mask on host. Avoid nents_owned(),
    // which reads the device-owned array from host in some Omega_h/CUDA builds.
    auto owned_mask_omega = Omega_h::HostRead<Omega_h::I8>(mesh.owned(0));
    LO n_owned_omega = 0;
    for (LO i = 0; i < mesh.nents(0); ++i) {
      if (owned_mask_omega[i] != 0)
        ++n_owned_omega;
    }
    require(n_owned == n_owned_omega, "owned count != owned mask count");
    require(n_local == mesh.nents(0), "local count != nents(0)");
    require(layout.GetNumGlobalDofHolder() == mesh.nglobal_ents(0),
            "global count != nglobal_ents(0)");

    // Total owned vertices < total local (owned + ghost) vertices.
    GO g_owned = 0, g_local = 0;
    const GO owned64 = static_cast<GO>(n_owned);
    const GO local64 = static_cast<GO>(n_local);
    MPI_Allreduce(&owned64, &g_owned, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local64, &g_local, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    std::cout << "g_owned: " << g_owned << ", g_local: " << g_local
              << std::endl;
    require(g_owned < g_local, "global owned !< global local (no ghosts?)");

    // Owned accessors are sized by the owned count.
    auto owned_gids = layout.GetOwnedGidsHost();
    auto local_gids = layout.GetGidsHost();
    auto owned_to_local = layout.GetOwnedToLocalHost();
    std::cout << "owned_gids.size(): " << owned_gids.size()
              << ", local_gids.size(): " << local_gids.size()
              << ", owned_to_local.size(): " << owned_to_local.size()
              << std::endl;
    require(owned_gids.size() == static_cast<size_t>(n_owned),
            "owned gids size");
    require(local_gids.size() == static_cast<size_t>(n_local),
            "local gids size");
    require(owned_to_local.size() == static_cast<size_t>(n_owned),
            "owned_to_local size");

    auto owned_coords = layout.GetOwnedDOFHolderCoordinates().GetValues();
    auto local_coords = layout.GetDOFHolderCoordinates().GetValues();
    std::cout << "owned_coords.extent(0): " << owned_coords.extent(0)
              << ", local_coords.extent(0): " << local_coords.extent(0)
              << std::endl;
    require(owned_coords.extent(0) == static_cast<size_t>(n_owned),
            "owned coords rows");
    require(local_coords.extent(0) == static_cast<size_t>(n_local),
            "local coords rows");

    // owned_to_local maps each owned holder to an owned local holder, and the
    // owned GID equals the local GID at that index.
    auto owned_mask = layout.GetOwnedHost();
    require(owned_mask.size() == static_cast<size_t>(n_local),
            "owned mask size");
    for (LO o = 0; o < n_owned; ++o) {
      const LO local = owned_to_local(o);
      require(local >= 0 && local < n_local, "owned_to_local out of range");
      require(owned_mask[static_cast<size_t>(local)],
              "owned_to_local maps to a non-owned holder");
      require(owned_gids(o) == local_gids(local), "owned gid != local gid");
    }

    // Owned coordinates equal the local coordinates at the mapped index.
    auto owned_coords_h =
      copy_coords_to_host(owned_coords, static_cast<int>(n_owned), mesh.dim());
    auto local_coords_h =
      copy_coords_to_host(local_coords, static_cast<int>(n_local), mesh.dim());
    for (LO o = 0; o < n_owned; ++o) {
      const LO local = owned_to_local(o);
      for (int d = 0; d < mesh.dim(); ++d) {
        require(owned_coords_h(o, d) == local_coords_h(local, d),
                "owned coord != local coord");
      }
    }

    // Field data: the owned data view is the owned subset of the local view.
    auto field = space->CreateFunction<Real>();
    std::vector<Real> local_data(static_cast<size_t>(n_local));
    for (LO i = 0; i < n_local; ++i) {
      local_data[static_cast<size_t>(i)] = static_cast<Real>(i);
    }
    field.SetDOFHolderDataHost(
      pcms::Rank2View<const Real, pcms::HostMemorySpace>(local_data.data(),
                                                         n_local, 1));

    auto local_field = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
    auto owned_field =
      pcms::FlattenToRank1View(field.GetOwnedDOFHolderDataHost());
    std::cout << "local_field.size(): " << local_field.size()
              << ", owned_field.size(): " << owned_field.size() << std::endl;
    require(local_field.size() == static_cast<size_t>(n_local),
            "local field size");
    require(owned_field.size() == static_cast<size_t>(n_owned),
            "owned field size");
    for (LO o = 0; o < n_owned; ++o) {
      const LO local = owned_to_local(o);
      require(owned_field[static_cast<size_t>(o)] ==
                local_field[static_cast<size_t>(local)],
              "owned field value != local field value");
    }

    // Ghost synchronization: set each rank's owned values to its rank id,
    // synchronize, and verify every ghost value matches its owner's rank id.
    auto sync_field = space->CreateFunction<Real>();
    std::vector<Real> sync_vals(static_cast<size_t>(n_local), Real(0));
    for (LO i = 0; i < n_local; ++i) {
      if (owned_mask[static_cast<size_t>(i)]) {
        sync_vals[static_cast<size_t>(i)] = static_cast<Real>(rank);
      }
    }
    sync_field.SetDOFHolderDataHost(
      pcms::Rank2View<const Real, pcms::HostMemorySpace>(sync_vals.data(),
                                                         n_local, 1));
    sync_field.SynchronizeGhosts();

    auto remotes = mesh.ask_owners(0);
    auto owners = Omega_h::HostRead<Omega_h::I32>(remotes.ranks);
    auto sync_data =
      pcms::FlattenToRank1View(sync_field.GetDOFHolderDataHost());
    for (LO i = 0; i < n_local; ++i) {
      if (!owned_mask[static_cast<size_t>(i)]) {
        require(sync_data[static_cast<size_t>(i)] ==
                  static_cast<Real>(owners[static_cast<size_t>(i)]),
                "ghost value != owner rank after sync");
      }
    }

    if (rank == 0) {
      std::printf("distributed field test passed (nproc=%d)\n", nproc);
    }
  }
  MPI_Finalize();
  return result;
}
