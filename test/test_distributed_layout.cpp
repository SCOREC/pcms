// Checks the owned/ghost split of a FieldLayout on a ghosted, distributed
// Omega_h mesh: the DOF holders with ghosts are the owned ones plus their ghost
// copies, and the compact owned views
// (owned indices, owned GIDs, owned coordinates, owned field data) agree
// with the ones with ghosts.
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include "distributed_test_utils.h"
#include "pcms/configuration.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

using pcms::GO;
using pcms::LO;
using pcms::Real;

namespace
{

void require(bool condition, const std::string& message)
{
  if (!condition) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] %s\n", rank, message.c_str());
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

template <typename View>
auto ToHost(const View& view)
{
  return Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), view);
}

void CheckLayout(Omega_h::Mesh& mesh, const pcms::LagrangeFunctionSpace& space,
                 const std::string& label)
{
  const auto& layout = *space.GetLayout();
  // Order-1 layouts: the DOF holders with ghosts are the mesh vertices.
  const LO n_local = mesh.nverts();
  const LO n_owned = layout.GetNumDofHolder();
  const int nc = layout.GetNumComponents();
  const int dim = layout.GetDimension();

  require(layout.IsDistributed(), label + ": layout is not distributed");
  require(n_owned <= n_local, label + ": more owned than local holders");
  require(layout.Size() == n_owned * nc, label + ": Size mismatch");
  require(layout.HasGhosts() == (n_owned != n_local),
          label + ": HasGhosts mismatch");

  // Every holder is owned by exactly one rank, and a ghosted mesh has ghosts.
  const long long n_owned_ll = n_owned;
  const long long n_ghost_ll = n_local - n_owned;
  long long owned_sum = 0;
  long long ghost_sum = 0;
  MPI_Allreduce(&n_owned_ll, &owned_sum, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&n_ghost_ll, &ghost_sum, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  require(owned_sum == static_cast<long long>(layout.GetNumGlobalDofHolder()),
          label + ": owned counts do not sum to the global count");
  require(ghost_sum > 0, label + ": ghosted mesh produced no ghost holders");

  // owned indices: one entry per owned holder, increasing, owned only.
  const auto owned = layout.GetOwnedHost();
  const auto owned_indices = layout.GetOwnedIndicesHost();
  require(static_cast<LO>(owned_indices.size()) == n_owned,
          label + ": owned->local map has the wrong size");
  LO n_owned_flags = 0;
  for (LO i = 0; i < n_local; ++i) {
    n_owned_flags += owned[i] ? 1 : 0;
  }
  require(n_owned_flags == n_owned, label + ": owned flags != owned count");
  for (LO o = 0; o < n_owned; ++o) {
    const LO l = owned_indices(o);
    require(l >= 0 && l < n_local && owned[l],
            label + ": owned->local maps to a non-owned holder");
    require(o == 0 || owned_indices(o - 1) < l,
            label + ": owned->local map is not increasing");
  }
  const auto owned_to_local_device = ToHost(layout.GetOwnedIndices());
  require(static_cast<LO>(owned_to_local_device.size()) == n_owned,
          label + ": device owned->local map has the wrong size");
  for (LO o = 0; o < n_owned; ++o) {
    require(owned_to_local_device(o) == owned_indices(o),
            label + ": device and host owned->local maps differ");
  }

  // Owned GIDs, coordinates and classification are the mesh's per-vertex
  // arrays (ghosts included) at the mapped index. The layouts here are order-1,
  // so DOF holder i is mesh vertex i.
  const auto gids =
    Omega_h::HostRead<Omega_h::GO>(mesh.get_array<Omega_h::GO>(0, "global"));
  require(gids.size() == n_local, label + ": local count != mesh vertices");
  const auto owned_gids = layout.GetGidsHost();
  require(static_cast<LO>(owned_gids.size()) == n_owned,
          label + ": owned GIDs have the wrong size");
  for (LO o = 0; o < n_owned; ++o) {
    require(owned_gids[o] == gids[owned_indices(o)],
            label + ": owned GID != local GID at the mapped index");
  }

  const auto coords = Omega_h::HostRead<Real>(mesh.coords());
  const auto owned_coords = layout.GetDOFHolderCoordinates().GetValues();
  require(static_cast<LO>(owned_coords.extent(0)) == n_owned,
          label + ": owned coordinates have the wrong size");
  Kokkos::View<Real**, pcms::DeviceMemorySpace> owned_coords_copy(
    "owned_coords", n_owned, dim);
  pcms::CopyDeviceRank2ViewToRank2View(pcms::MakeRank2View(owned_coords_copy),
                                       owned_coords);
  const auto owned_coords_host = ToHost(owned_coords_copy);
  for (LO o = 0; o < n_owned; ++o) {
    for (int d = 0; d < dim; ++d) {
      require(owned_coords_host(o, d) == coords[owned_indices(o) * dim + d],
              label + ": owned coordinate != local coordinate");
    }
  }

  // Owned entity offsets cover exactly the owned holders.
  const auto owned_offsets = layout.GetEntOffsets();
  require(owned_offsets.back() == static_cast<size_t>(n_owned),
          label + ": owned entity offsets do not end at the owned count");

  // Field data: write the owned values and read them back; the backend keeps
  // them at their rows among the DOF holders with ghosts and leaves the ghosts
  // at their initial zero.
  auto field = space.CreateFunction<Real>();
  std::vector<Real> values(static_cast<size_t>(n_owned) * nc);
  for (LO o = 0; o < n_owned; ++o) {
    for (int c = 0; c < nc; ++c) {
      values[static_cast<size_t>(o) * nc + c] =
        static_cast<Real>(owned_gids[o]) + 0.25 * c;
    }
  }
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(values.data(), n_owned,
                                                       nc));

  const auto owned_host = field.GetDOFHolderDataHost().GetValues();
  require(static_cast<LO>(owned_host.extent(0)) == n_owned &&
            static_cast<int>(owned_host.extent(1)) == nc,
          label + ": owned host data has the wrong shape");
  const auto owned_device = field.GetDOFHolderData().GetValues();
  require(static_cast<LO>(owned_device.extent(0)) == n_owned &&
            static_cast<int>(owned_device.extent(1)) == nc,
          label + ": owned device data has the wrong shape");
  Kokkos::View<Real**, pcms::DeviceMemorySpace> owned_device_copy(
    "owned_device", n_owned, nc);
  pcms::CopyDeviceRank2ViewToRank2View(pcms::MakeRank2View(owned_device_copy),
                                       owned_device);
  const auto owned_device_host = ToHost(owned_device_copy);
  for (LO o = 0; o < n_owned; ++o) {
    for (int c = 0; c < nc; ++c) {
      const Real expected = static_cast<Real>(owned_gids[o]) + 0.25 * c;
      require(owned_host(o, c) == expected,
              label + ": owned host data != written value");
      require(owned_device_host(o, c) == expected,
              label + ": owned device data != written value");
    }
  }

  const auto local_values = distributed_test::ValuesWithGhostsHost(field);
  require(static_cast<LO>(local_values.extent(0)) == n_local,
          label + ": backend does not store every local holder");
  for (LO i = 0; i < n_local; ++i) {
    for (int c = 0; c < nc; ++c) {
      const Real expected =
        owned[i] ? static_cast<Real>(gids[i]) + 0.25 * c : Real(0);
      require(local_values(i, c) == expected,
              label + ": owned write missed an owned row or changed a ghost");
    }
  }

  // Owned classification is the mesh classification at the mapped index.
  const auto class_dims =
    Omega_h::HostRead<Omega_h::I8>(mesh.get_array<Omega_h::I8>(0, "class_dim"));
  const auto class_ids = Omega_h::HostRead<Omega_h::ClassId>(
    mesh.get_array<Omega_h::ClassId>(0, "class_id"));
  const auto owned_class_dims =
    layout.GetDOFHolderClassificationDimensionsHost();
  const auto owned_class_ids = layout.GetDOFHolderClassificationIdsHost();
  require(static_cast<LO>(owned_class_dims.size()) == n_owned &&
            static_cast<LO>(owned_class_ids.size()) == n_owned,
          label + ": owned classification has the wrong size");
  for (LO o = 0; o < n_owned; ++o) {
    require(owned_class_dims[o] == class_dims[owned_indices(o)] &&
              owned_class_ids[o] == class_ids[owned_indices(o)],
            label + ": owned classification != local classification");
  }
}

// Report `what` and abort every rank: a rank that just returned would leave
// the others blocked in their next collective.
[[noreturn]] void AbortAllRanks(const char* what) noexcept
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::fprintf(stderr, "[rank %d] uncaught exception: %s\n", rank, what);
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort(); // MPI_Abort does not return, but is not declared noreturn
}

} // namespace

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    {
      Omega_h::Library library(&argc, &argv);
      auto world = library.world();
      require(world->size() >= 2, "test requires at least 2 MPI ranks");

      const int n = 4 * world->size();
      auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, n,
                                     n, 0, false);
      mesh.set_parting(OMEGA_H_GHOSTED, 1, false);

      for (int nc : {1, 2}) {
        const std::string suffix = " (" + std::to_string(nc) + " components)";
        CheckLayout(mesh,
                    *pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 1, nc, pcms::csys::Cartesian::Deferred(), "global",
                      pcms::LagrangeFunctionSpace::Backend::OmegaH),
                    "Omega_h backend" + suffix);
#ifdef PCMS_ENABLE_MESHFIELDS
        CheckLayout(mesh,
                    *pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 1, nc, pcms::csys::Cartesian::Deferred(), "global",
                      pcms::LagrangeFunctionSpace::Backend::MeshFields),
                    "MeshFields backend" + suffix);
#endif
      }

      if (world->rank() == 0) {
        std::printf("distributed layout test passed\n");
      }
    }
    MPI_Finalize();
    return 0;
  } catch (const std::exception& e) {
    AbortAllRanks(e.what());
  } catch (...) {
    AbortAllRanks("unknown exception");
  }
}
