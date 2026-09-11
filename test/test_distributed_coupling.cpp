#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/coupler/field_serializer.h"
#include "pcms/field/data/simple.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

using pcms::LO;
using pcms::Real;

namespace
{
void require(bool cond, const char* msg)
{
  if (!cond) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] distributed coupling check failed: %s\n",
                 rank, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}
} // namespace

// Exercises the owned-only serialization contract used by distributed coupling:
// only owned (rank-exclusive) DOF holders are put on the wire, and received
// values are scattered back into the local (owned + ghost) field.
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

    // Distributed 2D simplex box with a ghost layer (owned < local).
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
    const LO n_comp = layout.GetNumComponents();
    const auto owned_to_local = layout.GetOwnedToLocalHost();
    require(n_owned > 0, "rank owns at least one DOF holder");
    require(n_local > n_owned, "ghosted mesh has local > owned");

    // Set local field data: value = local index.
    auto field = space->CreateFunction<Real>();
    std::vector<Real> local_data(static_cast<size_t>(n_local));
    for (LO i = 0; i < n_local; ++i) {
      local_data[static_cast<size_t>(i)] = static_cast<Real>(i);
    }
    field.SetDOFHolderDataHost(
      pcms::Rank2View<const Real, pcms::HostMemorySpace>(local_data.data(),
                                                         n_local, n_comp));

    // Owned-indexed permutation: the last owned holder is "outside the overlap
    // region" (perm = -1); the rest get compact buffer slots
    // 0..n_participating.
    const LO n_participating = n_owned - 1;
    std::vector<LO> permutation(static_cast<size_t>(n_owned), -1);
    for (LO o = 0; o < n_participating; ++o) {
      permutation[static_cast<size_t>(o)] = o;
    }

    // Serialize: only participating owned data is written to the buffer.
    pcms::FieldSerializer<Real> serializer;
    std::vector<Real> buffer(static_cast<size_t>(n_participating) * n_comp);
    const int sent = serializer.Serialize(
      field.GetData(), layout, pcms::make_array_view(buffer),
      pcms::make_const_array_view(permutation));
    require(sent == static_cast<int>(n_participating) * n_comp,
            "serialized size matches participating holders");

    auto owned_field =
      pcms::FlattenToRank1View(field.GetOwnedDOFHolderDataHost());
    require(owned_field.size() == static_cast<size_t>(n_owned),
            "owned field size");
    for (LO o = 0; o < n_participating; ++o) {
      require(buffer[static_cast<size_t>(o)] ==
                owned_field[static_cast<size_t>(o)],
              "serialized value == owned value");
    }

    // Deserialize into a fresh field: participating owned values land at their
    // local positions; ghost and non-participating owned slots stay zero.
    auto target = space->CreateFunction<Real>();
    serializer.Deserialize(target.GetData(), layout,
                           pcms::make_const_array_view(buffer),
                           pcms::make_const_array_view(permutation));
    auto target_data = pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    require(target_data.size() == static_cast<size_t>(n_local),
            "target local size");

    for (LO o = 0; o < n_participating; ++o) {
      const LO local = owned_to_local(o);
      require(target_data[static_cast<size_t>(local)] ==
                local_data[static_cast<size_t>(local)],
              "participating owned value restored");
    }

    const LO non_participating_local = owned_to_local(n_owned - 1);
    require(target_data[static_cast<size_t>(non_participating_local)] ==
              Real(0),
            "non-participating owned holder is zero");

    auto owned_mask = layout.GetOwnedHost();
    for (LO i = 0; i < n_local; ++i) {
      if (!owned_mask[static_cast<size_t>(i)]) {
        require(target_data[static_cast<size_t>(i)] == Real(0),
                "ghost holder is zero");
      }
    }

    if (rank == 0) {
      std::printf("distributed coupling test passed (nproc=%d)\n", nproc);
    }
  }
  MPI_Finalize();
  return result;
}
