/**
 * @file test_distributed_field_eval_client.cpp
 * @brief Client side of the distributed field evaluation cross-rank test.
 *
 * It must be paired with a running distributed_field_eval_server (see
 * pcms/tools/distributed_field_eval_server.cpp) using the same channel name.
 *
 * ctest wires this up automatically as `test_distributed_field_eval` (see
 * pcms/test/CMakeLists.txt): a single-rank server job plus this job, launched
 * concurrently. Both sides declare the same round count, and
 * DistributedEvaluationChannel aborts loudly on a mismatch.
 *
 * No partition input is needed: DistributedPointEvaluator publishes each
 * rank's owned bounding box over the channel, and the server assembles its
 * routing partition from those boxes (see
 * distributed_evaluation_channel.h's class comment). The box this program
 * prints to stderr is purely diagnostic.
 */
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <redev.h>
#include "pcms/coupler/distributed_evaluation_channel.h"
#include "pcms/coupler/distributed_evaluation_partition.h"
#include "pcms/coupler/distributed_point_evaluator.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

using pcms::LO;
using pcms::Real;

namespace
{
KOKKOS_INLINE_FUNCTION Real linear_f(Real x, Real y)
{
  return x + 2.0 * y;
}

void require(bool condition, const char* message)
{
  if (!condition) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] %s\n", rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

// A 5x5 grid spanning the whole [0,1]^2 domain, identical on every rank: some
// of these points are local to a given rank (found by its own mesh/ghost
// layer), others belong to a different rank and can only be answered through
// the distributed evaluation server.
pcms::CoordinateView<pcms::DeviceMemorySpace> MakeGlobalSampleCoords(
  Kokkos::View<Real**, pcms::DeviceMemorySpace>& storage)
{
  constexpr int kSamplesPerAxis = 5;
  Kokkos::View<Real**, pcms::HostMemorySpace> host(
    "sample_coords", kSamplesPerAxis * kSamplesPerAxis, 2);
  int k = 0;
  for (int i = 0; i < kSamplesPerAxis; ++i) {
    for (int j = 0; j < kSamplesPerAxis; ++j, ++k) {
      host(k, 0) = (i + 0.5) / kSamplesPerAxis;
      host(k, 1) = (j + 0.5) / kSamplesPerAxis;
    }
  }
  storage = Kokkos::View<Real**, pcms::DeviceMemorySpace>(
    "sample_coords_device", host.extent(0), host.extent(1));
  pcms::DeepCopyMismatchLayouts(storage, host);
  return {pcms::CoordinateSystem::Cartesian, pcms::MakeRank2View(storage)};
}
} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  require(argc >= 2, "usage: test_distributed_field_eval_client <channel_name>");
  const std::string channel_name = argv[1];
  {
    Omega_h::Library library(&argc, &argv);
    auto world = library.world();
    int rank = world->rank();

    const int n = 4 * world->size();
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, n, n, 0, false);
    mesh.set_parting(OMEGA_H_GHOSTED, 1, false);
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::CoordinateSystem::Cartesian, "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    const auto& layout = *space->GetLayout();

    // Diagnostic only: report this rank's own box. The server does not need
    // this printout -- the evaluator publishes the same box over the channel
    // (DistributedBoundingBoxPartition::LocalOwnedBounds) -- but it is handy
    // when reading a test log.
    std::array<Real, 3> box_min{};
    std::array<Real, 3> box_max{};
    bool owns_nothing = false;
    pcms::DistributedBoundingBoxPartition::LocalOwnedBounds(
      layout, mesh.dim(), box_min, box_max, owns_nothing);
    std::fprintf(stderr, "[rank %d] owns [%.6f, %.6f] x [%.6f, %.6f]%s\n",
                 rank, box_min[0], box_max[0], box_min[1], box_max[1],
                 owns_nothing ? " (owns nothing)" : "");

    // Initialize the field with an affine function on owned holders, then
    // synchronize ghosts, exactly as test_distributed_field.cpp does.
    const LO n_local = layout.GetNumLocalDofHolder();
    auto owned = layout.GetOwnedHost();
    std::vector<Real> values_init(static_cast<size_t>(n_local));
    auto coords_mirror =
      pcms::detail::CopyToHostView(layout.GetDOFHolderCoordinates().GetValues());
    for (LO i = 0; i < n_local; ++i) {
      values_init[static_cast<size_t>(i)] =
        owned(i) ? linear_f(coords_mirror(i, 0), coords_mirror(i, 1)) : 0.0;
    }
    auto field = space->CreateFunction<Real>();
    field.SetDOFHolderDataHost(pcms::Rank2View<const Real, pcms::HostMemorySpace>(
      values_init.data(), n_local, 1));
    field.SynchronizeGhosts();

    // Set up the redev Client role and the matching distributed-evaluation
    // channel; "channel_name" must match what the paired server was launched
    // with.
    redev::Redev redev(MPI_COMM_WORLD);
    adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
    auto channel = redev.CreateAdiosChannel(channel_name, params,
                                            redev::TransportType::BP4);
    pcms::DistributedEvaluationChannel<Real> eval_channel(
      redev, channel, MPI_COMM_WORLD, channel_name, mesh.dim(),
      /*num_components=*/1,
      /*expected_rounds=*/1);

    Kokkos::View<Real**, pcms::DeviceMemorySpace> sample_storage;
    auto sample_coords = MakeGlobalSampleCoords(sample_storage);

    pcms::DistributedPointEvaluator<Real> evaluator(
      space->GetEvaluatorFactory(), sample_coords, eval_channel,
      pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::ERROR});

    const LO n_samples = static_cast<LO>(sample_storage.extent(0));
    Kokkos::View<Real**, pcms::DeviceMemorySpace> result("result", n_samples,
                                                         1);
    evaluator.Evaluate(field, pcms::MakeRank2View(result));

    auto result_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), result);
    auto sample_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), sample_storage);
    for (LO i = 0; i < n_samples; ++i) {
      const Real expected = linear_f(sample_host(i, 0), sample_host(i, 1));
      require(std::fabs(result_host(i, 0) - expected) < 1.0e-9,
              "distributed evaluation (with cross-rank fallback) returned "
              "the wrong value");
    }
    if (rank == 0) {
      std::printf("distributed field evaluation (cross-rank) passed\n");
    }
  }
  MPI_Finalize();
  return 0;
}
