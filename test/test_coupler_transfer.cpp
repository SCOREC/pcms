#include <Omega_h_mesh.hpp>
#include <cmath>
#include <iostream>
#include <pcms.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include "pcms/coupler/coupler.hpp"
#include "pcms/transfer/transfer_method.hpp"
#include "pcms/field/function_space/lagrange.h"

using pcms::Real;
namespace ts = test_support;

static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 4;

namespace
{

std::shared_ptr<pcms::LagrangeFunctionSpace> MakeSpace(Omega_h::Mesh& mesh)
{
  // Both ends name the layout "field" so their GID exchanges connect.
  return pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::DefaultBackend, "field");
}

void SetFieldToGids(const pcms::FieldLayout& layout,
                    pcms::FieldData<Real>* field)
{
  auto gids = layout.GetGidsHost();
  const auto n = layout.GetNumOwnedDofHolder();
  Omega_h::HostWrite<Real> values(n);
  Kokkos::parallel_for(
    "set_gids",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { values[i] = gids[i]; });
  field->SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(values.data(), n, 1));
}

bool FieldEqualsGids(const pcms::FieldLayout& layout,
                     pcms::FieldData<Real>* field, int rank)
{
  auto gids = layout.GetGidsHost();
  auto owned = layout.GetOwnedHost();
  auto got = pcms::FlattenToRank1View(field->GetDOFHolderDataHost());
  const auto n = layout.GetNumOwnedDofHolder();

  int expected = 0, matched = 0;
  Kokkos::parallel_reduce(
    "count_owned",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0,
                                                                owned.size()),
    KOKKOS_LAMBDA(int i, int& s) { s += owned[i] != 0; }, expected);
  Kokkos::parallel_reduce(
    "count_match",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(int i, int& s) {
      // Identity interpolation of an arbitrary per-vertex GID field: in exact
      // arithmetic each target vertex recovers its source value, but the
      // point-in-element search leaves a tiny (~1e-7) barycentric perturbation
      // that, scaled by the large unordered GID differences, can nudge a value
      // just below its integer. Round rather than truncate to compare.
      if (owned[i])
        s += std::llround(got[i]) == (long long)gids[i];
    },
    matched);

  const bool ok = matched == expected;
  std::cerr << "Rank " << rank << " - target field validation "
            << (ok ? "PASSED" : "FAILED") << " (" << matched << "/" << expected
            << ")\n";
  return ok;
}

} // namespace

// --------------------------------------------------------------------------
// The coupling server: receives a field from the source app, transfers it onto
// the target app's space, and sends it on. This is the whole coupler-integrated
// transfer story.
// --------------------------------------------------------------------------
void transfer_server(MPI_Comm comm, Omega_h::Mesh& mesh,
                     std::string_view cpn_file)
{
  pcms::Coupler cpl("coupler_transfer", comm, /*isServer=*/true,
                    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  auto* source_app = cpl.AddApplication("coupler_transfer_source");
  auto* target_app = cpl.AddApplication("coupler_transfer_target");

  auto source_space = MakeSpace(mesh);
  auto target_space = MakeSpace(mesh);

  auto source =
    source_app->AddFunction(source_space->CreateFunction<Real>("field"));
  auto target =
    target_app->AddFunction(target_space->CreateFunction<Real>("field"));

  auto transfer =
    cpl.AddTransfer(source, target, pcms::method::Interpolation<Real>{});

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      source_app->ReceivePhase([&] { source.Receive(); });
      transfer.Run(); // source field -> target field
      target_app->SendPhase([&] { target.Send(); });
    }
  } while (!done);
}

void transfer_source_client(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  pcms::Coupler coupler("coupler_transfer", comm, false, {});
  auto* app = coupler.AddApplication("coupler_transfer_source");
  auto factory = MakeSpace(mesh);

  auto field = factory->CreateFunction<Real>("field");
  auto* field_ptr = &field.GetData();
  SetFieldToGids(*factory->GetLayout(), field_ptr);
  app->AddField(std::move(field));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->SendPhase([&] { app->SendField("field"); });
    }
  } while (!done);
}

void transfer_target_client(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  pcms::Coupler coupler("coupler_transfer", comm, false, {});
  auto* app = coupler.AddApplication("coupler_transfer_target");
  auto factory = MakeSpace(mesh);

  auto field = factory->CreateFunction<Real>("field");
  auto* field_ptr = &field.GetData();
  app->AddField(std::move(field));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->ReceivePhase([&] { app->ReceiveField("field"); });
      // Identity interpolation on the same mesh: the target must receive the
      // exact field the source sent.
      if (!FieldEqualsGids(*factory->GetLayout(), field_ptr, rank)) {
        exit(EXIT_FAILURE);
      }
    }
  } while (!done);
}

int main(int argc, char** argv)
{
  try {
    auto lib = Omega_h::Library(&argc, &argv);
    const int rank = lib.world()->rank();
    if (argc != 4) {
      if (!rank) {
        std::cerr
          << "Usage: " << argv[0]
          << " <clientId=-1|0|1> /path/to/mesh /path/to/partition.cpn\n";
      }
      exit(EXIT_FAILURE);
    }
    const auto clientId = atoi(argv[1]);
    REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
    Omega_h::Mesh mesh(&lib);
    Omega_h::binary::read(argv[2], lib.world(), &mesh);
    MPI_Comm mpi_comm = lib.world()->get_impl();
    switch (clientId) {
      case -1: transfer_server(mpi_comm, mesh, argv[3]); break;
      case 0: transfer_source_client(mpi_comm, mesh); break;
      case 1: transfer_target_client(mpi_comm, mesh); break;
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  }
}
