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

// Fill a field with scale*gid at each dof holder. The scale lets two fields on
// the same space carry distinct data through one shared transfer operator.
void SetFieldToGids(const pcms::FieldLayout& layout,
                    pcms::FieldData<Real>* field, Real scale = 1.0)
{
  auto gids = layout.GetGidsHost();
  const auto n = layout.GetNumOwnedDofHolder();
  Omega_h::HostWrite<Real> values(n);
  Kokkos::parallel_for(
    "set_gids",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { values[i] = scale * gids[i]; });
  field->SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(values.data(), n, 1));
}

bool FieldEqualsGids(const pcms::FieldLayout& layout,
                     pcms::FieldData<Real>* field, int rank, Real scale,
                     const char* label)
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
        s += std::llround(got[i]) == std::llround(scale * gids[i]);
    },
    matched);

  const bool ok = matched == expected;
  std::cerr << "Rank " << rank << " - target field '" << label
            << "' validation " << (ok ? "PASSED" : "FAILED") << " (" << matched
            << "/" << expected << ")\n";
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

  // Two fields (think displacement and velocity) that live on the same source
  // and target spaces. Each pair is registered from the same space objects, so
  // both share a function space -- and therefore one transfer operator can
  // serve both.
  auto source =
    source_app->AddFunction(source_space->CreateFunction<Real>("field"));
  auto target =
    target_app->AddFunction(target_space->CreateFunction<Real>("field"));
  auto source2 =
    source_app->AddFunction(source_space->CreateFunction<Real>("field2"));
  auto target2 =
    target_app->AddFunction(target_space->CreateFunction<Real>("field2"));

  // Build the interpolation operator ONCE from the shared (source, target)
  // space pair -- this is the expensive localization step -- then bind that one
  // operator to each field pair. transfer and transfer2 reuse the same cached
  // localization; there is no second Build().
  std::shared_ptr<const pcms::TransferOperator<Real>> interp =
    pcms::method::Interpolation<Real>{}.Build(source.GetSpace(),
                                              target.GetSpace());
  // Named transfers: the returned handle and the name are interchangeable ways
  // to drive them (see the two Run styles below).
  auto transfer = cpl.AddTransfer("field", interp, source, target);
  cpl.AddTransfer("field2", interp, source2, target2);

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      source_app->ReceivePhase([&] {
        source.Receive();
        source2.Receive();
      });
      transfer.Run();            // via handle: field  source -> target
      cpl.RunTransfer("field2"); // via name:   field2 source -> target
      target_app->SendPhase([&] {
        target.Send();
        target2.Send();
      });
    }
  } while (!done);
}

void transfer_source_client(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  pcms::Coupler coupler("coupler_transfer", comm, false, {});
  auto* app = coupler.AddApplication("coupler_transfer_source");
  auto factory = MakeSpace(mesh);

  auto field = factory->CreateFunction<Real>("field");
  SetFieldToGids(*factory->GetLayout(), &field.GetData());
  app->AddField(std::move(field));

  // A second field on the same space, with a distinct scale so it cannot be
  // confused with the first as it flows through the shared operator.
  auto field2 = factory->CreateFunction<Real>("field2");
  SetFieldToGids(*factory->GetLayout(), &field2.GetData(), 3.0);
  app->AddField(std::move(field2));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->SendPhase([&] {
        app->SendField("field");
        app->SendField("field2");
      });
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

  auto field2 = factory->CreateFunction<Real>("field2");
  auto* field2_ptr = &field2.GetData();
  app->AddField(std::move(field2));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->ReceivePhase([&] {
        app->ReceiveField("field");
        app->ReceiveField("field2");
      });
      // Identity interpolation on the same mesh: each target must receive the
      // exact field the source sent. The two fields carry different scales, so
      // if the shared operator crossed their data the mismatch would show here.
      if (!FieldEqualsGids(*factory->GetLayout(), field_ptr, rank, 1.0,
                           "field") ||
          !FieldEqualsGids(*factory->GetLayout(), field2_ptr, rank, 3.0,
                           "field2")) {
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
      default: break; // unreachable: clientId asserted to [-1, 1] above
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  }
}
