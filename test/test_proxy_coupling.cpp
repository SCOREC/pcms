#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include "pcms/coupler/coupler.hpp"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include <chrono>
#include <thread>

using pcms::GO;
using pcms::make_array_view;

using namespace std::chrono_literals;

static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 4;
namespace ts = test_support;

void initializeFieldWithGids(const pcms::FieldLayout& layout,
                             pcms::FieldData<pcms::Real>* field,
                             pcms::Real multiplier = 1.0)
{
  auto gids = layout.GetGidsHost();
  const auto n = layout.GetNumOwnedDofHolder();

  Omega_h::HostWrite<pcms::Real> ids(n);
  PCMS_ALWAYS_ASSERT(n == gids.size());
  Kokkos::parallel_for(
    "init gid",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { ids[i] = gids[i] * multiplier; });

  field->SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(ids.data(), n, 1));
}

bool validateField(const pcms::FieldLayout& layout,
                   pcms::FieldData<pcms::Real>* field,
                   const std::string& field_name, int rank,
                   pcms::Real multiplier = 1.0)
{
  auto gids = layout.GetGidsHost();
  auto copied_array = pcms::FlattenToRank1View(field->GetDOFHolderDataHost());
  auto owned = layout.GetOwnedHost();
  const auto n = layout.GetNumOwnedDofHolder();

  PCMS_ALWAYS_ASSERT(copied_array.size() == gids.size());
  PCMS_ALWAYS_ASSERT(owned.size() == gids.size());

  int expected = 0;
  int sum = 0;
  Kokkos::parallel_reduce(
    "count_owned",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0,
                                                                owned.size()),
    KOKKOS_LAMBDA(int i, int& local_sum) { local_sum += owned[i] != 0; },
    expected);
  Kokkos::parallel_reduce(
    "count_match",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(int i, int& local_sum) {
      if (owned[i])
        local_sum += (int)copied_array[i] == (int)(gids[i] * multiplier);
    },
    sum);

  if (sum != expected) {
    std::cerr << "Rank " << rank << " - Field '" << field_name
              << "' validation FAILED: expected " << expected
              << " matches, got " << sum << std::endl;
    return false;
  } else {
    std::cerr << "Rank " << rank << " - Field '" << field_name
              << "' validation PASSED (" << sum << " values)\n";
    return true;
  }
}

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_delta_f");

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  app->AddLayout("gids", factory.GetLayout());

  auto gids_field = factory.CreateField<pcms::Real>(pcms::FieldMetadata{});
  auto gids2_field = factory.CreateField<pcms::Real>(pcms::FieldMetadata{});

  auto* gids_ptr = &gids_field.GetData();
  auto* gids2_ptr = &gids2_field.GetData();

  initializeFieldWithGids(*factory.GetLayout(), gids_ptr, 1.0);
  initializeFieldWithGids(*factory.GetLayout(), gids2_ptr, 2.0);

  app->AddField("gids", std::move(gids_field));
  app->AddField("gids2", std::move(gids2_field));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids"); //(Alt) df_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) df_gid_field->Receive();
      app->EndReceivePhase();

      if (!validateField(*factory.GetLayout(), gids_ptr, "gids", rank, 1.0)) {
        std::cerr << "xgc_delta_f: Field validation failed at round " << i
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } while (!done);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_total_f");

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  app->AddLayout("gids", factory.GetLayout());

  auto gids_field = factory.CreateField<pcms::Real>(pcms::FieldMetadata{});

  auto* gids_ptr = &gids_field.GetData();

  initializeFieldWithGids(*factory.GetLayout(), gids_ptr, 10.0);

  app->AddField("gids", std::move(gids_field));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids"); //(Alt) tf_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) tf_gid_field->Receive();
      app->EndReceivePhase();

      if (!validateField(*factory.GetLayout(), gids_ptr, "gids", rank, 10.0)) {
        std::cerr << "xgc_total_f: Field validation failed at round " << i
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } while (!done);
}
void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl("proxy_couple", comm, true,
                    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");
  auto factory_total = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  total_f->AddLayout("gids", factory_total.GetLayout());

  auto factory_delta = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  delta_f->AddLayout("gids", factory_delta.GetLayout());
  // TODO, fields should have a transfer policy rather than parameters
  auto total_gids_field =
    factory_total.CreateField<pcms::Real>(pcms::FieldMetadata{});
  auto delta_gids_field =
    factory_delta.CreateField<pcms::Real>(pcms::FieldMetadata{});
  auto delta_gids2_field =
    factory_delta.CreateField<pcms::Real>(pcms::FieldMetadata{});

  auto* total_gids_ptr = &total_gids_field.GetData();
  auto* delta_gids_ptr = &delta_gids_field.GetData();
  auto* delta_gids2_ptr = &delta_gids2_field.GetData();

  total_f->AddField("gids", std::move(total_gids_field));
  delta_f->AddField("gids", std::move(delta_gids_field));
  delta_f->AddField("gids2", std::move(delta_gids2_field));

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      total_f->BeginReceivePhase();
      total_f->ReceiveField("gids");
      total_f->EndReceivePhase();

      if (!validateField(*factory_total.GetLayout(), total_gids_ptr, "gids",
                         rank, 10.0)) {
        std::cerr << "xgc_coupler: total_f field validation failed at round "
                  << i << std::endl;
        exit(EXIT_FAILURE);
      }

      delta_f->BeginReceivePhase();
      delta_f->ReceiveField("gids");
      delta_f->EndReceivePhase();

      if (!validateField(*factory_delta.GetLayout(), delta_gids_ptr, "gids",
                         rank, 1.0)) {
        std::cerr << "xgc_coupler: delta_f field validation failed at round "
                  << i << std::endl;
        exit(EXIT_FAILURE);
      }

      total_f->BeginSendPhase();
      total_f->SendField("gids");
      total_f->EndSendPhase();

      delta_f->BeginSendPhase();
      delta_f->SendField("gids", pcms::Mode::Deferred);
      delta_f->EndSendPhase();
    }
  } while (!done);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}

int main(int argc, char** argv)
{
  try {
    auto lib = Omega_h::Library(&argc, &argv);
    auto world = lib.world();
    const int rank = world->rank();
    if (argc != 4) {
      if (!rank) {
        std::cerr << "Usage: " << argv[0]
                  << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                     "/path/to/partitionFile.cpn\n";
      }
      exit(EXIT_FAILURE);
    }
    OMEGA_H_CHECK(argc == 4);
    const auto clientId = atoi(argv[1]);
    REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
    const auto meshFile = argv[2];
    const auto classPartitionFile = argv[3];
    Omega_h::Mesh mesh(&lib);
    Omega_h::binary::read(meshFile, lib.world(), &mesh);
    MPI_Comm mpi_comm = lib.world()->get_impl();
    switch (clientId) {
      case -1: xgc_coupler(mpi_comm, mesh, classPartitionFile); break;
      case 0: xgc_delta_f(mpi_comm, mesh); break;
      case 1: xgc_total_f(mpi_comm, mesh); break;
      default:
        std::cerr << "Unhandled client id (should be -1, 0,1)\n";
        exit(EXIT_FAILURE);
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception caught in main" << std::endl;
    return 1;
  }
}
