#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_for.hpp>
#include <redev.h>
#include <vector>
#include "pcms/coupler/field_communicator.hpp"
#include "pcms/coupler/coupler.hpp"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/field/data/simple.h"
#include "test_support.h"

namespace ts = test_support;
using pcms::Real;

ts::ClassificationPartition getClassificationPartition(Omega_h::Mesh& mesh)
{
  ts::ClassificationPartition cp;
  std::map<redev::ClassPtn::ModelEnt, int> ent_to_rank;

  for (int e = 0; e <= 0; ++e) {
    auto ids_d = mesh.get_array<Omega_h::ClassId>(e, "class_id");
    auto dims_d = mesh.get_array<Omega_h::I8>(e, "class_dim");
    auto remotes = mesh.ask_owners(e);

    auto ids = Omega_h::HostRead(ids_d);
    auto dims = Omega_h::HostRead(dims_d);
    auto remotes_ranks = Omega_h::HostRead(remotes.ranks);

    for (int i = 0; i < mesh.nents(e); ++i) {
      redev::ClassPtn::ModelEnt ent({dims[i], ids[i]});
      int rank = remotes_ranks[i];
      auto it = ent_to_rank.find(ent);
      if (it == ent_to_rank.end()) {
        ent_to_rank[ent] = rank;
        cp.ranks.push_back(rank);
        cp.modelEnts.push_back(ent);
      } else {
        PCMS_ALWAYS_ASSERT(rank == it->second);
      }
    }
  }
  return cp;
}

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName)
{
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank()
                               ? ts::readClassPartitionFile(cpnFileName)
                               : ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  REDEV_ALWAYS_ASSERT(mesh.nelems()); // all ranks should have elements
  auto ptn = getClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

// Test that two fields sharing a layout via AddLayout use the same
// FieldLayoutCommunicator rather than creating separate ones.
static void test_shared_layout(Omega_h::Library& lib,
                               std::string_view mesh_file,
                               std::string_view cpn_file, bool is_server)
{
  auto world = lib.world();
  MPI_Comm mpi_comm = world->get_impl();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(std::string(mesh_file).c_str(), lib.world(), &mesh);

  if (is_server) {
    auto partition = setupServerPartition(mesh, cpn_file);
    pcms::Coupler cpl("shared_layout_server", mpi_comm, true,
                      redev::Partition{partition});
    auto* app = cpl.AddApplication("shared_layout");
    auto factory = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
    auto layout = factory.GetLayout();
    app->AddLayout("shared", layout);
    PCMS_ALWAYS_ASSERT(app->GetLayoutCommunicatorCount() == 1);
    auto f1 = app->AddField("field_a", factory.CreateField<Real>());
    PCMS_ALWAYS_ASSERT(app->GetLayoutCommunicatorCount() ==
                       1); // still 1 after adding field
    auto f2 = app->AddField("field_b", factory.CreateField<Real>());
    PCMS_ALWAYS_ASSERT(app->GetLayoutCommunicatorCount() == 1); // still 1
    app->ReceivePhase([&]() {
      f1.Receive();
      f2.Receive();
    });
    app->SendPhase([&]() {
      f1.Send();
      f2.Send();
    });
  } else {
    pcms::Coupler cpl("shared_layout_client", mpi_comm, false,
                      redev::Partition{});
    auto* app = cpl.AddApplication("shared_layout");
    auto factory = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
    auto layout = factory.GetLayout();
    app->AddLayout("shared", layout);
    auto f1 = app->AddField("field_a", factory.CreateField<Real>());
    auto f2 = app->AddField("field_b", factory.CreateField<Real>());
    app->SendPhase([&]() {
      f1.Send();
      f2.Send();
    });
    app->ReceivePhase([&]() {
      f1.Receive();
      f2.Receive();
    });
  }
}

void client1(MPI_Comm comm, Omega_h::Mesh& mesh, std::string comm_name,
             int order, const adios2::Params& params)
{
  redev::Redev rdv(MPI_COMM_WORLD);
  auto channel =
    rdv.CreateAdiosChannel("field2_chan1", params, redev::TransportType::BP4);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, 1, pcms::CoordinateSystem::Cartesian);
  auto layout = factory.GetLayout();
  auto gids = layout->GetGidsHost();
  const auto n = layout->GetNumOwnedDofHolder();
  Omega_h::HostWrite<Real> ids(n);
  PCMS_ALWAYS_ASSERT(n == gids.size());
  Kokkos::parallel_for(
    "id gid", Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { ids[i] = gids[i]; });

  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  field.SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(ids.data(), n, 1));

  pcms::FieldLayoutCommunicator layout_comm(comm_name + "1", comm, rdv, channel,
                                            *layout);
  pcms::FieldCommunicator<pcms::Real> field_comm(layout_comm.GetName(),
                                                 layout_comm, field);

  channel.BeginSendCommunicationPhase();
  field_comm.Send();
  channel.EndSendCommunicationPhase();
}

void client2(MPI_Comm comm, Omega_h::Mesh& mesh, std::string comm_name,
             int order, const adios2::Params& params)
{
  redev::Redev rdv(MPI_COMM_WORLD);
  auto channel =
    rdv.CreateAdiosChannel("field2_chan2", params, redev::TransportType::BP4);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, 1, pcms::CoordinateSystem::Cartesian);
  auto layout = factory.GetLayout();
  auto gids = layout->GetGidsHost();
  const auto n = layout->GetNumOwnedDofHolder();

  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  pcms::FieldLayoutCommunicator layout_comm(comm_name + "2", comm, rdv, channel,
                                            *layout);
  pcms::FieldCommunicator<pcms::Real> field_comm(layout_comm.GetName(),
                                                 layout_comm, field);

  channel.BeginReceiveCommunicationPhase();
  field_comm.Receive();
  channel.EndReceiveCommunicationPhase();

  auto copied_array = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  auto owned = layout->GetOwnedHost();

  PCMS_ALWAYS_ASSERT(copied_array.size() == gids.size());
  PCMS_ALWAYS_ASSERT(owned.size() == gids.size());

  int expected = 0;
  int sum = 0;
  Kokkos::parallel_reduce(
    "reduce1",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0,
                                                                owned.size()),
    KOKKOS_LAMBDA(int i, int& local_sum) { local_sum += owned[i] != 0; },
    expected);
  Kokkos::parallel_reduce(
    "reduce2",
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(int i, int& local_sum) {
      if (owned[i])
        local_sum += (int)copied_array[i] == gids[i];
    },
    sum);

  if (sum != expected) {
    std::cerr << "Field Communication Failed, expected " << expected << " got "
              << sum << std::endl;
    exit(EXIT_FAILURE);
  } else {
    std::cerr << "Field Communication Passed\n";
  }
}

void server(MPI_Comm comm, Omega_h::Mesh& mesh, std::string comm_name,
            int order, const adios2::Params& params,
            const std::string& cpn_filename)
{
  redev::Redev rdv(MPI_COMM_WORLD,
                   redev::Partition{setupServerPartition(mesh, cpn_filename)},
                   redev::ProcessType::Server);
  auto channel1 =
    rdv.CreateAdiosChannel("field2_chan1", params, redev::TransportType::BP4);
  auto channel2 =
    rdv.CreateAdiosChannel("field2_chan2", params, redev::TransportType::BP4);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, 1, pcms::CoordinateSystem::Cartesian);
  auto layout = factory.GetLayout();
  const auto n = layout->GetNumOwnedDofHolder();
  Omega_h::HostWrite<Real> ids(n);
  Kokkos::parallel_for(
    "id 0", Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, n),
    [=](int i) { ids[i] = 0; });

  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  pcms::FieldLayoutCommunicator layout_comm1(comm_name + "1", comm, rdv,
                                             channel1, *layout);
  pcms::FieldLayoutCommunicator layout_comm2(comm_name + "2", comm, rdv,
                                             channel2, *layout);
  pcms::FieldCommunicator<pcms::Real> field_comm1(layout_comm1.GetName(),
                                                  layout_comm1, field);
  pcms::FieldCommunicator<pcms::Real> field_comm2(layout_comm2.GetName(),
                                                  layout_comm2, field);

  channel1.BeginReceiveCommunicationPhase();
  field_comm1.Receive();
  channel1.EndReceiveCommunicationPhase();

  channel2.BeginSendCommunicationPhase();
  field_comm2.Send();
  channel2.EndSendCommunicationPhase();
}

int main(int argc, char** argv)
{
  try {
    auto lib = Omega_h::Library(&argc, &argv);
    auto world = lib.world();
    int rank = world->rank();
    if (argc != 4) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1|2|3> /path/to/omega_h/mesh"
                << " /path/to/partitionFile.cpn\n";
      exit(EXIT_FAILURE);
    }
    int clientId = atoi(argv[1]);
    REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 3);
    const auto meshFile = argv[2];
    const auto classPartitionFile = argv[3];

    Omega_h::Mesh mesh = Omega_h::binary::read(meshFile, world);
    adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
    MPI_Comm mpi_comm = lib.world()->get_impl();

    switch (clientId) {
      case -1:
        server(mpi_comm, mesh, "lin_field_comm", 1, params, classPartitionFile);
        break;
      case 0: client1(mpi_comm, mesh, "lin_field_comm", 1, params); break;
      case 1: client2(mpi_comm, mesh, "lin_field_comm", 1, params); break;
      case 2:
        test_shared_layout(lib, meshFile, classPartitionFile,
                           /*is_server=*/true);
        break;
      case 3:
        test_shared_layout(lib, meshFile, classPartitionFile,
                           /*is_server=*/false);
        break;
      default:
        std::cerr << "Unhandled client id (should be -1,0,1,2,3)\n";
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
