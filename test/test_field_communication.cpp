#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_for.hpp>
#include <redev.h>
#include <vector>
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/field_communicator2.h"
#include "pcms/field_communicator.h"
#include "pcms/create_field.h"
#include "test_support.h"

namespace ts = test_support;

ts::ClassificationPartition getClassificationPartition(Omega_h::Mesh& mesh) {
  ts::ClassificationPartition cp;
  std::map<redev::ClassPtn::ModelEnt, int> ent_to_rank;

  for (int e = 0; e <= 0; ++e) {
    auto ids = mesh.get_array<Omega_h::ClassId>(e, "class_id");
    auto dims = mesh.get_array<Omega_h::I8>(e, "class_dim");
    auto remotes = mesh.ask_owners(e);

    for (int i = 0; i < mesh.nents(e); ++i) {
      redev::ClassPtn::ModelEnt ent({dims[i], ids[i]});
      int rank = remotes.ranks[i];
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

void client1(MPI_Comm comm, Omega_h::Mesh& mesh, std::string comm_name,
             int order, const adios2::Params& params)
{
  redev::Redev rdv(MPI_COMM_WORLD);
  auto channel =
    rdv.CreateAdiosChannel("field2_chan1", params, redev::TransportType::BP4);

  auto layout = pcms::CreateLagrangeLayout(mesh, order, 1,
                                           pcms::CoordinateSystem::Cartesian);
  auto gids = layout->GetGids();
  const auto n = layout->GetNumOwnedDofHolder();
  Omega_h::Write<double> ids(n);
  PCMS_ALWAYS_ASSERT(n == gids.size());
  Omega_h::parallel_for(
    n, OMEGA_H_LAMBDA(int i) { ids[i] = gids[i]; });

  auto field = pcms::CreateField("", *layout);
  pcms::Rank1View<const double, pcms::HostMemorySpace> array_view{
    std::data(ids), std::size(ids)};
  pcms::FieldDataView<const double, pcms::HostMemorySpace> field_data_view{
    array_view, field->GetCoordinateSystem()};
  field->SetDOFHolderData(field_data_view);

  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm(comm_name + "1", comm, rdv, channel, *layout);
  pcms::FieldCommunicator2<pcms::Real> field_comm(layout_comm, *field);

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

  auto layout = pcms::CreateLagrangeLayout(mesh, order, 1,
                                           pcms::CoordinateSystem::Cartesian);
  auto gids = layout->GetGids();
  const auto n = layout->GetNumOwnedDofHolder();

  auto field = pcms::CreateField("", *layout);
  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm(comm_name + "2", comm, rdv, channel, *layout);
  pcms::FieldCommunicator2<pcms::Real> field_comm(layout_comm, *field);

  channel.BeginReceiveCommunicationPhase();
  field_comm.Receive();
  channel.EndReceiveCommunicationPhase();

  auto copied_array = field->GetDOFHolderData().GetValues();
  auto owned = layout->GetOwned();

  PCMS_ALWAYS_ASSERT(copied_array.size() == gids.size());
  PCMS_ALWAYS_ASSERT(owned.size() == gids.size());

  int expected = 0;
  int sum = 0;
  Kokkos::parallel_reduce(
    owned.size(),
    KOKKOS_LAMBDA(int i, int& local_sum) {
      local_sum += owned[i] != 0;
    },
    expected);
  Kokkos::parallel_reduce(
    n,
    KOKKOS_LAMBDA(int i, int& local_sum) {
      if (owned[i])
        local_sum += (int) copied_array[i] == gids[i];
    },
    sum);

  if (sum != expected) {
    std::cerr << "Field Communication Failed, expected " << expected << " got "
              << sum << std::endl;
    exit(EXIT_FAILURE);
  }
  else {
    std::cerr << "Field Communication Passed\n";
  }
}

void server2(MPI_Comm comm, Omega_h::Mesh& mesh, std::string comm_name,
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

  auto layout = pcms::CreateLagrangeLayout(mesh, order, 1,
                                           pcms::CoordinateSystem::Cartesian);
  const auto n = layout->GetNumOwnedDofHolder();
  Omega_h::Write<double> ids(n);
  Omega_h::parallel_for(
    n, OMEGA_H_LAMBDA(int i) { ids[i] = 0; });

  auto field = pcms::CreateField("", *layout);
  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm1(comm_name + "1", comm, rdv, channel1, *layout);
  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm2(comm_name + "2", comm, rdv, channel2, *layout);
  pcms::FieldCommunicator2<pcms::Real> field_comm1(layout_comm1, *field);
  pcms::FieldCommunicator2<pcms::Real> field_comm2(layout_comm2, *field);

  channel1.BeginReceiveCommunicationPhase();
  field_comm1.Receive();
  channel1.EndReceiveCommunicationPhase();

  channel2.BeginSendCommunicationPhase();
  field_comm2.Send();
  channel2.EndSendCommunicationPhase();
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  int rank = world->rank();
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
              << " <clientId=-1|0|1> /path/to/omega_h/mesh"
              << "/path/to/partitionFile.cpn\n";
    exit(EXIT_FAILURE);
  }
  int clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];

  Omega_h::Mesh mesh;
  adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  MPI_Comm mpi_comm = lib.world()->get_impl();

  switch (clientId) {
    case -1:
      mesh = Omega_h::binary::read(meshFile, world);
      server2(mpi_comm, mesh, "lin_field_comm", 1, params, classPartitionFile);
      break;
    case 0:
      mesh = Omega_h::binary::read(meshFile, world);
      client1(mpi_comm, mesh, "lin_field_comm", 1, params);
      break;
    case 1:
      mesh = Omega_h::binary::read(meshFile, world);
      client2(mpi_comm, mesh, "lin_field_comm", 1, params);
      break;
    default:
      std::cerr << "Unhandled client id (should be -1,0,1)\n";
      exit(EXIT_FAILURE);
  }

  return 0;
}
