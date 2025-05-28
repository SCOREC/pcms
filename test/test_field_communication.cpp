#include <iostream>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <redev.h>
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/field_communicator2.h"
#include "test_support.h"

namespace ts = test_support;

void client(MPI_Comm comm, redev::Redev& rdv, redev::Channel& channel, Omega_h::Mesh& mesh)
{
  const auto nverts = mesh.nents(0);
  Omega_h::Write<double> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = i; });
  mesh.add_tag<double>(0, "test", 1, Omega_h::Read(ids));

  auto layout = pcms::OmegaHFieldLayout(
    mesh, pcms::OmegaHFieldLayoutLocation::PieceWise, 2);
  pcms::OmegaHField2 test("test", pcms::CoordinateSystem::Cartesian, layout,
                          mesh);
  pcms::FieldCommunicator2<pcms::Real> field_comm("test_comm", comm, rdv, channel, layout);
  channel.BeginSendCommunicationPhase();
  field_comm.Send(test);
  channel.EndSendCommunicationPhase();
}

void server(MPI_Comm comm, redev::Redev& rdv, redev::Channel& channel, Omega_h::Mesh& mesh) {
  const auto nverts = mesh.nents(0);
  Omega_h::Write<double> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = 0; });
  mesh.add_tag<double>(0, "test", 1, Omega_h::Read(ids));

  auto layout = pcms::OmegaHFieldLayout(
    mesh, pcms::OmegaHFieldLayoutLocation::PieceWise, 2);
  pcms::OmegaHField2 test("test", pcms::CoordinateSystem::Cartesian, layout,
                          mesh);
  pcms::FieldCommunicator2<pcms::Real> field_comm("test_comm", comm, rdv, channel, layout);
  channel.BeginReceiveCommunicationPhase();
  field_comm.Receive(test);
  channel.EndReceiveCommunicationPhase();

  auto copied_array = test.GetDOFHolderData().GetValues();

  int expected = (nverts * (nverts - 1)) / 2;
  int sum = 0;
  Kokkos::parallel_reduce(
    nverts,
    KOKKOS_LAMBDA(int i, int& local_sum) {
      local_sum += (int) copied_array[i];
    },
    sum);

  if (sum != expected) {
    std::cerr << "Field Communication Failed, expected " << expected << " got "
              << sum << std::endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 2) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=0|1>";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 2);
  const auto clientId = atoi(argv[1]);

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 100, 100, 0, false);
  bool isRdv = clientId == 0;
  const auto classPartition =
    isRdv ? ts::CreateClassificationPartition(mesh) : ts::ClassificationPartition();
  MPI_Comm mpi_comm = lib.world()->get_impl();
  auto partition = redev::ClassPtn(MPI_COMM_WORLD, classPartition.ranks,
                                   classPartition.modelEnts);
  redev::Redev rdv(MPI_COMM_WORLD, redev::Partition{std::move(partition)},
                   static_cast<redev::ProcessType>(isRdv));
  adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  const std::string name = "test_ids";
  auto channel =
    rdv.CreateAdiosChannel(name, params, redev::TransportType::BP4);

  switch (clientId) {
    case 0: server(mpi_comm, rdv, channel, mesh); break;
    case 1: client(mpi_comm, rdv, channel, mesh); break;
    default:
      std::cerr << "Unhandled client id (should be 0,1)\n";
      exit(EXIT_FAILURE);
  }

  return 0;
}
