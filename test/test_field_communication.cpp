#include <chrono>
#include <iostream>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <redev.h>
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/field_communicator2.h"
#include "pcms/field_communicator.h"
#include "test_support.h"

namespace ts = test_support;

struct Timer {
  Timer() {
    start_ = std::chrono::high_resolution_clock::now();
  }

  void stop() {
    end_ = std::chrono::high_resolution_clock::now();
  }

  auto elapsed() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end_ - start_).count();
  }

  std::chrono::system_clock::time_point start_;
  std::chrono::system_clock::time_point end_;
};

void client(MPI_Comm comm, redev::Redev& rdv, redev::Channel& channel, Omega_h::Mesh& mesh)
{
  const auto nverts = mesh.nents(0);
  Omega_h::Write<double> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = i; });
  mesh.add_tag<double>(0, "field", 1, Omega_h::Read(ids));

  auto layout = pcms::OmegaHFieldLayout(
    mesh, pcms::OmegaHFieldLayoutLocation::PieceWise, 2);
  pcms::OmegaHField2 new_field("field", pcms::CoordinateSystem::Cartesian, layout,
                          mesh);
  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm("new_comm", comm, rdv, channel, layout);
  pcms::FieldCommunicator2<pcms::Real> new_comm(layout_comm, new_field);

  pcms::OmegaHFieldAdapter<double> old_field("field", mesh);
  pcms::FieldCommunicator<pcms::OmegaHFieldAdapter<double>> old_comm("old_comm", comm, rdv, channel, old_field);

  MPI_Barrier(MPI_COMM_WORLD);
  Timer old_time{};
  channel.BeginSendCommunicationPhase();
  old_comm.Send();
  channel.EndSendCommunicationPhase();
  old_time.stop();

  MPI_Barrier(MPI_COMM_WORLD);
  Timer new_time{};
  channel.BeginSendCommunicationPhase();
  new_comm.Send();
  channel.EndSendCommunicationPhase();
  new_time.stop();

  std::cerr << "Old fields: " << old_time.elapsed() << " ns\n";
  std::cerr << "New fields: " << new_time.elapsed() << " ns\n";
}

void server(MPI_Comm comm, redev::Redev& rdv, redev::Channel& channel, Omega_h::Mesh& mesh) {
  const auto nverts = mesh.nents(0);
  Omega_h::Write<double> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = 0; });
  mesh.add_tag<double>(0, "field", 1, Omega_h::Read(ids));

  auto layout = pcms::OmegaHFieldLayout(
    mesh, pcms::OmegaHFieldLayoutLocation::PieceWise, 2);
  pcms::OmegaHField2 new_field("field", pcms::CoordinateSystem::Cartesian, layout,
                          mesh);
  pcms::FieldLayoutCommunicator<pcms::Real> layout_comm("new_comm", comm, rdv, channel, layout);
  pcms::FieldCommunicator2<pcms::Real> new_comm(layout_comm, new_field);


  pcms::OmegaHFieldAdapter<double> old_field("field", mesh);
  pcms::FieldCommunicator<pcms::OmegaHFieldAdapter<double>> old_comm("old_comm", comm, rdv, channel, old_field);

  MPI_Barrier(MPI_COMM_WORLD);
  channel.BeginReceiveCommunicationPhase();
  old_comm.Receive();
  channel.EndReceiveCommunicationPhase();

  MPI_Barrier(MPI_COMM_WORLD);
  channel.BeginReceiveCommunicationPhase();
  new_comm.Receive();
  channel.EndReceiveCommunicationPhase();

  auto copied_array = new_field.GetDOFHolderData().GetValues();

  int expected = nverts;
  int sum = 0;
  Kokkos::parallel_reduce(
    nverts,
    KOKKOS_LAMBDA(int i, int& local_sum) {
      local_sum += (int) copied_array[i] == i;
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
