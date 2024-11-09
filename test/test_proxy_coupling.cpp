#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include <pcms/omega_h_field.h>
#include <chrono>
#include <thread>

using pcms::Copy;
using pcms::CouplerClient;
using pcms::CouplerServer;
using pcms::FieldEvaluationMethod;
using pcms::FieldTransferMethod;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::OmegaHField;
using pcms::OmegaHFieldAdapter;

using namespace std::chrono_literals;

static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 4;
namespace ts = test_support;

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  CouplerClient cpl("proxy_couple_xgc_delta_f", comm);

  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  cpl.AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  cpl.AddField("gids2",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      cpl.BeginSendPhase();
      cpl.SendField("gids");  //(Alt) df_gid_field->Send();
      cpl.SendField("gids2"); //(Alt) df_gid_field->Send();
      cpl.EndSendPhase();
      cpl.BeginReceivePhase();
      cpl.ReceiveField("gids"); //(Alt) df_gid_field->Receive();
      cpl.EndReceivePhase();
      // cpl.ReceiveField("gids2"); //(Alt) df_gid_field->Receive();
    }
  } while (!done);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  pcms::CouplerClient cpl("proxy_couple_xgc_total_f", comm);
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  cpl.AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      cpl.BeginSendPhase();
      cpl.SendField("gids"); //(Alt) tf_gid_field->Send();
      cpl.EndSendPhase();
      cpl.BeginReceivePhase();
      cpl.ReceiveField("gids"); //(Alt) tf_gid_field->Receive();
      cpl.EndReceivePhase();
    }
  } while (!done);
}
void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::CouplerServer cpl(
    "proxy_couple", comm,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)}, mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");
  // TODO, fields should have a transfer policy rather than parameters
  auto* total_f_gids = total_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("total_f_gids", mesh, is_overlap),
    FieldTransferMethod::Copy, // to Omega_h
    FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, // from Omega_h
    FieldEvaluationMethod::None, is_overlap);
  auto* delta_f_gids = delta_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("delta_f_gids", mesh, is_overlap),
    FieldTransferMethod::Copy, FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, FieldEvaluationMethod::None, is_overlap);
  auto* delta_f_gids2 = delta_f->AddField(
    "gids2", OmegaHFieldAdapter<GO>("delta_f_gids2", mesh, is_overlap),
    FieldTransferMethod::Copy, FieldEvaluationMethod::None,
    FieldTransferMethod::Copy, FieldEvaluationMethod::None, is_overlap);
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      //  Gather OHField
      // 1. receives any member fields .Receive()
      // 2. field_transfer native to internal
      // 3. combine internal fields into combined internal field
      total_f->ReceivePhase([&]() { total_f_gids->Receive(); });
      delta_f->ReceivePhase([&]() {
        delta_f_gids->Receive();
        delta_f_gids2->Receive();
      });
      // Scatter OHField
      // 1. OHField transfer internal to native
      // 2. Send data to members
      total_f->SendPhase([&]() { total_f_gids->Send(); });
      delta_f->SendPhase([&]() {
        delta_f_gids->Send(pcms::Mode::Deferred);
        delta_f_gids2->Send(pcms::Mode::Deferred);
      });
    }
  } while (!done);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}

int main(int argc, char** argv)
{
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
  const std::string name = "meshVtxIds";
  switch (clientId) {
    case -1: xgc_coupler(mpi_comm, mesh, classPartitionFile); break;
    case 0: xgc_delta_f(mpi_comm, mesh); break;
    case 1: xgc_total_f(mpi_comm, mesh); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  return 0;
}
