#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include "pcms/adapter/omega_h/omega_h_field.h"
#include <chrono>
#include <thread>

using pcms::Copy;
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
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_delta_f");

  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  app->AddField("gids2",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  PCMS_FUNCTION_TIMER
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids");  //(Alt) df_gid_field->Send();
      app->SendField("gids2"); //(Alt) df_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) df_gid_field->Receive();
      app->EndReceivePhase();
      // cpl.ReceiveField("gids2"); //(Alt) df_gid_field->Receive();
    }
  } while (!done);
  MPI_Barrier(comm);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_total_f");
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  PCMS_FUNCTION_TIMER
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids"); //(Alt) tf_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) tf_gid_field->Receive();
      app->EndReceivePhase();
    }
  } while (!done);
  MPI_Barrier(comm);
}
void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  PCMS_FUNCTION_TIMER
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl(
    "proxy_couple", comm, true,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");
  // TODO, fields should have a transfer policy rather than parameters
  auto* total_f_gids = total_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("total_f_gids", mesh, is_overlap));
  auto* delta_f_gids = delta_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("delta_f_gids", mesh, is_overlap));
  auto* delta_f_gids2 = delta_f->AddField(
    "gids2", OmegaHFieldAdapter<GO>("delta_f_gids2", mesh, is_overlap));
  {
  PCMS_FUNCTION_TIMER
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      total_f->ReceivePhase([&]() { total_f_gids->Receive(); });
      delta_f->ReceivePhase([&]() {
        delta_f_gids->Receive();
        delta_f_gids2->Receive();
      });
      total_f->SendPhase([&]() { total_f_gids->Send(); });
      delta_f->SendPhase([&]() {
        delta_f_gids->Send(pcms::Mode::Deferred);
        delta_f_gids2->Send(pcms::Mode::Deferred);
      });
    }
  } while (!done);
  MPI_Barrier(comm);
  }
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
