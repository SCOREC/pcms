#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include "pcms/coupler2.h"
#include "pcms/create_field.h"
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
  pcms::Coupler2 coupler("proxy_couple", comm, false, {});
  pcms::Application2* app = coupler.AddApplication("proxy_couple_xgc_delta_f");

  auto& layout = app->AddLayout(
    "gids",
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian));

  app->AddField("gids", layout.CreateFieldReal());
  app->AddField("gids2", layout.CreateFieldReal());

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids");  //(Alt) df_gid_field->Send();
      // app->SendField("gids2"); //(Alt) df_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) df_gid_field->Receive();
      // app->ReceiveField("gids2"); //(Alt) df_gid_field->Receive();
      app->EndReceivePhase();
      // cpl.ReceiveField("gids2"); //(Alt) df_gid_field->Receive();
    }
  } while (!done);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  pcms::Coupler2 coupler("proxy_couple", comm, false, {});
  pcms::Application2* app = coupler.AddApplication("proxy_couple_xgc_total_f");

  auto& layout = app->AddLayout(
    "gids",
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian));

  app->AddField("gids", layout.CreateFieldReal());

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
}
void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler2 cpl(
    "proxy_couple", comm, true,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());

  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");

  auto& layout_total = total_f->AddLayout(
    "gids",
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian));

  auto& layout_delta = delta_f->AddLayout(
    "gids",
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian));

  // TODO, fields should have a transfer policy rather than parameters
  total_f->AddField("gids", layout_total.CreateFieldReal());
  delta_f->AddField("gids", layout_delta.CreateFieldReal());
  delta_f->AddField("gids2", layout_delta.CreateFieldReal());

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      total_f->BeginReceivePhase();
      total_f->ReceiveField("gids");
      total_f->EndReceivePhase();

      delta_f->BeginReceivePhase();
      delta_f->ReceiveField("gids");
      // delta_f->ReceiveField("gids2");
      delta_f->EndReceivePhase();

      total_f->BeginSendPhase();
      total_f->SendField("gids");
      total_f->EndSendPhase();

      delta_f->BeginSendPhase();
      delta_f->SendField("gids", pcms::Mode::Deferred);
      // delta_f->SendField("gids2", pcms::Mode::Deferred);
      delta_f->EndSendPhase();
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
