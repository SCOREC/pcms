#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include "test_support.h"
#include "pcms/adapter/omega_h/omega_h_field.h"

static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 4;

void xgc_delta_f(MPI_Comm comm)
{
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_delta_f");

  const auto GDI = app->Add_GDI<pcms::GO>("global_comm", comm);
  auto mean = std::vector<pcms::GO>(1);
  mean[0] = 16;
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      GDI->Send(mean.data(), "mean", mean.size());
      app->EndSendPhase();
      app->BeginReceivePhase();
      mean = GDI->Receive("mean", mean.size());
      app->EndReceivePhase();
      mean[0] = mean[0]/2;
    }
  } while (!done);
  assert(mean[0]==1);
  printf("GDI test successful.\n");
}
void xgc_total_f(MPI_Comm comm)
{
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_total_f");

  auto GDI = app->Add_GDI<pcms::GO>("global_comm", comm);
  auto mean = std::vector<pcms::GO>(1);
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginReceivePhase();
      mean = GDI->Receive("mean", mean.size());
      app->EndReceivePhase();
      mean[0] = mean[0]/2;
      app->BeginSendPhase();
      GDI->Send(mean.data(), "mean", mean.size());
      app->EndSendPhase();
    }
  } while (!done);
}
void xgc_coupler(MPI_Comm comm)
{
  pcms::Coupler cpl("proxy_couple", comm, true,
                    {});
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");

  auto GDI_f = total_f->Add_GDI<pcms::GO>("global_comm", comm);
  auto GDI_d = delta_f->Add_GDI<pcms::GO>("global_comm", comm);
  auto mean = std::vector<pcms::GO>(1);
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      total_f->BeginReceivePhase();
      mean = GDI_f->Receive("mean", 1);
      total_f->EndReceivePhase();
      mean[0] = mean[0]/2;
      const auto  msg_size = mean.size();
      delta_f->BeginSendPhase();
      GDI_d->Send(mean.data(), "mean", msg_size);
      delta_f->EndSendPhase();
      delta_f->BeginReceivePhase();
      mean = GDI_d->Receive("mean", msg_size);
      delta_f->EndReceivePhase();
      mean[0] = mean[0]/2;
      total_f->BeginSendPhase();
      GDI_f->Send(mean.data(), "mean", msg_size);
      total_f->EndSendPhase();
    }
  } while (!done);
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 2) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh ";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 2);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);

  Omega_h::Mesh mesh(&lib);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  const std::string name = "meshVtxIds";
  switch (clientId) {
    case -1: xgc_coupler(mpi_comm); break;
    case 0: xgc_delta_f(mpi_comm); break;
    case 1: xgc_total_f(mpi_comm); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  return 0;
}
