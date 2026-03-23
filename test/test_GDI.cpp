
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include "test_support.h"
#include "pcms/adapter/meshfields/mesh_fields_adapter.h"
#include "pcms/coupler.h"
#include <pcms/utility/types.h>
static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 1;

void xgc_delta_f(MPI_Comm comm)
{
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_delta_f");

  const auto GDI = app->Add_GDI<pcms::GO>("global_comm", comm);
  auto mean = std::vector<long>(1);
  mean[0] = 16;
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      GDI->Send(mean.data(), "mean", mean.size());
      app->EndSendPhase();
      printf("delta Sent mean:%d\n", mean[0]);
      app->BeginReceivePhase();
      mean = GDI->Receive("mean", mean.size());
      app->EndReceivePhase();
      mean[0] = mean[0]/2;
    }
  } while (!done);
  printf("final Mean = %d\n", mean[0]);
  assert(std::fabs(mean[0] - 1.0) < 1e-12);
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
      printf("total Recieved mean:%d\n", mean[0]);
      mean[0] = mean[0]/2;
      app->BeginSendPhase();
      GDI->Send(mean.data(), "mean", mean.size());
      app->EndSendPhase();
      printf("total Sent mean:%d\n", mean[0]);
    }
  } while (!done);
}
void xgc_coupler(MPI_Comm comm)
{
  // Define Partition
  redev::LO dim = 3;
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  auto partition = redev::Partition{redev::RCBPtn{dim, ranks, cuts}};

  pcms::Coupler cpl("proxy_couple", comm, true,
                   partition);
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");

  auto GDI_total = total_f->Add_GDI<pcms::GO>("global_comm", comm);
  auto GDI_delta = delta_f->Add_GDI<pcms::GO>("global_comm", comm);
  auto mean = std::vector<pcms::GO>(1);
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      delta_f->BeginReceivePhase();
      mean = GDI_delta->Receive("mean", 1);
      delta_f->EndReceivePhase();
      printf("delta Received mean:%d\n", mean[0]);
      mean[0] = mean[0]/2;
      const auto  msg_size = mean.size();
      total_f->BeginSendPhase();
      GDI_total->Send(mean.data(), "mean", msg_size);
      total_f->EndSendPhase();
      printf("total sent mean:%d\n", mean[0]);
      total_f->BeginReceivePhase();
      mean = GDI_total->Receive("mean", msg_size);
      total_f->EndReceivePhase();
      printf("delta Received mean:%d\n", mean[0]);
      mean[0] = mean[0]/2;
      delta_f->BeginSendPhase();
      GDI_delta->Send(mean.data(), "mean", msg_size);
      delta_f->EndSendPhase();
      printf("detla sent mean:%d\n", mean[0]);
    }
  } while (!done);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv); // MPI init

  OMEGA_H_CHECK(argc == 2);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);

  int color;
  if (clientId == -1)
    color = 0; // coupler
  else if (clientId == 0)
    color = 1; // client A
  else if (clientId == 1)
    color = 2; // client B
  else
    color = MPI_UNDEFINED;

  MPI_Comm subcomm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &subcomm);

  switch (clientId) {
    case -1: xgc_coupler(subcomm); break;
    case 0: xgc_delta_f(subcomm); break;
    case 1: xgc_total_f(subcomm); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  MPI_Finalize();
  return 0;
}
