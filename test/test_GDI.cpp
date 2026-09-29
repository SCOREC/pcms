#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include "test_support.h"
#include "pcms/coupler/coupler.hpp"
#include <pcms/utility/types.h>
#include <vector>

static constexpr bool done = true;
static constexpr int COMM_ROUNDS = 1;

void xgc_delta_f(MPI_Comm comm)
{
  pcms::CouplerComm coupler("proxy_couple", comm, false, {});
  pcms::ApplicationComm* app =
    coupler.AddApplication("proxy_couple_xgc_delta_f");

  std::vector<pcms::GO> mean_storage(1);
  auto mean = pcms::make_array_view(mean_storage);

  app->AddData<pcms::GO>("mean", mean, comm);

  mean[0] = 16;

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendData("mean");
      app->EndSendPhase();
      printf("delta Sent mean:%ld\n", mean[0]);
      app->BeginReceivePhase();
      app->ReceiveData("mean");
      app->EndReceivePhase();
      mean[0] = mean[0] / 2;
    }
  } while (!done);
  printf("final Mean = %ld\n", mean[0]);
  assert(std::fabs(mean[0] - 1.0) < 1e-12);
  printf("GDI test successful.\n");
}
void xgc_total_f(MPI_Comm comm)
{
  pcms::CouplerComm coupler("proxy_couple", comm, false, {});
  pcms::ApplicationComm* app =
    coupler.AddApplication("proxy_couple_xgc_total_f");

  std::vector<pcms::GO> mean_storage(1);
  auto mean = pcms::make_array_view(mean_storage);
  app->AddData<pcms::GO>("mean", mean, comm);

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginReceivePhase();
      app->ReceiveData("mean");
      app->EndReceivePhase();
      printf("total Recieved mean:%ld\n", mean[0]);
      mean[0] = mean[0] / 2;
      app->BeginSendPhase();
      app->SendData("mean");
      app->EndSendPhase();
      printf("total Sent mean:%ld\n", mean[0]);
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

  pcms::CouplerComm cpl("proxy_couple", comm, true, partition);
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f");
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f");

  std::vector<pcms::GO> mean_storage(1);
  auto mean = pcms::make_array_view(mean_storage);
  total_f->AddData<pcms::GO>("mean", mean, comm);
  delta_f->AddData<pcms::GO>("mean", mean, comm);

  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      delta_f->BeginReceivePhase();
      delta_f->ReceiveData("mean");
      delta_f->EndReceivePhase();
      printf("delta Received mean:%ld\n", mean[0]);
      mean[0] = mean[0] / 2;
      total_f->BeginSendPhase();
      total_f->SendData("mean");
      total_f->EndSendPhase();
      printf("total sent mean:%ld\n", mean[0]);
      total_f->BeginReceivePhase();
      total_f->ReceiveData("mean");
      total_f->EndReceivePhase();
      printf("delta Received mean:%ld\n", mean[0]);
      mean[0] = mean[0] / 2;
      delta_f->BeginSendPhase();
      delta_f->SendData("mean");
      delta_f->EndSendPhase();
      printf("delta sent mean:%ld\n", mean[0]);
    }
  } while (!done);
}

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    OMEGA_H_CHECK(argc == 2);

    const auto clientId = std::atoi(argv[1]);
    REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);

    MPI_Comm comm = MPI_COMM_WORLD;

    switch (clientId) {
      case -1: xgc_coupler(comm); break;

      case 0: xgc_delta_f(comm); break;

      case 1: xgc_total_f(comm); break;
      default:
        std::cerr << "Unhandled client id; expected -1, 0, or 1\n";
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception caught in main\n";
    return 1;
  }
}
