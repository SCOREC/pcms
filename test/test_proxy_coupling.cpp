#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/print.h"
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
static constexpr int COMM_ROUNDS = 10;
static std::string sstDataTransport;
static std::string adiosEngine;
namespace ts = test_support;

auto getAdiosEngine() {
  if( adiosEngine == "SST" )
    return redev::TransportType::SST;
  else if( adiosEngine == "BP4" )
    return redev::TransportType::BP4;
  else
    exit(EXIT_FAILURE);
}

adios2::Params getAdiosParams(const redev::TransportType engine) {
  if( engine != redev::TransportType::SST ) {
    return {{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  }

  adios2::Params params;
  if( sstDataTransport == "MPI" )
    params = {{"DataTransport", "MPI"}, {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  else if( sstDataTransport == "RDMA" )
    params = {{"DataTransport", "RDMA"}, {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  else if( sstDataTransport == "WAN" )
    params = {{"DataTransport", "WAN"}, {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  else
    exit(EXIT_FAILURE);
  return params;
}

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_delta_f", "", adiosEngine, adiosParams);
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  app->AddField("gids2",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  PCMS_FUNCTION_TIMER
  const auto start{std::chrono::steady_clock::now()};
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
      if(!rank) fprintf(stderr, "round %d is done\n", i);
    }
  } while (!done);
  MPI_Barrier(comm);
  const auto finish{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{finish - start};
  if(!rank) std::cerr << "xgc_delta_f " << elapsed_seconds.count() << "\n";
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  pcms::Coupler coupler("proxy_couple", comm, false, {});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  pcms::Application* app = coupler.AddApplication("proxy_couple_xgc_total_f", "", adiosEngine, adiosParams);
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{});
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("global", mesh, is_overlap));
  const auto numOverlapVerts = Omega_h::get_sum(is_overlap);
  const auto hasOverlapVerts = (numOverlapVerts > 0) ? 1 : 0;
  const auto numRanksWithOverlapVerts = mesh.comm()->allreduce(hasOverlapVerts, OMEGA_H_SUM);
  if(rank == 0) {
    pcms::printInfo("numOverlapVerts %d numRanksWithOverlapVerts %d\n", numOverlapVerts, numRanksWithOverlapVerts);
  }
  PCMS_FUNCTION_TIMER
  const auto start{std::chrono::steady_clock::now()};
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      app->BeginSendPhase();
      app->SendField("gids"); //(Alt) tf_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) tf_gid_field->Receive();
      app->EndReceivePhase();
      if(!rank) fprintf(stderr, "round %d is done\n", i);
    }
  } while (!done);
  MPI_Barrier(comm);
  const auto finish{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{finish - start};
  if(!rank) std::cerr << "xgc_total_f " << elapsed_seconds.count() << "\n";
}
void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  //redirect stdout to one file per rank
  auto fname = std::string("coupler_p") + std::to_string(rank) + ".log";
  FILE* fhandle = fopen(fname.c_str(), "w");
  pcms::setStdout(fhandle);
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl(
    "proxy_couple", comm, true,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  pcms::printInfo("mesh numVtx %d\n", mesh.nverts());
  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  auto* total_f = cpl.AddApplication("proxy_couple_xgc_total_f", "", adiosEngine, adiosParams);
  auto* delta_f = cpl.AddApplication("proxy_couple_xgc_delta_f", "", adiosEngine, adiosParams);
  // TODO, fields should have a transfer policy rather than parameters
  auto* total_f_gids = total_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("total_f_gids", mesh, is_overlap));
  auto* delta_f_gids = delta_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("delta_f_gids", mesh, is_overlap));
  auto* delta_f_gids2 = delta_f->AddField(
    "gids2", OmegaHFieldAdapter<GO>("delta_f_gids2", mesh, is_overlap));
  const auto numServerOverlapVerts = Omega_h::get_sum(is_overlap);

  pcms::printInfo("numServerOverlapVerts %d\n", numServerOverlapVerts);
  pcms::printInfo("round, total_f, delta_f_gids, delta_f_gids2, local_total, global_total\n");
  {
  PCMS_FUNCTION_TIMER
  auto start{std::chrono::steady_clock::now()};
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      total_f->ReceivePhase([&]() { total_f_gids->Receive(); });
      delta_f->ReceivePhase([&]() {
        delta_f_gids->Receive();
        delta_f_gids2->Receive();
      });
      // Get bytes received after receive phase, don't include in timing, only
      // need one round of data
      if(i == 0) {
        const auto getBytes_start{std::chrono::steady_clock::now()};
        size_t total_f_bytes = total_f_gids->GetBytesReceived();
        size_t delta_f_gids_bytes = delta_f_gids->GetBytesReceived();
        size_t delta_f_gids2_bytes = delta_f_gids2->GetBytesReceived();
        size_t total_bytes = total_f_bytes + delta_f_gids_bytes + delta_f_gids2_bytes;
        size_t global_total_bytes = redev::GetTotalBytesReceived(total_bytes, comm);
        pcms::printInfo("%d, %zu, %zu, %zu, %zu, %zu\n",
                  i, total_f_bytes, delta_f_gids_bytes,
                  delta_f_gids2_bytes, total_bytes, global_total_bytes);
        start -= std::chrono::steady_clock::now() - getBytes_start;
      }
      total_f->SendPhase([&]() { total_f_gids->Send(); });
      delta_f->SendPhase([&]() {
        delta_f_gids->Send(pcms::Mode::Deferred);
        delta_f_gids2->Send(pcms::Mode::Deferred);
      });
    }
  } while (!done);
  MPI_Barrier(comm);
  const auto finish{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{finish - start};
  if(!rank) std::cerr << "xgc_coupler " << elapsed_seconds.count() << "\n";
  }
  fclose(fhandle);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}

int main(int argc, char** argv)
{
  int provide;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provide);
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(!rank) std::cerr << "mpi thread level: " << provide << "\n";

  if (argc != 6) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                   "/path/to/partitionFile.cpn "
                   "sstDataTransport=[RDMA|WAN|MPI] "
                   "adiosEngine=[BP4|SST]\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 6);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  sstDataTransport = argv[4];
  adiosEngine = argv[5];
  if(!rank) {
    std::cerr << "inputs: " << clientId << ", " << meshFile << " " << classPartitionFile << " "
              << sstDataTransport << " " << adiosEngine << "\n";
  }
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
  if(!rank) std::cerr << "done " << "\n";
  return 0;
}
