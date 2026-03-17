#include <Omega_h_mesh.hpp>
#include <iostream>
#include <pcms.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/utility/print.h"
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
static std::string overlapSize;
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
    params = {{"DataTransport", "RDMA"}, {"Streaming", "On"}, {"OpenTimeoutSecs", "360"}};
  else if( sstDataTransport == "WAN" )
    params = {{"DataTransport", "WAN"}, {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
  else
    exit(EXIT_FAILURE);
  return params;
}

void validate_received_gids(const std::string& prefix, const std::string& field_name, Omega_h::Mesh& mesh,
                           const Omega_h::Read<Omega_h::I8>& is_overlap, MPI_Comm comm)
{
  PERFSTUBS_SCOPED_TIMER("validate_received_gids");
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get the received field data (which contains gids)
  auto received_gids = mesh.get_array<GO>(0, field_name);

  // Get the mesh's vertex globals
  auto mesh_globals = mesh.globals(0);

  OMEGA_H_CHECK_OP(received_gids.size(),==, mesh_globals.size());
  if (received_gids.size() != mesh_globals.size()) {
    std::cerr << "ERROR: rank " << rank << " received_gids.size() != mesh_globals.size()\n";
  }

  // Check for mismatches on device using parallel_for
  const auto n = received_gids.size();
  Omega_h::Write<Omega_h::I8> mismatch_flags(n, 0);

  Omega_h::parallel_for(n, OMEGA_H_LAMBDA(Omega_h::LO i) {
    // Only check vertices that are in the overlap region (where data is received)
    if (is_overlap[i]) {
      if (received_gids[i] != mesh_globals[i]) {
        mismatch_flags[i] = 1;
      }
    }
  });

  // Count total mismatches
  const int mismatches = Omega_h::get_sum(Omega_h::Read<Omega_h::I8>(mismatch_flags));

  if (mismatches > 0) {
    mesh.add_tag(Omega_h::VERT, "recvGidDoesNotMatch", 1, Omega_h::read(mismatch_flags));
    Omega_h::vtk::write_parallel(prefix + field_name + std::string(".vtk"), &mesh, mesh.dim());

    auto received_gids_h = Omega_h::HostRead<GO>(received_gids);
    auto mesh_globals_h = Omega_h::HostRead<GO>(mesh_globals);
    auto mismatch_flags_h = Omega_h::HostRead<Omega_h::I8>(mismatch_flags);

    int printed = 0;
    for (int i = 0; i < n && printed < 10; ++i) {
      if (mismatch_flags_h[i]) {
        std::stringstream ss;
        ss << "Rank " << rank << " field '" << field_name
           << "' mismatch at vertex " << i
           << ": received=" << received_gids_h[i]
           << " expected=" << mesh_globals_h[i] << "\n";
        std::string str = ss.str();
        pcms::printInfo("%s", str.c_str());
        printed++;
      }
    }
  }

  int global_mismatches;
  MPI_Allreduce(&mismatches, &global_mismatches, 1, MPI_INT, MPI_SUM, comm);

  if (rank == 0 && global_mismatches) {
    std::cerr << "Field " << field_name << " validation FAILED: "
              << global_mismatches << " total mismatches\n";
  }
  if (global_mismatches) {
    MPI_Abort(comm, 1);
  }
}

Omega_h::Read<GO> createGlobalsCopy(Omega_h::Mesh& mesh) {
  Omega_h::Write<GO> dup(mesh.nverts());
  Omega_h::copy_into(mesh.globals(Omega_h::VERT), dup);
  return Omega_h::read(dup);
}

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  pcms::Coupler coupler("coupler", comm, false, {});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  pcms::Application* app = coupler.AddApplication("coupler_xgc_delta_f", "", adiosEngine, adiosParams);
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{overlapSize});
  auto deltaf_gids_r = createGlobalsCopy(mesh);
  mesh.add_tag(Omega_h::VERT, "deltaf_gids", 1, deltaf_gids_r);
  auto deltaf_gids2_r = createGlobalsCopy(mesh);
  mesh.add_tag(Omega_h::VERT, "deltaf_gids2", 1, deltaf_gids2_r);
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("deltaf_gids", mesh, is_overlap));
  app->AddField("gids2",
               OmegaHFieldAdapter<GO>("deltaf_gids2", mesh, is_overlap));

  const auto numOverlapVerts = Omega_h::get_sum(is_overlap);
  const auto hasOverlapVerts = (numOverlapVerts > 0) ? 1 : 0;
  const auto numGlobalOverlapVerts = mesh.comm()->allreduce(numOverlapVerts, OMEGA_H_SUM);
  const auto numRanksWithOverlapVerts = mesh.comm()->allreduce(hasOverlapVerts, OMEGA_H_SUM);
  if(rank == 0) {
    pcms::printInfo("numGlobalOverlapVerts %d numRanksWithOverlapVerts %d\n", numGlobalOverlapVerts, numRanksWithOverlapVerts);
  }

  Omega_h::vtk::write_parallel("xgc_delta_f_init.vtk", &mesh, mesh.dim());
  PCMS_FUNCTION_TIMER
  auto start{std::chrono::steady_clock::now()};
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      auto round_start{std::chrono::steady_clock::now()};
      app->BeginSendPhase();
      app->SendField("gids");  //(Alt) df_gid_field->Send();
      app->SendField("gids2"); //(Alt) df_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) df_gid_field->Receive();
      app->ReceiveField("gids2"); //(Alt) df_gid_field->Receive();
      app->EndReceivePhase();

      const auto validate_start{std::chrono::steady_clock::now()};
      validate_received_gids("deltaf_validate_", "deltaf_gids", mesh, is_overlap, comm);
      validate_received_gids("deltaf_validate_", "deltaf_gids2", mesh, is_overlap, comm);
      const auto validate_elapsed = std::chrono::steady_clock::now() - validate_start;
      start -= validate_elapsed;
      round_start -= validate_elapsed;

      const auto round_finish{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> round_elapsed_seconds{round_finish - round_start};
      if(!rank) std::cerr << "round " << i << " done in " << round_elapsed_seconds.count() << " seconds\n";
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
  pcms::Coupler coupler("coupler", comm, false, {});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  pcms::Application* app = coupler.AddApplication("coupler_xgc_total_f", "", adiosEngine, adiosParams);
  auto is_overlap = ts::markOverlapMeshEntities(mesh, ts::IsModelEntInOverlap{overlapSize});
  auto totalf_gids_r = createGlobalsCopy(mesh);
  mesh.add_tag(Omega_h::VERT, "totalf_gids", 1, totalf_gids_r);
  app->AddField("gids",
               OmegaHFieldAdapter<GO>("totalf_gids", mesh, is_overlap));
  const auto numOverlapVerts = Omega_h::get_sum(is_overlap);
  const auto hasOverlapVerts = (numOverlapVerts > 0) ? 1 : 0;
  const auto numGlobalOverlapVerts = mesh.comm()->allreduce(numOverlapVerts, OMEGA_H_SUM);
  const auto numRanksWithOverlapVerts = mesh.comm()->allreduce(hasOverlapVerts, OMEGA_H_SUM);
  pcms::printInfo("numGlobalOverlapVerts %d numRanksWithOverlapVerts %d numLocalOverlapVerts %d\n", numGlobalOverlapVerts, numRanksWithOverlapVerts, numOverlapVerts);

  Omega_h::vtk::write_parallel("xgc_total_f_init.vtk", &mesh, mesh.dim());
  PCMS_FUNCTION_TIMER
  auto start{std::chrono::steady_clock::now()};
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      auto round_start{std::chrono::steady_clock::now()};
      app->BeginSendPhase();
      app->SendField("gids"); //(Alt) tf_gid_field->Send();
      app->EndSendPhase();
      app->BeginReceivePhase();
      app->ReceiveField("gids"); //(Alt) tf_gid_field->Receive();
      app->EndReceivePhase();

      const auto validate_start{std::chrono::steady_clock::now()};
      validate_received_gids("totalf_validate_", "totalf_gids", mesh, is_overlap, comm);
      const auto validate_elapsed = std::chrono::steady_clock::now() - validate_start;
      start -= validate_elapsed;
      round_start -= validate_elapsed;

      const auto round_finish{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> round_elapsed_seconds{round_finish - round_start};
      if(!rank) std::cerr << "round " << i << " done in " << round_elapsed_seconds.count() << " seconds\n";
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

  auto setup_start{std::chrono::steady_clock::now()};
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl(
    "coupler", comm, true,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  pcms::printInfo("mesh numVtx %d\n", mesh.nverts());
  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{overlapSize});
  const auto adiosEngine = getAdiosEngine();
  const auto adiosParams = getAdiosParams(adiosEngine);
  auto* total_f = cpl.AddApplication("coupler_xgc_total_f", "", adiosEngine, adiosParams);
  auto* delta_f = cpl.AddApplication("coupler_xgc_delta_f", "", adiosEngine, adiosParams);
  // TODO, fields should have a transfer policy rather than parameters
  auto* total_f_gids = total_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("total_f_gids", mesh, is_overlap));
  auto* delta_f_gids = delta_f->AddField(
    "gids", OmegaHFieldAdapter<GO>("delta_f_gids", mesh, is_overlap));
  auto* delta_f_gids2 = delta_f->AddField(
    "gids2", OmegaHFieldAdapter<GO>("delta_f_gids2", mesh, is_overlap));
  const auto numServerOverlapVerts = Omega_h::get_sum(is_overlap);
  const auto setup_finish{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> setup_elapsed_seconds{setup_finish - setup_start};
  MPI_Barrier(comm);
  if(!rank) std::cerr << "setup done in " << setup_elapsed_seconds.count() << " seconds\n";

  pcms::printInfo("numServerOverlapVerts %d\n", numServerOverlapVerts);
  pcms::printInfo("round, total_f, delta_f_gids, delta_f_gids2, local_total, global_total\n");
  Omega_h::vtk::write_parallel("xgc_coupler_init.vtk", &mesh, mesh.dim());
  {
  PCMS_FUNCTION_TIMER
  auto start{std::chrono::steady_clock::now()};
  do {
    for (int i = 0; i < COMM_ROUNDS; ++i) {
      auto round_start{std::chrono::steady_clock::now()};
      std::string timerName = std::string("CommRound") + std::to_string(i);
      PERFSTUBS_SCOPED_TIMER(timerName.c_str());
      total_f->ReceivePhase([&]() { total_f_gids->Receive(); });
      delta_f->ReceivePhase([&]() {
        delta_f_gids->Receive();
        delta_f_gids2->Receive();
      });

      const auto validate_start{std::chrono::steady_clock::now()};
      validate_received_gids("coupler_validate_", "total_f_gids", mesh, is_overlap, comm);
      validate_received_gids("coupler_validate_", "delta_f_gids", mesh, is_overlap, comm);
      validate_received_gids("coupler_validate_", "delta_f_gids2", mesh, is_overlap, comm);
      auto validate_elapsed = std::chrono::steady_clock::now() - validate_start;
      start -= validate_elapsed;
      round_start -= validate_elapsed;

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

      std::string sendTotfTimerName = std::string("SendTotalf") + std::to_string(i);
      {
      PERFSTUBS_SCOPED_TIMER(sendTotfTimerName.c_str());
      total_f->SendPhase([&]() { total_f_gids->Send(); });
      }

      std::string sendDelfTimerName = std::string("SendDeltaf") + std::to_string(i);
      {
      PERFSTUBS_SCOPED_TIMER(sendDelfTimerName.c_str());
      delta_f->SendPhase([&]() {
        delta_f_gids->Send(pcms::Mode::Deferred);
        delta_f_gids2->Send(pcms::Mode::Deferred);
      });
      }

      const auto round_finish{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> round_elapsed_seconds{round_finish - round_start};
      if(!rank) std::cerr << "round " << i << " done in " << round_elapsed_seconds.count() << " seconds\n";
    }
  } while (!done);
  MPI_Barrier(comm);
  const auto finish{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{finish - start};
  if(!rank) std::cerr << "xgc_coupler " << elapsed_seconds.count() << "\n";
  }
  fclose(fhandle);
  Omega_h::vtk::write_parallel("coupler.vtk", &mesh, mesh.dim());
}

int main(int argc, char** argv)
{
  int provide;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provide);
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(!rank) std::cerr << "mpi thread level: " << provide << "\n";

  if (argc != 7) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                   "/path/to/partitionFile.cpn "
                   "sstDataTransport=[RDMA|WAN|MPI] "
                   "adiosEngine=[BP4|SST] "
                   "overlap=[small|large]\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 7);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  sstDataTransport = argv[4];
  adiosEngine = argv[5];
  overlapSize = argv[6];
  assert(overlapSize == "small" || overlapSize == "large");
  if(!rank) {
    std::cerr << "inputs: " << clientId << ", " << meshFile << " " << classPartitionFile << " "
              << sstDataTransport << " " << adiosEngine << " " << overlapSize << "\n";
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
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
