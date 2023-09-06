#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include "pcms.h"
#include "test_support.h"

namespace ts = test_support;

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  int rank = world->rank();
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh\n";
    std::cerr << "WARNING: this test is currently hardcoded for the "
                 "xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh mesh\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  if (!rank)
    REDEV_ALWAYS_ASSERT(
      mesh.nelems() ==
      23); // sanity check that the loaded mesh is the expected one
  // partition the omegah mesh by classification and return the rank-to-classid
  // array
  const auto classPartition =
    isRdv ? ts::migrateAndGetPartition(mesh) : ts::ClassificationPartition();
  if (isRdv) {
    ts::writeVtk(mesh, "rdvSplit", 0);
  }
  auto partition = redev::ClassPtn(MPI_COMM_WORLD, classPartition.ranks,
                                   classPartition.modelEnts);
  redev::Redev rdv(MPI_COMM_WORLD, redev::Partition{std::move(partition)},
                   static_cast<redev::ProcessType>(isRdv));
  const std::string name = "meshVtxIds";
  adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "2"}};
  auto channel =
    rdv.CreateAdiosChannel(name, params, redev::TransportType::BP4);
  auto commPair = channel.CreateComm<redev::GO>(name,rdv.GetMPIComm());

  // build dest, offsets, and permutation arrays
  ts::OutMsg appOut =
    !isRdv ? ts::prepareAppOutMessage(
               mesh, std::get<decltype(partition)>(rdv.GetPartition()))
           : ts::OutMsg();
  if (!isRdv) {
    redev::LOs expectedDest = {0, 1};
    REDEV_ALWAYS_ASSERT(appOut.dest == expectedDest);
    redev::LOs expectedOffset = {0, 6, 19};
    REDEV_ALWAYS_ASSERT(appOut.offset == expectedOffset);
    redev::LOs expectedPermute = {0,  6,  1,  2,  3,  4,  5,  7,  8, 9,
                                  10, 11, 12, 13, 14, 15, 16, 17, 18};
    REDEV_ALWAYS_ASSERT(appOut.permute == expectedPermute);

    commPair.SetOutMessageLayout(appOut.dest, appOut.offset);
  }

  redev::GOs rdvInPermute;

  for (int iter = 0; iter < 3; iter++) {
    if (!rank)
      fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //////////////////////////////////////////////////////
    // the non-rendezvous app sends global vtx ids to rendezvous
    //////////////////////////////////////////////////////
    if (!isRdv) {
      // fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(gids_h.size(), 0);
      for (size_t i = 0; i < msgs.size(); i++) {
        msgs[appOut.permute[i]] = gids_h[i];
      }
      auto start = std::chrono::steady_clock::now();
      channel.SendPhase([&]() { commPair.Send(msgs.data()); });
      ts::getAndPrintTime(start, name + " write", rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = channel.ReceivePhase([&]() { return commPair.Recv(); });
      const auto rdvIn = commPair.GetInMessageLayout();
      REDEV_ALWAYS_ASSERT(rdvIn.offset == redev::GOs({0, 6, 19}));
      REDEV_ALWAYS_ASSERT(rdvIn.srcRanks == redev::GOs({0, 0}));
      if (!rank) {
        REDEV_ALWAYS_ASSERT(rdvIn.start == 0 && rdvIn.count == 6);
        REDEV_ALWAYS_ASSERT(msgs == redev::GOs({0, 2, 3, 4, 5, 6}));
      } else {
        REDEV_ALWAYS_ASSERT(rdvIn.start == 6 && rdvIn.count == 13);
        REDEV_ALWAYS_ASSERT(
          msgs == redev::GOs({1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}));
      }
      ts::getAndPrintTime(start, name + " read", rank);
      // attach the ids to the mesh
      if (iter == 0)
        rdvInPermute = ts::getRdvPermutation(mesh, msgs);
      ts::checkAndAttachIds(mesh, "inVtxGids", msgs, rdvInPermute);
      ts::writeVtk(mesh, "rdvInGids", iter);
    } // end non-rdv -> rdv
  }   // end iter loop
  return 0;
}
