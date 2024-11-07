#include <numeric> // std::exclusive_scan
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_scalar.hpp> // divide_no_remainder
#include "pcms.h"
#include "test_support.h"

namespace ts = test_support;

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh\n";
    std::cerr << "WARNING: this test is currently hardcoded for the "
                 "xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh\n";
    std::cerr << "mesh for the rendezvous processes and "
                 "xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh mesh/2p.osh "
                 "for the\n";
    std::cerr << "for the non-rendezvous processes\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  // partition the omegah mesh by classification and return the rank-to-classid
  // array
  const auto classPartition =
    isRdv ? ts::migrateAndGetPartition(mesh) : ts::ClassificationPartition();
  if (isRdv) {
    ts::writeVtk(mesh, "rdvSplit", 0);
  } else {
    REDEV_ALWAYS_ASSERT(world->size() == 2);
    if (!rank)
      REDEV_ALWAYS_ASSERT(mesh.nelems() == 11);
    ts::writeVtk(mesh, "appSplit", 0);
  }
  auto partition = redev::ClassPtn(MPI_COMM_WORLD, classPartition.ranks,
                                   classPartition.modelEnts);
  redev::Redev rdv(MPI_COMM_WORLD, redev::Partition{std::move(partition)},
                   static_cast<redev::ProcessType>(isRdv));

  const std::string name = "meshVtxIds";
  const int rdvRanks = 2;
  const int appRanks = 2;

  adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "12"}};
  auto channel =
    rdv.CreateAdiosChannel(name, params, redev::TransportType::BP4);
  auto commPair = channel.CreateComm<redev::GO>(name,rdv.GetMPIComm());

  // Build the dest, offsets, and permutation arrays for the forward
  // send from non-rendezvous to rendezvous.
  ts::OutMsg appOut =
    !isRdv ? ts::prepareAppOutMessage(
               mesh, std::get<decltype(partition)>(rdv.GetPartition()))
           : ts::OutMsg();
  if (!isRdv) {
    commPair.SetOutMessageLayout(appOut.dest, appOut.offset);
  }

  redev::GOs rdvInPermute;
  ts::CSR rdvOutPermute;
  ts::OutMsg rdvOut;

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
      ts::getAndPrintTime(start, name + " appWrite", rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = channel.ReceivePhase(
        [&]() { return commPair.Recv(redev::Mode::Synchronous); });
      ts::getAndPrintTime(start, name + " rdvRead", rank);
      // attach the ids to the mesh
      if (iter == 0) {
        // We have received the first input message in the rendezvous
        // processes.  Using the meta data of the incoming message we will:
        //- compute the permutation from the incoming vertex global ids to the
        //   on-process global ids
        //- set the message layout for the reverse (rendezvous->non-rendezvous)
        // send by
        //   building the dest and offsets array.
        //- compute the reverse send's permutation array using the layout of
        //   global vertex ids in 'msgs'.
        // These operations only need to be done once per coupling as long as
        // the topology and partition of the rendezvous and non-rendezvous
        // meshes remains the same.
        auto rdvIn = commPair.GetInMessageLayout();
        rdvInPermute = ts::getRdvPermutation(mesh, msgs);
        rdvOut = ts::prepareRdvOutMessage(mesh, rdvIn);
        REDEV_ALWAYS_ASSERT(rdvOut.dest == redev::LOs({0, 1}));
        if (!rank)
          REDEV_ALWAYS_ASSERT(rdvOut.offset == redev::LOs({0, 4, 9}));
        if (rank)
          REDEV_ALWAYS_ASSERT(rdvOut.offset == redev::LOs({0, 8, 15}));
        commPair.SetOutMessageLayout(rdvOut.dest, rdvOut.offset);
        rdvOutPermute = ts::getRdvOutPermutation(mesh, msgs);
      }
      ts::checkAndAttachIds(mesh, "inVtxGids", msgs, rdvInPermute);
      ts::writeVtk(mesh, "rdvInGids", iter);
    } // end non-rdv -> rdv
    //////////////////////////////////////////////////////
    // the rendezvous app sends global vtx ids to non-rendezvous
    //////////////////////////////////////////////////////
    if (isRdv) {
      // fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(rdvOutPermute.off.back());
      for (int i = 0; i < gids_h.size(); i++) {
        for (int j = rdvOutPermute.off[i]; j < rdvOutPermute.off[i + 1]; j++) {
          msgs[rdvOutPermute.val[j]] = gids_h[i];
        }
      }
      auto start = std::chrono::steady_clock::now();
      channel.SendPhase([&]() { commPair.Send(msgs.data()); });
      ts::getAndPrintTime(start, name + " rdvWrite", rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = channel.ReceivePhase(
        [&]() { return commPair.Recv(redev::Mode::Synchronous); });
      ts::getAndPrintTime(start, name + " appRead", rank);
      { // check incoming messages are in the correct order
        auto gids = mesh.globals(0);
        auto gids_h = Omega_h::HostRead(gids);
        REDEV_ALWAYS_ASSERT(msgs.size() == static_cast<size_t>(gids_h.size()));
        for (size_t i = 0; i < msgs.size(); i++) {
          REDEV_ALWAYS_ASSERT(gids_h[i] == msgs[appOut.permute[i]]);
        }
      }
    } // end rdv -> non-rdv
  }   // end iter loop
  return 0;
}
