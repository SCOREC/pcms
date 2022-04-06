#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <redev_comm.h>
#include "wdmcpl.h"
#include "test_support.h"

namespace ts = test_support;

void prepareMsg(Omega_h::Mesh& mesh, redev::ClassPtn& partition,
    ts::OutMsg& out, redev::LOs& permute) {
  //transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classIds_h = Omega_h::HostRead(classIds);
  //count number of vertices going to each destination process by calling getRank - degree array
  std::map<int,int> destRankCounts;
  for(auto rank : partition.GetRanks() ) {
    destRankCounts[rank] = 0;
  }
  for(auto i=0; i<classIds_h.size(); i++) {
    auto dr = partition.GetRank(classIds_h[i]);
    assert(destRankCounts.count(dr));
    destRankCounts[dr]++;
  }
  REDEV_ALWAYS_ASSERT(destRankCounts[0] == 6);
  REDEV_ALWAYS_ASSERT(destRankCounts[1] == 13);
  //create dest and offsets arrays from degree array
  out.offset.resize(destRankCounts.size()+1);
  out.dest.resize(destRankCounts.size());
  out.offset[0] = 0;
  int i = 1;
  for(auto rankCount : destRankCounts) {
    out.dest[i-1] = rankCount.first;
    out.offset[i] = out.offset[i-1]+rankCount.second;
    i++;
  }
  redev::LOs expectedDest = {0,1};
  REDEV_ALWAYS_ASSERT(out.dest == expectedDest);
  redev::LOs expectedOffset = {0,6,19};
  REDEV_ALWAYS_ASSERT(out.offset == expectedOffset);
  //fill permutation array such that for vertex i permute[i] contains the
  //  position of vertex i's data in the message array
  std::map<int,int> destRankIdx;
  for(size_t i=0; i<out.dest.size(); i++) {
    auto dr = out.dest[i];
    destRankIdx[dr] = out.offset[i];
  }
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  permute.resize(classIds_h.size());
  for(auto i=0; i<classIds_h.size(); i++) {
    auto dr = partition.GetRank(classIds_h[i]);
    auto idx = destRankIdx[dr]++;
    permute[i] = idx;
  }
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  int rank = world->rank();
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh\n";
    std::cerr << "WARNING: this test is currently hardcoded for the xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh mesh\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  if(!rank) REDEV_ALWAYS_ASSERT(mesh.nelems() == 23); //sanity check that the loaded mesh is the expected one
  //partition the omegah mesh by classification and return the rank-to-classid array
  const auto classPartition = isRdv ? ts::migrateAndGetPartition(mesh) :
                                      ts::ClassificationPartition();
  if(isRdv) {
    ts::writeVtk(mesh,"rdvSplit",0);
  }
  auto partition = redev::ClassPtn(classPartition.ranks,classPartition.classIds);
  redev::Redev rdv(MPI_COMM_WORLD,partition,isRdv);
  rdv.Setup();
  const std::string name = "meshVtxIds";
  const int rdvRanks = 2;
  redev::AdiosComm<redev::GO> comm(MPI_COMM_WORLD, rdvRanks, rdv.getToEngine(), rdv.getToIO(), name);

  //build dest, offsets, and permutation arrays
  ts::OutMsg appOut = !isRdv ? ts::prepareAppOutMessage(mesh, partition) : ts::OutMsg();
  if(!isRdv) {
    redev::LOs expectedDest = {0,1};
    REDEV_ALWAYS_ASSERT(appOut.dest == expectedDest);
    redev::LOs expectedOffset = {0,6,19};
    REDEV_ALWAYS_ASSERT(appOut.offset == expectedOffset);
    redev::LOs expectedPermute = {0,6,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18};
    REDEV_ALWAYS_ASSERT(appOut.permute == expectedPermute);

    comm.SetOutMessageLayout(appOut.dest, appOut.offset);
  }

  redev::GOs rdvInPermute;

  for(int iter=0; iter<3; iter++) {
    if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //////////////////////////////////////////////////////
    //the non-rendezvous app sends global vtx ids to rendezvous
    //////////////////////////////////////////////////////
    if(!isRdv) {
      //fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(gids_h.size(),0);
      for(size_t i=0; i<msgs.size(); i++) {
        msgs[appOut.permute[i]] = gids_h[i];
      }
      auto start = std::chrono::steady_clock::now();
      comm.Send(msgs.data());
      ts::getAndPrintTime(start,name + " write",rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = comm.Unpack();
      const auto rdvIn = comm.GetInMessageLayout();
      REDEV_ALWAYS_ASSERT(rdvIn.offset == redev::GOs({0,6,19}));
      REDEV_ALWAYS_ASSERT(rdvIn.srcRanks == redev::GOs({0,0}));
      if(!rank) {
        REDEV_ALWAYS_ASSERT(rdvIn.start==0 && rdvIn.count==6);
        REDEV_ALWAYS_ASSERT(msgs == redev::GOs({0,2,3,4,5,6}));
      } else {
        REDEV_ALWAYS_ASSERT(rdvIn.start==6 && rdvIn.count==13);
        REDEV_ALWAYS_ASSERT(msgs == redev::GOs({1,7,8,9,10,11,12,13,14,15,16,17,18}));
      }
      ts::getAndPrintTime(start,name + " read",rank);
      //attach the ids to the mesh
      if(iter==0) ts::getRdvPermutation(mesh, msgs, rdvInPermute);
      ts::checkAndAttachIds(mesh, "inVtxGids", msgs, rdvInPermute);
      ts::writeVtk(mesh,"rdvInGids",iter);
    } //end non-rdv -> rdv
  } //end iter loop
  return 0;
}
