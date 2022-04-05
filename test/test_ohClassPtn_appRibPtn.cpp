#include <numeric> // std::exclusive_scan
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_scalar.hpp> // divide_no_remainder
#include <redev_comm.h>
#include "wdmcpl.h"
#include "test_support.h"

namespace ts = test_support;

struct CSR {
  redev::GOs off;
  redev::GOs val;
};

//creates the outbound (rdv->non-rdv) permutation CSR given inGids and the rdv mesh instance
//this only needs to be computed once for each topological dimension
void getOutboundRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids, CSR& perm) {
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  auto iGids = ts::sort_indexes(gids_h);
  auto iInGids = ts::sort_indexes(inGids);
  //count the number of times each gid is included in inGids
  perm.off.resize(gids_h.size()+1);
  int j=0;
  for(size_t i=0; i<inGids.size(); i++) {
    while(gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j!=gids_h.size()); //found
    perm.off[iGids[j]]++;
  }
  //create the offsets array from the counts
  std::exclusive_scan(perm.off.begin(), perm.off.end(), perm.off.begin(), 0);
  //fill the permutation array
  perm.val.resize(perm.off.back());
  redev::LOs count(gids_h.size()); //how many times each gid was written
  j=0;
  for(size_t i=0; i<inGids.size(); i++) {
    while(gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j!=gids_h.size()); //found
    const auto subIdx = count[iGids[j]]++;
    const auto startIdx = perm.off[iGids[j]];
    const auto offIdx = startIdx + subIdx;
    perm.val[offIdx] = iInGids[i];
  }
}

void prepareRdvOutMessage(Omega_h::Mesh& mesh, ts::InMsg const& in, ts::OutMsg& out, CSR& permute){
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  const auto nproc = ohComm->size();
  auto nAppProcs = Omega_h::divide_no_remainder(in.srcRanks.size(),static_cast<size_t>(nproc));
  REDEV_ALWAYS_ASSERT(nAppProcs==2);
  //build dest and offsets arrays from incoming message metadata
  redev::LOs senderDeg(nAppProcs);
  for(size_t i=0; i<nAppProcs-1; i++) {
    senderDeg[i] = in.srcRanks[(i+1)*nproc+rank] - in.srcRanks[i*nproc+rank];
  }
  const auto totInMsgs = in.offset[rank+1]-in.offset[rank];
  senderDeg[nAppProcs-1] = totInMsgs - in.srcRanks[(nAppProcs-1)*nproc+rank];
  if(!rank) REDEV_ALWAYS_ASSERT( senderDeg == redev::LOs({4,5}) );
  if(rank) REDEV_ALWAYS_ASSERT( senderDeg == redev::LOs({8,7}) );
  for(size_t i=0; i<nAppProcs; i++) {
    if(senderDeg[i] > 0) {
      out.dest.push_back(i);
    }
  }
  REDEV_ALWAYS_ASSERT( out.dest == redev::LOs({0,1}) );
  redev::GO sum = 0;
  for(auto deg : senderDeg) { //exscan over values > 0
    if(deg>0) {
      out.offset.push_back(sum);
      sum+=deg;
    }
  }
  out.offset.push_back(sum);
  if(!rank) REDEV_ALWAYS_ASSERT( out.offset == redev::LOs({0,4,9}) );
  if(rank) REDEV_ALWAYS_ASSERT( out.offset == redev::LOs({0,8,15}) );
  getOutboundRdvPermutation(mesh, in.msgs, permute);
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh\n";
    std::cerr << "WARNING: this test is currently hardcoded for the xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh\n";
    std::cerr << "mesh for the rendezvous processes and xgc1_data/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh mesh/2p.osh for the\n";
    std::cerr << "for the non-rendezvous processes\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  redev::LOs ranks;
  redev::LOs classIds;
  if(isRdv) {
    //partition the omegah mesh by classification and return the
    //rank-to-classid array
    ts::getClassPartition(mesh, ranks, classIds);
    REDEV_ALWAYS_ASSERT(ranks.size()==3);
    REDEV_ALWAYS_ASSERT(ranks.size()==classIds.size());
    ts::writeVtk(mesh,"rdvSplit",0);
  } else {
    REDEV_ALWAYS_ASSERT(world->size()==2);
    if(!rank) REDEV_ALWAYS_ASSERT(mesh.nelems()==11);
    ts::writeVtk(mesh,"appSplit",0);
  }
  auto partition = redev::ClassPtn(ranks,classIds);
  redev::Redev rdv(MPI_COMM_WORLD,partition,isRdv);
  rdv.Setup();

  const std::string name = "meshVtxIds";
  const int rdvRanks = 2;
  const int appRanks = 2;
  redev::AdiosComm<redev::GO> commA2R(MPI_COMM_WORLD, rdvRanks, rdv.getToEngine(), rdv.getToIO(), name+"_A2R");
  redev::AdiosComm<redev::GO> commR2A(MPI_COMM_WORLD, appRanks, rdv.getFromEngine(), rdv.getFromIO(), name+"_R2A");

  redev::LOs appOutPermute;
  ts::OutMsg appOut;
  ts::InMsg appIn;

  redev::GOs rdvInPermute;
  CSR rdvOutPermute;
  ts::OutMsg rdvOut;
  ts::InMsg rdvIn;

  for(int iter=0; iter<3; iter++) {
    if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //////////////////////////////////////////////////////
    //the non-rendezvous app sends global vtx ids to rendezvous
    //////////////////////////////////////////////////////
    if(!isRdv) {
      //build dest, offsets, and permutation arrays
      if(iter==0) ts::prepareAppOutMessage(mesh, partition, appOut, appOutPermute);
      //fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(gids_h.size(),0);
      for(size_t i=0; i<msgs.size(); i++) {
        msgs[appOutPermute[i]] = gids_h[i];
      }
      auto start = std::chrono::steady_clock::now();
      commA2R.Pack(appOut.dest, appOut.offset, msgs.data());
      commA2R.Send();
      ts::getAndPrintTime(start,name + " appWrite",rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const bool knownSizes = (iter == 0) ? false : true;
      ts::unpack(commA2R,knownSizes,rdvIn);
      ts::getAndPrintTime(start,name + " rdvRead",rank);
      //attach the ids to the mesh
      if(iter==0) ts::getRdvPermutation(mesh, rdvIn.msgs, rdvInPermute);
      ts::checkAndAttachIds(mesh, "inVtxGids", rdvIn.msgs, rdvInPermute);
      ts::writeVtk(mesh,"rdvInGids",iter);
    } //end non-rdv -> rdv
    //////////////////////////////////////////////////////
    //the rendezvous app sends global vtx ids to non-rendezvous
    //////////////////////////////////////////////////////
    if(isRdv) {
      //build dest, offsets, and permutation arrays
      if( iter==0 ) prepareRdvOutMessage(mesh,rdvIn,rdvOut,rdvOutPermute);
      //fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(rdvOutPermute.off.back());
      for(int i=0; i<gids_h.size(); i++) {
        for(int j=rdvOutPermute.off[i]; j<rdvOutPermute.off[i+1]; j++) {
          msgs[rdvOutPermute.val[j]] = gids_h[i];
        }
      }
      auto start = std::chrono::steady_clock::now();
      commR2A.Pack(rdvOut.dest, rdvOut.offset, msgs.data());
      commR2A.Send();
      ts::getAndPrintTime(start,name + " rdvWrite",rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const bool knownSizes = (iter == 0) ? false : true;
      ts::unpack(commR2A,knownSizes,appIn);
      ts::getAndPrintTime(start,name + " appRead",rank);
      { //check incoming messages are in the correct order
        auto gids = mesh.globals(0);
        auto gids_h = Omega_h::HostRead(gids);
        REDEV_ALWAYS_ASSERT(appIn.count == static_cast<size_t>(gids_h.size()));
        for(size_t i=0; i<appIn.msgs.size(); i++) {
          REDEV_ALWAYS_ASSERT(gids_h[i] == appIn.msgs[appOutPermute[i]]);
        }
      }
    } //end rdv -> non-rdv
  } //end iter loop
  return 0;
}
