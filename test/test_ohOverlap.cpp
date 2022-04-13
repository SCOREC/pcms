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

//creates the rdv->non-rdv permutation CSR given inGids and the rdv mesh instance
CSR getRdvOutPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids) {
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  auto iGids = ts::sortIndexes(gids_h);
  auto iInGids = ts::sortIndexes(inGids);
  //count the number of times each gid is included in inGids
  CSR perm;
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
  return perm;
}

ts::OutMsg prepareRdvOutMessage(Omega_h::Mesh& mesh, const redev::InMessageLayout& in) {
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
  ts::OutMsg out;
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
  return out;
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh /path/to/partitionFile.cpn\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  std::string_view cpnFileName(argv[3]);
  ts::ClassificationPartition classPartition;
  if(isRdv) {
    const auto facePartition = !rank ? ts::readClassPartitionFile(cpnFileName) : ts::ClassificationPartition();
    ts::migrateMeshElms(mesh, facePartition);
    classPartition = ts::CreateClassificationPartition(mesh);
    ts::writeVtk(mesh,"rdvClassPtn",0);
  } else {
    ts::writeVtk(mesh,"appPartition",0);
  }
  auto partition = redev::ClassPtn(classPartition.ranks,classPartition.modelEnts);
  partition.Gather(MPI_COMM_WORLD); //TODO should be included in rdv.Setup - faild here for 37k d3d - split ownership on part boundary
  redev::Redev rdv(MPI_COMM_WORLD,partition,isRdv);
  rdv.Setup();

  return 0;
}
