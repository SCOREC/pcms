#include "test_support.h"
#include <Omega_h_file.hpp> //vtk::write_parallel
#include <algorithm> // std::sort, std::stable_sort
#include <fstream> // ifstream
#include <mpi.h>

namespace test_support {

void writeVtk(Omega_h::Mesh& mesh, std::string_view name, int step) {
  std::stringstream ss;
  ss << name << step << ".vtk";
  Omega_h::vtk::write_parallel(ss.str(), &mesh, mesh.dim());
}

void timeMinMaxAvg(double time, double& min, double& max, double& avg) {
  const auto comm = MPI_COMM_WORLD;
  int nproc;
  MPI_Comm_size(comm, &nproc);
  double tot = 0;
  MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&time, &tot, 1, MPI_DOUBLE, MPI_SUM, comm);
  avg = tot / nproc;
}

void printTime(std::string_view mode, double min, double max, double avg) {
  std::stringstream ss;
  ss << mode << " elapsed time min, max, avg (s): "
     << min << " " << max << " " << avg << "\n";
  std::cout << ss.str();
}

void readClassPtnFromCpn(std::string cpnFileName, redev::LOs& ranks, redev::LOs& classIds) {
  std::ifstream in_file(cpnFileName);
  if(!in_file) {
    std::cerr << "cannot read " << cpnFileName << "... exiting\n";
    exit(EXIT_FAILURE);
  }

  int numLines;
  in_file >> numLines;
  ranks.reserve(numLines);
  classIds.reserve(numLines);
  int rank, classId;
  while(in_file >> classId >> rank) {
    ranks.push_back(rank);
    ranks.push_back(classId);
  }
  in_file.close();
}

void migrateMeshElms(Omega_h::Mesh& mesh, const redev::LOs& ranks, const redev::LOs& classIds) {
  auto ohComm = mesh.comm();
  auto mpiComm = ohcomm.get_impl();
  const auto rank = ohComm->rank();
  if(rank) REDEV_ALWAYS_ASSERT(mesh.nelems()==0); //only rank zero should have elements
  //loop over elements
  if(!rank) {
    const auto dim = mesh.dim();
    auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
    auto class_ids_h = Omega_h::HostRead(class_ids);
    typedef std::map<int,int> mii;
    mii classIdToRank;
    for(int i=0; i<ranks.size(); i++)
      classIdToRank[classIds[i]] = ranks[i];
    typedef std::map<int,vector<int>> miv;
    miv elmsPerRank;
    for(int i=0; i<mesh.nelems(); i++) {
      const auto dest = classIdToRank[class_ids_h[i]];
      elmsPerRank[dest].push_back(i);
    }
    for(auto iter = elmsPerRank.begin(); iter != elmsPerRank.end(); iter++) {
      //send the vector to the dest rank
      auto dest = iter.first;
      auto elms = iter.second;
      auto numElms = elms.size();
      MPI_Send(&numElms, 1, MPI_INT, dest, 0, mpiComm);
      MPI_Send(elms.data(), elms.size(), MPI_INT, dest, 0, mpiComm);
    }
  } else {
     const int src=0;
     int numElms;
     MPI_Recv(&numElms,1,MPI_INT,src,0,mpiComm,NULL);
     redev::LOs elms(numElms);
     MPI_Recv(elms.data(),numElms,MPI_INT,src,0,mpiComm,NULL);
  }
  //HERE - not tested
}

ClassificationPartition migrateAndGetPartition(Omega_h::Mesh& mesh) {
  auto ohComm = mesh.comm();
  const auto dim = mesh.dim();
  auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  Omega_h::ClassId max_class = Omega_h::get_max(class_ids);
  auto max_class_g = ohComm->allreduce(max_class,Omega_h_Op::OMEGA_H_MAX);
  REDEV_ALWAYS_ASSERT(ohComm->size() == max_class_g);
  auto class_ids_h = Omega_h::HostRead(class_ids);
  if(!ohComm->rank()) {
    //send ents with classId=2 to rank 1
    Omega_h::Write<Omega_h::I32> partitionRanks(5,0);
    Omega_h::Write<Omega_h::LO> partitionIdxs(5);
    Omega_h::fill_linear(partitionIdxs,0,1);
    auto owners = Omega_h::Remotes(partitionRanks, partitionIdxs);
    mesh.migrate(owners);
  } else {
    int firstElm = 5;
    const int elms = 18; // migrating elements [5:22]
    Omega_h::Write<Omega_h::I32> partitionRanks(elms,0);
    Omega_h::Write<Omega_h::LO> partitionIdxs(elms);
    Omega_h::fill_linear(partitionIdxs,firstElm,1);
    auto owners = Omega_h::Remotes(partitionRanks, partitionIdxs);
    mesh.migrate(owners);
  }

  //the hardcoded assignment of classids to ranks
  ClassificationPartition cp;
  cp.ranks.resize(3);
  cp.classIds.resize(3);
  cp.classIds[0] = 1; cp.ranks[0] = 0;
  cp.classIds[1] = 2; cp.ranks[1] = 1;
  cp.classIds[2] = 3; cp.ranks[2] = 0;  //center ('O point') model vertex
  return cp;
}

void checkAndAttachIds(Omega_h::Mesh& mesh, std::string_view name, const redev::GOs& vtxData, const redev::GOs& rdvPermute) {
  REDEV_ALWAYS_ASSERT(rdvPermute.size() == vtxData.size());
  auto gids_h = Omega_h::HostRead(mesh.globals(0));
  Omega_h::HostWrite<Omega_h::GO> inVtxData_h(mesh.nverts());
  for(int i=0; i<mesh.nverts(); i++)
    inVtxData_h[i] = -1;
  for(size_t i=0; i<rdvPermute.size(); i++) {
    inVtxData_h[rdvPermute[i]] = vtxData[i];
    REDEV_ALWAYS_ASSERT(gids_h[rdvPermute[i]] == vtxData[i]);
  }
  Omega_h::Write inVtxData(inVtxData_h);
  auto ohStr = std::string(name); //omegah wants a std::string const reference
  mesh.add_tag(0,ohStr,1,Omega_h::read(inVtxData));
  mesh.sync_tag(0,ohStr);
}

OutMsg prepareAppOutMessage(Omega_h::Mesh& mesh, const redev::ClassPtn& partition) {
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

  OutMsg out;
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

  //fill permutation array such that for vertex i permute[i] contains the
  //  position of vertex i's data in the message array
  std::map<int,int> destRankIdx;
  for(size_t i=0; i<out.dest.size(); i++) {
    auto dr = out.dest[i];
    destRankIdx[dr] = out.offset[i];
  }
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  out.permute.resize(classIds_h.size());
  for(auto i=0; i<classIds_h.size(); i++) {
    auto dr = partition.GetRank(classIds_h[i]);
    auto idx = destRankIdx[dr]++;
    out.permute[i] = idx;
  }
  return out;
}

redev::GOs getRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids) {
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  auto iGids = sortIndexes(gids_h);
  auto iInGids = sortIndexes(inGids);
  redev::GOs rdvPermute(inGids.size());
  int j=0;
  for(size_t i=0; i<inGids.size(); i++) {
    while(gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j!=gids_h.size()); //not found
    rdvPermute[iInGids[i]] = iGids[j];
  }
  return rdvPermute;
}

}
