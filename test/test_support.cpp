#include "test_support.h"
#include <Omega_h_file.hpp>    //vtk::write_parallel
#include <Omega_h_atomics.hpp> //atomic_increment
#include <Omega_h_for.hpp>     //parallel_for
#include <algorithm>           // std::sort, std::stable_sort
#include <fstream>             // ifstream
#include <mpi.h>
#include <redev_exclusive_scan.h>

namespace test_support
{

void writeVtk(Omega_h::Mesh& mesh, std::string_view name, int step)
{
  std::stringstream ss;
  ss << name << step << ".vtk";
  Omega_h::vtk::write_parallel(ss.str(), &mesh, mesh.dim());
}

void timeMinMaxAvg(double time, double& min, double& max, double& avg)
{
  const auto comm = MPI_COMM_WORLD;
  int nproc;
  MPI_Comm_size(comm, &nproc);
  double tot = 0;
  MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&time, &tot, 1, MPI_DOUBLE, MPI_SUM, comm);
  avg = tot / nproc;
}

void printTime(std::string_view mode, double min, double max, double avg)
{
  std::stringstream ss;
  ss << mode << " elapsed time min, max, avg (s): " << min << " " << max << " "
     << avg << "\n";
  std::cout << ss.str();
}

ClassificationPartition readClassPartitionFile(std::string_view cpnFileName)
{
  std::string cpnFileNameStr(cpnFileName);
  std::ifstream in_file(cpnFileNameStr); // doesn't take a string_view
  if (!in_file) {
    std::stringstream ss;
    ss << "cannot read " << cpnFileName << "... exiting\n";
    std::cout << ss.str();
    exit(EXIT_FAILURE);
  }

  int numLines;
  in_file >> numLines;
  ClassificationPartition cp;
  cp.ranks.reserve(numLines);
  cp.modelEnts.reserve(numLines);
  int rank, classId;
  while (in_file >> classId >> rank) {
    cp.ranks.push_back(rank);
    const redev::ClassPtn::ModelEnt ent(
      2, classId); // cpn files only contain model faces!
    cp.modelEnts.push_back(ent);
  }
  in_file.close();
  return cp;
}

/**
 * return the permutation array of mesh entities of dimension 'dim' that orders
 * the list of model entities they are classified on (defined by a pair of
 * integers for the model entity id and dimension) in ascending order.
 */
Omega_h::LOs getModelEntityPermutation(Omega_h::Mesh& mesh, const int dim)
{
  auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  auto class_dims = mesh.get_array<Omega_h::I8>(dim, "class_dim");
  auto ids = Omega_h::HostRead(class_ids);
  auto dims = Omega_h::HostRead(class_dims);
  Omega_h::HostWrite<Omega_h::LO> idx(ids.size());
  std::iota(&idx[0], (&idx[ids.size() - 1]) + 1, 0);
  std::stable_sort(&idx[0], (&idx[idx.size() - 1]) + 1,
                   [&](const size_t lhs, const size_t rhs) {
                     const auto ldim = dims[lhs];
                     const auto lid = ids[lhs];
                     const auto rdim = dims[rhs];
                     const auto rid = ids[rhs];
                     const auto dimLess = (ldim < rdim);
                     const auto dimEqIdLess = (ldim == rdim) && (lid < rid);
                     const auto res = dimLess || dimEqIdLess;
                     return res;
                   });
  Omega_h::LOs perm(idx);
  return perm;
}

Omega_h::LO countModelEnts(Omega_h::Mesh& mesh, const Omega_h::LOs permutation,
                           const int dim)
{
  auto ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  auto dims = mesh.get_array<Omega_h::I8>(dim, "class_dim");
  REDEV_ALWAYS_ASSERT(ids.size() == permutation.size());
  Omega_h::Write<Omega_h::LO> numModelEnts(1, 0);
  auto isSameModelEnt = OMEGA_H_LAMBDA(int i, int j)
  {
    const auto ip = permutation[i];
    const auto jp = permutation[j];
    return (ids[ip] == ids[jp]) && (dims[ip] == dims[jp]);
  };
  // only looking for equal order classification
  //  (mesh ent dimension == model ent dimension)
  auto isCurrentModelEntDim = OMEGA_H_LAMBDA(int i)
  {
    const auto ip = permutation[i];
    return (dims[ip] == dim);
  };
  auto countEnts = OMEGA_H_LAMBDA(int i)
  {
    const auto currentDim = isCurrentModelEntDim(i);
    const auto prevSame = (i == 0) ? false : isSameModelEnt(i, i - 1);
    if (!currentDim || prevSame)
      return;
    Omega_h::atomic_increment(&numModelEnts[0]);
  };
  Omega_h::parallel_for(ids.size(), countEnts);
  Omega_h::HostRead numModelEnts_h(Omega_h::read(numModelEnts));
  return numModelEnts_h.last();
}

struct ModelEntityOwners
{
  Omega_h::LOs ids;
  Omega_h::Read<Omega_h::I8> dims;
  Omega_h::LOs owners;
};

/**
 * For each model entity find the minimum rank that has a mesh entity
 * classified on it.  This rank is defined as the owner of that model entity in
 * the classification partition.
 */
ModelEntityOwners getModelEntityOwners(Omega_h::Mesh& mesh,
                                       const Omega_h::LOs permutation,
                                       const int dim, const int numModelEnts)
{
  auto ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  auto dims = mesh.get_array<Omega_h::I8>(dim, "class_dim");
  auto remotes = mesh.ask_owners(dim);
  const auto numEnts = ids.size();

  auto mdlEntIds = Omega_h::Write<Omega_h::LO>(numModelEnts);
  auto mdlEntDims = Omega_h::Write<Omega_h::I8>(numModelEnts);
  auto mdlEntOwners = Omega_h::Write<Omega_h::LO>(numModelEnts);

  Omega_h::Write<Omega_h::LO> count(1, 0);

  REDEV_ALWAYS_ASSERT(ids.size() == permutation.size());
  auto isSameModelEnt = OMEGA_H_LAMBDA(int i, int j)
  {
    const auto ip = permutation[i];
    const auto jp = permutation[j];
    return (ids[ip] == ids[jp]) && (dims[ip] == dims[jp]);
  };
  // only looking for equal order classification
  //  (mesh ent dimension == model ent dimension)
  auto isCurrentModelEntDim = OMEGA_H_LAMBDA(int i)
  {
    const auto ip = permutation[i];
    return (dims[ip] == dim);
  };
  // Loop through the model entities and then find the minimum
  // rank of mesh entities classified on it (the inner while loop).
  // Note, this function only has N concurrent indices doing work
  // where N is the number of model entities of the active dimension.
  // It is essentially a min-reduction over a bunch of small arrays
  //(one array for each model entity).  Could be written as a segment-sort.
  auto getEnts = OMEGA_H_LAMBDA(int i)
  {
    const auto currentDim = isCurrentModelEntDim(i);
    const auto prevSame = (i == 0) ? false : isSameModelEnt(i, i - 1);
    if (!currentDim || prevSame)
      return;
    auto minOwner = remotes.ranks[permutation[i]];
    auto next = i + 1;
    while ((next < numEnts) && isSameModelEnt(i, next)) {
      auto nextOwner = remotes.ranks[permutation[next]];
      if (minOwner > nextOwner)
        minOwner = nextOwner;
      next++;
    }
    const int mdlEntIdx = Omega_h::atomic_fetch_add(&count[0], 1);
    mdlEntIds[mdlEntIdx] = ids[permutation[i]];
    mdlEntDims[mdlEntIdx] = dims[permutation[i]];
    mdlEntOwners[mdlEntIdx] = minOwner;
  };
  Omega_h::parallel_for(numEnts, getEnts);
  ModelEntityOwners meow{mdlEntIds, mdlEntDims, mdlEntOwners};
  return meow;
}

/**
 * I don't feel like writing the array merge... so we dump it into a map
 */
void append(const ModelEntityOwners& meow,
            redev::ClassPtn::ModelEntToRank& entToRank)
{
  auto ids = Omega_h::HostRead(meow.ids);
  auto dims = Omega_h::HostRead(meow.dims);
  auto owners = Omega_h::HostRead(meow.owners);
  for (size_t i = 0; i < ids.size(); i++) {
    const auto ent = redev::ClassPtn::ModelEnt(dims[i], ids[i]);
    if (entToRank.count(ent)) {
      REDEV_ALWAYS_ASSERT(entToRank[ent] == owners[i]);
    } else {
      entToRank[ent] = owners[i];
    }
  }
}

ClassificationPartition fromMap(redev::ClassPtn::ModelEntToRank& entToRank)
{
  ClassificationPartition cp;
  cp.ranks.reserve(entToRank.size());
  cp.modelEnts.reserve(entToRank.size());
  for (auto iter = entToRank.begin(); iter != entToRank.end(); iter++) {
    cp.modelEnts.push_back(iter->first);
    cp.ranks.push_back(iter->second);
  }
  return cp;
}

/**
 * given a mesh partitioned based on model face ownership, derive model vertex
 * and edge assignment to mpi ranks from the ownership of mesh entities on the
 * partition boundaries
 */
ClassificationPartition CreateClassificationPartition(Omega_h::Mesh& mesh)
{
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  redev::ClassPtn::ModelEntToRank m2r;
  for (int dim = 0; dim <= mesh.dim(); dim++) {
    auto perm = getModelEntityPermutation(mesh, dim);
    auto numModelEnts = countModelEnts(mesh, perm, dim);
    auto modelEntityOwners =
      getModelEntityOwners(mesh, perm, dim, numModelEnts);
    append(modelEntityOwners, m2r);
  }
  return fromMap(m2r);
}

void migrateMeshElms(Omega_h::Mesh& mesh,
                     const ClassificationPartition& partition)
{
  auto ohComm = mesh.comm();
  auto mpiComm = ohComm->get_impl();
  const auto rank = ohComm->rank();
  if (rank)
    REDEV_ALWAYS_ASSERT(mesh.nelems() ==
                        0); // only rank zero should have elements
  if (!rank) {
    const auto dim = mesh.dim();
    auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
    auto class_ids_h = Omega_h::HostRead(class_ids);
    auto class_dims = mesh.get_array<Omega_h::I8>(dim, "class_dim");
    auto class_dims_h = Omega_h::HostRead(class_dims);
    using ModelEnt = redev::ClassPtn::ModelEnt;
    std::map<ModelEnt, int> modelEntToRank;
    for (int i = 0; i < partition.ranks.size(); i++)
      modelEntToRank[partition.modelEnts[i]] = partition.ranks[i];
    typedef std::map<int, std::vector<int>> miv;
    miv elemsPerRank;
    for (int i = 0; i < mesh.nelems(); i++) {
      const ModelEnt ent({class_dims_h[i], class_ids_h[i]});
      REDEV_ALWAYS_ASSERT(modelEntToRank.count(ent));
      const auto dest = modelEntToRank[ent];
      elemsPerRank[dest].push_back(i);
    }
    // make sure we are not sending elements to ranks that don't exist
    REDEV_ALWAYS_ASSERT(elemsPerRank.size() == ohComm->size());
    for (auto iter = elemsPerRank.begin(); iter != elemsPerRank.end(); iter++) {
      const auto dest = iter->first;
      REDEV_ALWAYS_ASSERT(dest <
                          ohComm->size()); // the destination rank must exist
      if (dest) {
        const auto elms = iter->second;
        const auto numElms = elms.size();
        // send the array of element ids each destination rank will be the owner
        // of
        MPI_Send(&numElms, 1, MPI_INT, dest, 0, mpiComm);
        MPI_Send(elms.data(), elms.size(), MPI_INT, dest, 0, mpiComm);
      }
    }
    // fill a device array with the element ids that remain on rank 0
    const auto elems = elemsPerRank[0];
    const auto numElems = elems.size();
    Omega_h::HostWrite<Omega_h::LO> elemIdxs(numElems);
    for (size_t i = 0; i < numElems; i++)
      elemIdxs[i] = elems[i];
    Omega_h::HostWrite<Omega_h::I32> elemRanks(numElems);
    std::fill(elemRanks.data(), elemRanks.data() + numElems, 0);
    // copy to device
    Omega_h::Write elemIdxs_d(elemIdxs);
    Omega_h::Write elemRanks_d(elemRanks);
    // create remotes
    auto owners = Omega_h::Remotes(elemRanks_d, elemIdxs_d);
    // call migrate
    mesh.migrate(owners);
  } else {
    const int src = 0;
    int numElems;
    MPI_Status stat;
    MPI_Recv(&numElems, 1, MPI_INT, src, 0, mpiComm, &stat);
    Omega_h::HostWrite<Omega_h::LO> elemIdxs(numElems);
    MPI_Recv(elemIdxs.data(), numElems, MPI_INT, src, 0, mpiComm, &stat);
    Omega_h::HostWrite<Omega_h::I32> elemRanks(numElems);
    std::fill(elemRanks.data(), elemRanks.data() + numElems, 0);
    // copy to device
    Omega_h::Write elemIdxs_d(elemIdxs);
    Omega_h::Write elemRanks_d(elemRanks);
    // create remotes
    auto owners = Omega_h::Remotes(elemRanks_d, elemIdxs_d);
    // call migrate
    mesh.migrate(owners);
  }
}

ClassificationPartition migrateAndGetPartition(Omega_h::Mesh& mesh)
{
  auto ohComm = mesh.comm();
  const auto dim = mesh.dim();
  auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  Omega_h::ClassId max_class = Omega_h::get_max(class_ids);
  auto max_class_g = ohComm->allreduce(max_class, Omega_h_Op::OMEGA_H_MAX);
  REDEV_ALWAYS_ASSERT(ohComm->size() == max_class_g);
  auto class_ids_h = Omega_h::HostRead(class_ids);
  if (!ohComm->rank()) {
    // send ents with classId=2 to rank 1
    Omega_h::Write<Omega_h::I32> partitionRanks(5, 0);
    Omega_h::Write<Omega_h::LO> partitionIdxs(5);
    Omega_h::fill_linear(partitionIdxs, 0, 1);
    auto owners = Omega_h::Remotes(partitionRanks, partitionIdxs);
    mesh.migrate(owners);
  } else {
    int firstElm = 5;
    const int elms = 18; // migrating elements [5:22]
    Omega_h::Write<Omega_h::I32> partitionRanks(elms, 0);
    Omega_h::Write<Omega_h::LO> partitionIdxs(elms);
    Omega_h::fill_linear(partitionIdxs, firstElm, 1);
    auto owners = Omega_h::Remotes(partitionRanks, partitionIdxs);
    mesh.migrate(owners);
  }
  auto cp = CreateClassificationPartition(mesh);
  // check the hardcoded assignment of classids to ranks
  if (!ohComm->rank()) {
    redev::ClassPtn::ModelEntVec expectedEnts{{0, 1}, {0, 3}, {1, 1}, {2, 1}};
    redev::LOs expectedRanks{{0, 0, 0, 0}};
    REDEV_ALWAYS_ASSERT(cp.ranks == expectedRanks);
    REDEV_ALWAYS_ASSERT(cp.modelEnts == expectedEnts);
  } else {
    redev::ClassPtn::ModelEntVec expectedEnts{
      {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
    redev::LOs expectedRanks{{0, 1, 0, 1, 1}};
    REDEV_ALWAYS_ASSERT(cp.ranks == expectedRanks);
    REDEV_ALWAYS_ASSERT(cp.modelEnts == expectedEnts);
  }
  return cp;
}

void checkAndAttachIds(Omega_h::Mesh& mesh, std::string_view name,
                       const redev::GOs& vtxData, const redev::GOs& rdvPermute)
{
  REDEV_ALWAYS_ASSERT(rdvPermute.size() == vtxData.size());
  auto gids_h = Omega_h::HostRead(mesh.globals(0));
  Omega_h::HostWrite<Omega_h::GO> inVtxData_h(mesh.nverts());
  for (int i = 0; i < mesh.nverts(); i++)
    inVtxData_h[i] = -1;
  for (size_t i = 0; i < rdvPermute.size(); i++) {
    inVtxData_h[rdvPermute[i]] = vtxData[i];
    REDEV_ALWAYS_ASSERT(gids_h[rdvPermute[i]] == vtxData[i]);
  }
  Omega_h::Write inVtxData(inVtxData_h);
  auto ohStr = std::string(name); // omegah wants a std::string const reference
  mesh.add_tag(0, ohStr, 1, Omega_h::read(inVtxData));
  mesh.sync_tag(0, ohStr);
}

OutMsg prepareAppOutMessage(Omega_h::Mesh& mesh,
                            const redev::ClassPtn& partition)
{
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classIds_h = Omega_h::HostRead(classIds);
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto classDims_h = Omega_h::HostRead(classDims);
  auto isOverlap =
    mesh.has_tag(0, "isOverlap")
      ? mesh.get_array<Omega_h::I8>(0, "isOverlap")
      : Omega_h::Read<Omega_h::I8>(classIds.size(), 1,
                                   "isOverlap"); // no mask for overlap vertices
  auto isOverlap_h = Omega_h::HostRead(isOverlap);
  // count number of vertices going to each destination process by calling
  // getRank - degree array
  std::map<int, int> destRankCounts;
  for (auto rank : partition.GetRanks()) {
    destRankCounts[rank] = 0;
  }
  for (auto i = 0; i < classIds_h.size(); i++) {
    if (isOverlap_h[i]) {
      const auto ent =
        redev::ClassPtn::ModelEnt({classDims_h[i], classIds_h[i]});
      auto dr = partition.GetRank(ent);
      assert(destRankCounts.count(dr));
      destRankCounts[dr]++;
    }
  }

  OutMsg out;
  // create dest and offsets arrays from degree array
  out.offset.resize(destRankCounts.size() + 1);
  out.dest.resize(destRankCounts.size());
  out.offset[0] = 0;
  // isn't this just an exclusive scan for the offset?
  int i = 1;
  for (auto rankCount : destRankCounts) {
    out.dest[i - 1] = rankCount.first;
    out.offset[i] = out.offset[i - 1] + rankCount.second;
    i++;
  }

  // fill permutation array such that for vertex i permute[i] contains the
  //   position of vertex i's data in the message array
  std::map<int, int> destRankIdx;
  for (size_t i = 0; i < out.dest.size(); i++) {
    auto dr = out.dest[i];
    destRankIdx[dr] = out.offset[i];
  }
  // auto gids = mesh.globals(0);
  // auto gids_h = Omega_h::HostRead(gids);
  out.permute.resize(out.offset.back());
  int j = 0;
  for (auto i = 0; i < classIds_h.size(); i++) {
    if (isOverlap_h[i]) {
      const auto ent =
        redev::ClassPtn::ModelEnt({classDims_h[i], classIds_h[i]});
      auto dr = partition.GetRank(ent);
      auto idx = destRankIdx[dr]++;
      out.permute[j] = idx;
      assert(j < out.permute.size());
      j++;
    }
  }
  return out;
}

redev::GOs getRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids)
{
  // global vertex ids
  auto gids = mesh.globals(0);
  // host side global vertex ids
  auto gids_h = Omega_h::HostRead(gids);
  // indexes of global vertex ids sorted by id
  auto iGids = sortIndexes(gids_h);
  // indexes of input
  auto iInGids = sortIndexes(inGids);
  redev::GOs rdvPermute(inGids.size());
  int j = 0;
  // for i in size of global ids
  for (size_t i = 0; i < inGids.size(); i++) {
    // loop through the global ids until get to the matching input global id
    while (gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j != gids_h.size()); // not found
    // at the
    rdvPermute[iInGids[i]] = iGids[j];
  }
  return rdvPermute;
}

CSR getRdvOutPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids)
{
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  auto iGids = sortIndexes(gids_h);
  auto iInGids = sortIndexes(inGids);
  // count the number of times each gid is included in inGids
  CSR perm;
  perm.off.resize(gids_h.size() + 1);
  int j = 0;
  for (size_t i = 0; i < inGids.size(); i++) {
    while (gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j != gids_h.size()); // found
    perm.off[iGids[j]]++;
  }
  //create the offsets array from the counts
  redev::exclusive_scan(perm.off.begin(), perm.off.end(), perm.off.begin(), 0);
  //fill the permutation array
  perm.val.resize(perm.off.back());
  redev::LOs count(gids_h.size()); // how many times each gid was written
  j = 0;
  for (size_t i = 0; i < inGids.size(); i++) {
    while (gids_h[iGids[j]] != inGids[iInGids[i]] && j < gids_h.size()) {
      j++;
    }
    REDEV_ALWAYS_ASSERT(j != gids_h.size()); // found
    const auto subIdx = count[iGids[j]]++;
    const auto startIdx = perm.off[iGids[j]];
    const auto offIdx = startIdx + subIdx;
    perm.val[offIdx] = iInGids[i];
  }
  return perm;
}

OutMsg prepareRdvOutMessage(Omega_h::Mesh& mesh,
                            const redev::InMessageLayout& in)
{
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  const auto nproc = ohComm->size();
  auto nAppProcs = Omega_h::divide_no_remainder(in.srcRanks.size(),
                                                static_cast<size_t>(nproc));
  // build dest and offsets arrays from incoming message metadata
  redev::LOs senderDeg(nAppProcs);
  for (size_t i = 0; i < nAppProcs - 1; i++) {
    senderDeg[i] =
      in.srcRanks[(i + 1) * nproc + rank] - in.srcRanks[i * nproc + rank];
  }
  const auto totInMsgs = in.offset[rank + 1] - in.offset[rank];
  senderDeg[nAppProcs - 1] =
    totInMsgs - in.srcRanks[(nAppProcs - 1) * nproc + rank];
  OutMsg out;
  for (size_t i = 0; i < nAppProcs; i++) {
    if (senderDeg[i] > 0) {
      out.dest.push_back(i);
    }
  }
  redev::GO sum = 0;
  for (auto deg : senderDeg) { // exscan over values > 0
    if (deg > 0) {
      out.offset.push_back(sum);
      sum += deg;
    }
  }
  out.offset.push_back(sum);
  return out;
}
//Omega_h::Read<Omega_h::I8> markServerOverlapRegion(
//  Omega_h::Mesh& mesh, const redev::ClassPtn& classPtn,
//  const EntInOverlapFunc& entInOverlap)
//{
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  // transfer vtx classification to host
//  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
//  auto classIds_h = Omega_h::HostRead(classIds);
//  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
//  auto classDims_h = Omega_h::HostRead(classDims);
//  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
//  Omega_h::parallel_for(
//    classIds.size(), OMEGA_H_LAMBDA(int i) {
//      isOverlap[i] = entInOverlap(classDims[i], classIds[i]);
//    });
//  auto owned_h = Omega_h::HostRead(mesh.owned(0));
//  auto isOverlap_h = Omega_h::HostRead<Omega_h::I8>(isOverlap);
//  // mask to only class partition owned entities
//  auto isOverlapOwned = Omega_h::HostWrite<Omega_h::I8>(
//    classIds.size(), "isOverlapAndOwnsModelEntInClassPartition");
//  for (int i = 0; i < mesh.nverts(); i++) {
//    redev::ClassPtn::ModelEnt ent(classDims_h[i], classIds_h[i]);
//    auto destRank = classPtn.GetRank(ent);
//    auto isModelEntOwned = (destRank == rank);
//    isOverlapOwned[i] = isModelEntOwned && isOverlap_h[i];
//    if (owned_h[i] && !isModelEntOwned) {
//      fprintf(stderr, "%d owner conflict %d ent (%d,%d) owner %d owned %d\n",
//              rank, i, classDims_h[i], classIds_h[i], destRank, owned_h[i]);
//    }
//  }
//  auto isOverlapOwned_dr = Omega_h::Read<Omega_h::I8>(isOverlapOwned);
//  // auto isOverlapOwned_hr = Omega_h::HostRead(isOverlapOwned_dr);
//  mesh.add_tag(0, "isOverlap", 1, isOverlapOwned_dr);
//  return isOverlapOwned_dr;
//}
//Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(
//  Omega_h::Mesh& mesh, const EntInOverlapFunc& entInOverlap)
//{
//  // transfer vtx classification to host
//  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
//  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
//  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
//  auto markOverlap = OMEGA_H_LAMBDA(int i)
//  {
//    isOverlap[i] = entInOverlap(classDims[i], classIds[i]);
//  };
//  Omega_h::parallel_for(classIds.size(), markOverlap);
//  auto isOwned = mesh.owned(0);
//  // try masking out to only owned entities
//  Omega_h::parallel_for(
//    isOverlap.size(),
//    OMEGA_H_LAMBDA(int i) { isOverlap[i] = (isOwned[i] && isOverlap[i]); });
//
//  auto isOverlap_r = Omega_h::read(isOverlap);
//  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
//  return isOverlap_r;
//}
redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName)
{
  namespace ts = test_support;
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank()
                               ? ts::readClassPartitionFile(cpnFileName)
                               : ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  REDEV_ALWAYS_ASSERT(mesh.nelems()); // all ranks should have elements
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

} // namespace test_support
