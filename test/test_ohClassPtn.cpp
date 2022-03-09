#include <cassert>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <redev_comm.h>
#include "wdmcpl.h"

void getClassPtn(Omega_h::Mesh& mesh, redev::LOs& ranks, redev::LOs& classIds) {
  auto ohComm = mesh.comm();
  const auto dim = mesh.dim();
  auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  Omega_h::ClassId max_class = Omega_h::get_max(class_ids);
  auto max_class_g = ohComm->allreduce(max_class,Omega_h_Op::OMEGA_H_MAX);
  REDEV_ALWAYS_ASSERT(ohComm->size() == max_class_g);
  auto class_ids_h = Omega_h::deep_copy(class_ids);
  if(!ohComm->rank()) {
    //send ents with classId=2 to rank 1
    Omega_h::Write<Omega_h::I32> ptnRanks(5,0);
    Omega_h::Write<Omega_h::LO> ptnIdxs(5);
    Omega_h::fill_linear(ptnIdxs,0,1);
    auto owners = Omega_h::Remotes(ptnRanks, ptnIdxs);
    mesh.migrate(owners);
  } else {
    int firstElm = 5;
    const int elms = 18; // migrating elements [5:22]
    Omega_h::Write<Omega_h::I32> ptnRanks(elms,0);
    Omega_h::Write<Omega_h::LO> ptnIdxs(elms);
    Omega_h::fill_linear(ptnIdxs,firstElm,1);
    auto owners = Omega_h::Remotes(ptnRanks, ptnIdxs);
    mesh.migrate(owners);
  }

  //the hardcoded assignment of classids to ranks
  ranks.resize(2);
  classIds.resize(2);
  ranks[0] = 0; classIds[0] = 1;
  ranks[1] = 1; classIds[1] = 2;
  //sanity check
  class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  if(!ohComm->rank()) {
    Omega_h::Read<Omega_h::ClassId> ones(class_ids.size(), 1);
    assert(class_ids == ones);
  } else {
    Omega_h::Read<Omega_h::ClassId> twos(class_ids.size(), 2);
    assert(class_ids == twos);
  }
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  int rank = world->rank();
  int nproc = world->size();
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
  redev::LOs ranks;
  redev::LOs classIds;
  if(isRdv) {
    //partition the omegah mesh by classification and return the
    //rank-to-classid array
    getClassPtn(mesh, ranks, classIds);
    REDEV_ALWAYS_ASSERT(ranks.size()==2);
    REDEV_ALWAYS_ASSERT(ranks.size()==classIds.size());
    Omega_h::vtk::write_parallel("rdvSplit.vtk", &mesh, mesh.dim());
  }
  auto ptn = redev::ClassPtn(ranks,classIds);
  redev::Redev rdv(MPI_COMM_WORLD,ptn,isRdv);
  rdv.Setup();
  std::string name = "foo";
  redev::AdiosComm<redev::LO> comm(MPI_COMM_WORLD, ranks.size(), rdv.getToEngine(), rdv.getIO(), name);
  return 0;
}
