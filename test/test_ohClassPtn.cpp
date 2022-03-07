#include <cassert>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include "wdmcpl.h"

void getClassPtn(Omega_h::Mesh& mesh, redev::LOs ranks, redev::LOs classIds) {
  auto ohComm = mesh.comm();
  auto class_ids = mesh.get_array<Omega_h::ClassId>(dim, "class_id");
  Omega_h::ClassId max_class = Omega_h::get_max(class_ids);
  auto max_class_g = ohComm->allreduce(max_class);
  fprintf(stderr, "max_class_g %d\n", max_class_g);
}

int main(int argc, char** argv) {
  int rank, nproc;
  auto lib = Omega_h::Library(&argc, &argv);
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh\n";
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3);
  auto isRdv = atoi(argv[1]);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  redev::LOs ranks;
  redev::LOs classIds;
  getClassPtn(mesh, ranks, classIds);
//  const auto isRdv = true;
//  redev::Redev rdv(MPI_COMM_WORLD,ptn,isRdv);
//  rdv.Setup();
//  std::string name = "foo";
//  redev::AdiosComm<redev::LO> comm(MPI_COMM_WORLD, ranks.size(), rdv.getToEngine(), rdv.getIO(), name);
  return 0;
}
