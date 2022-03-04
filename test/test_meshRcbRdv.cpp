#include <cassert>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include "wdmcpl.h"

int main(int argc, char** argv) {
  int rank, nproc;
  auto lib = Omega_h::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 2);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  const auto expectedRanks = redev::LOs({0,1,2,3});
  const auto expectedCuts = redev::Reals({0,0.5,0.75,0.25});
  auto ranks = rank==0 ? expectedRanks : redev::LOs(4);
  auto cuts = rank==0 ? expectedCuts : redev::Reals(4);
  const auto dim = 2;
  auto ptn = redev::RCBPtn(dim,ranks,cuts);
  ptn.Broadcast(MPI_COMM_WORLD);
  assert(ptn.GetRanks() == expectedRanks);
  assert(ptn.GetCuts() == expectedCuts);
  return 0;
}
