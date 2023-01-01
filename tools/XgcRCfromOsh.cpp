#include <Omega_h_file.hpp>
#include <wdmcpl/xgc_reverse_classification.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_tag.hpp>

int main(int argc, char** argv)
{
  if ((argc != 2) && (argc != 3)) {
    printf("Usage: %s <mesh file> <numbering>\n", argv[0]);
    std::abort();
  }
  Omega_h::Library lib(&argc, &argv);
  const auto* meshFile = argv[1];
  auto mesh = Omega_h::binary::read(meshFile, lib.world());
  std::string numbering;
  if (argc == 2) {
    numbering = "simNumbering";
  } else {
    numbering = std::string(argv[2]);
  }
  const auto* tag = mesh.get_tagbase(0, numbering);
  auto index_base = wdmcpl::IndexBase::Zero;
  if (numbering == "simNumbering") {
    index_base = wdmcpl::IndexBase::One;
  }
  if (Omega_h::is<Omega_h::LO>(tag)) {
    std::cout << wdmcpl::ConstructRCFromOmegaHMesh<Omega_h::LO>(mesh, numbering,
                                                                index_base);
  } else if (Omega_h::is<Omega_h::GO>(tag)) {
    std::cout << wdmcpl::ConstructRCFromOmegaHMesh<Omega_h::GO>(mesh, numbering,
                                                                index_base);
  } else {
    std::cerr << "IDs should be either be LO or GO\n";
    std::abort();
  }
  return 0;
}