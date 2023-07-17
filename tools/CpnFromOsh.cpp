#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_tag.hpp>
#include <array>
#include <map>
#include <numeric>

int main(int argc, char** argv)
{
  if (argc != 3) {
    printf("Usage: %s <mesh file> <nparts>\n", argv[0]);
    std::abort();
  }
  Omega_h::Library lib(&argc, &argv);
  const auto* meshFile = argv[1];
  const auto nparts = std::atoi(argv[2]);
  auto mesh = Omega_h::binary::read(meshFile, lib.world());
  auto classIds_h = Omega_h::HostRead<Omega_h::ClassId>(
    mesh.get_array<Omega_h::ClassId>(0, "class_id"));
  auto classDims_h =
    Omega_h::HostRead<Omega_h::I8>(mesh.get_array<Omega_h::I8>(0, "class_dim"));
  assert(classIds_h.size() == classDims_h.size());
  // note we use a map, because we want the geometric entities to be sorted
  // this is sensible for XGC since the entiteis are sorted.
  std::array<std::map<Omega_h::ClassId, int>, 3> geometric_ents;
  // we vouny the number of vertices in each geometric dimension/fase
  for (int i = 0; i < classDims_h.size(); ++i) {
    geometric_ents[classDims_h[i]][classIds_h[i]] += 1;
  }
  assert(geometric_ents[2].size() == 0);
  const auto nfaces = geometric_ents[1].size();
  // now lets scan the number of vertices in each face
  std::vector<Omega_h::ClassId> verts_per_face(nfaces);
  std::inclusive_scan(
    std::begin(geometric_ents[1]), std::end(geometric_ents[1]),
    std::begin(verts_per_face),
    [](auto first, auto second) { return first + second.second; },
    (Omega_h::ClassId)0);
  const auto total_verts_in_faces = verts_per_face.back();
  // roughly how many verts per part
  const auto verts_per_part = total_verts_in_faces / nparts;

  std::cout << nfaces <<"\n";
  int idx=0;
  for(int i=0; i<nfaces; ++i) {
      std::cout<<i+1<<" "<<idx<<"\n";
      if(verts_per_face[i] > verts_per_part*(idx+1)) {
        idx++;
      }
  }


  return 0;
}
