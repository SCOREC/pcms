#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array.hpp>


using Omega_h::Read;
using Omega_h::Write;
using Omega_h::Real;
using Omega_h::Reals;
using Omega_h::Mesh;
using Omega_h::LOs;
using Omega_h::LO;
using Omega_h::parallel_for;
using Omega_h::gather_verts;
using Omega_h::gather_vectors;
using Omega_h::gather_scalars;


// should have two interpolate overloads.
// 1) for scalar...all input data evaluated at same isoparametric point
// 1) for vector...all input data evaluated at corresponding isoparametric point
Reals interpolate(Mesh& mesh, Reals field, const LOs& elements, Omega_h::Vector<2>&& coord) {
  Write<Real> result(field.size());
  const auto tris2verts = mesh.ask_elem_verts();
  const auto coords = mesh.coords();
  parallel_for("", elements.size(), OMEGA_H_LAMBDA(LO i){
    const auto current_tri2verts =  gather_verts<3>(tris2verts, i);
    // 2d mesh with 2d coords, but 3 triangles
    const auto vertex_coords = gather_vectors<3, 2>(coords, current_tri2verts);
    const auto vertex_field = gather_scalars<3>(field, current_tri2verts);
    const auto xi0 = coord(0);
    const auto xi1 = coord(1);
    const auto xi2 = 1-coord[0]-coord[1];
    result[i] = vertex_field(0)*xi0+vertex_field(1)*xi1+vertex_field(2)*xi2;
  });
  return result;
}






int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::Mesh(&lib);
  const auto world = lib.world();
  LOs elements{1,2,5,7};
  Omega_h::Vector<2>  coord{0.1,0.1};
  std::string field_name = "field";
  auto field = mesh.get_array<Real>(2, field_name);
  auto result = interpolate(mesh, field, elements, {0.1,0.1});
}