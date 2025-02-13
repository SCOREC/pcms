#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/adj_search.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace Omega_h;

namespace {

void print_patches(Omega_h::Graph& patches) {
  auto offsets = HostRead(patches.a2ab);
  auto values = HostRead(patches.ab2b);
  std::cout << "num patches " << patches.nnodes() << "\n";
  for(int patch=0; patch<patches.nnodes(); patch++) {
    std::cout << "patch " << patch << " patchElms ";
    for (auto valIdx = offsets[patch]; valIdx < offsets[patch + 1]; ++valIdx) {
      auto patchElm = values[valIdx];
      std::cout << patchElm << " ";
    }
    std::cout << "\n";
  }

}

std::string getTestName(size_t interp_degree, size_t func_degree) {
  return "test interpolation degree " + std::to_string(interp_degree) +
    ", polynomial degree " + std::to_string(func_degree);
}

KOKKOS_INLINE_FUNCTION
double func(pcms::Coord& p, int degree)
{
  [[maybe_unused]] auto x = p.x;
  [[maybe_unused]] auto y = p.y;
  if (degree == 0) {
    return 3;
  } else if (degree == 1) {
    return x + y;
  } else if (degree == 2) {
    return pow(x, 2) + pow(y, 2);
  } else if (degree == 3) {

    return pow(x, 3) + pow(y, 3);
  } else {
    printf("No polynomials with degree = %d\n", degree);
  }
  return -1;
}

void test(Mesh& mesh, Omega_h::Graph& patches, int degree,
          Reals source_values, Reals exact_target_values,
          Reals source_coordinates, Reals target_coordinates)
{

  int dim = mesh.dim();
  Real tolerance = 0.0005;

  Omega_h::Write<Real> ignored(patches.ab2b.size(), 1);
  SupportResults support{patches.a2ab,patches.ab2b,ignored};

  auto approx_target_values =
    pcms::mls_interpolation (source_values, source_coordinates, target_coordinates,
        support, dim, degree, support.radii2, pcms::RadialBasisFunction::NO_OP);

  auto host_approx_target_values = HostRead<Real>(approx_target_values);
  mesh.add_tag(OMEGA_H_VERT, "approx_target_values", 1, read(approx_target_values));
  Omega_h::vtk::write_parallel("box.vtk", &mesh);

  auto host_exact_target_values = HostRead<Real>(exact_target_values);

  int m = exact_target_values.size();
  int n = approx_target_values.size();

  REQUIRE(m == n);

  for (size_t i = 0; i < m; ++i) {
    CHECK_THAT(
        host_exact_target_values[i],
        Catch::Matchers::WithinAbs(host_approx_target_values[i], tolerance));
  }
}

} // end anonymous namespace

// Test cases for centroid to node mapping using MLS
TEST_CASE("meshfields_spr_test")
{

  auto lib = Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();
  const auto boxSize = 1.0;
  const auto nElms = 6;
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, boxSize, boxSize, 0, nElms, nElms, 0, false);
  std::cout << "mesh: elms " << mesh.nelems() << " verts " << mesh.nverts() << "\n";

  const auto dim = mesh.dim();

  const auto& target_coordinates = mesh.coords();

  const auto& nfaces = mesh.nfaces();

  const auto& ntargets = mesh.nverts();

  Write<Real> source_coordinates(
    dim * nfaces, 0, "stores coordinates of cell centroid of each tri element");

  const auto& faces2nodes = mesh.ask_down(FACE, VERT).ab2b;

  Kokkos::parallel_for(
    "calculate the centroid in each tri element", nfaces,
    OMEGA_H_LAMBDA(const LO id) {
      const auto current_el_verts = gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
        gather_vectors<3, 2>(target_coordinates, current_el_verts);
      auto centroid = average(current_el_vert_coords);
      int index = 2 * id;
      source_coordinates[index] = centroid[0];
      source_coordinates[index + 1] = centroid[1];
    });

  pcms::Points source_points;

  source_points.coordinates =
    pcms::PointsViewType("Number of local source supports", nfaces);
  Kokkos::parallel_for(
    "target points", nfaces, KOKKOS_LAMBDA(int j) {
      source_points.coordinates(j).x = source_coordinates[j * dim];
      source_points.coordinates(j).y = source_coordinates[j * dim + 1];
    });

  pcms::Points target_points;

  target_points.coordinates =
    pcms::PointsViewType("Number of local source supports", mesh.nverts());
  Kokkos::parallel_for(
    "target points", mesh.nverts(), KOKKOS_LAMBDA(int j) {
      target_points.coordinates(j).x = target_coordinates[j * dim];
      target_points.coordinates(j).y = target_coordinates[j * dim + 1];
    });


  //Define tests
  for(size_t interp_degree=1; interp_degree<=3; interp_degree++) {
    for(int func_degree=interp_degree; func_degree>=0; func_degree--) {
      SECTION(getTestName(interp_degree, func_degree).c_str()) {
        std::cerr << "start " << interp_degree << ", " << func_degree << " \n";
        const auto minPatchSize = interp_degree+2;
        std::cerr << "minPatchSize " << minPatchSize << "\n";
        auto patches = mesh.get_vtx_patches(minPatchSize);
        print_patches(patches);

        Write<Real> source_values(nfaces, 0, "exact target values");

        Kokkos::parallel_for(
            nfaces, KOKKOS_LAMBDA(int i) {
            source_values[i] = func(source_points.coordinates(i), func_degree);
            });
        mesh.add_tag(mesh.dim(), "source_values", 1, read(source_values));

        Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

        Kokkos::parallel_for(
            mesh.nverts(), KOKKOS_LAMBDA(int i) {
            exact_target_values[i] = func(target_points.coordinates(i), func_degree);
            });
        mesh.add_tag(OMEGA_H_VERT, "target_values", 1, read(exact_target_values));

        test(mesh, patches, interp_degree, Reals(source_values),
            Reals(exact_target_values), Reals(source_coordinates),
            Reals(target_coordinates));
        std::cerr << "done " << interp_degree << ", " << func_degree << " \n";
      } //end SECTION
    }
  }
}
