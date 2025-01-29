#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/adj_search.hpp>
#include <pcms/interpolator/MLS_rbf_options.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <pcms/interpolator/points.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace Omega_h;

namespace {

KOKKOS_INLINE_FUNCTION
double func(Coord& p, int degree)
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

void test(Mesh& mesh, Real cutoffDistance, int degree, LO min_num_supports,
          Reals source_values, Reals exact_target_values,
          Reals source_coordinates, Reals target_coordinates)
{

  int dim = mesh.dim();
  Real tolerance = 0.0005;

  std::vector<RadialBasisFunction> rbf_types = {
    RadialBasisFunction::RBF_GAUSSIAN, RadialBasisFunction::RBF_C4,
    RadialBasisFunction::RBF_CONST

  };

  SupportResults support =
    searchNeighbors(mesh, cutoffDistance, min_num_supports);

  for (const auto& rbf : rbf_types) {
    auto approx_target_values =
      mls_interpolation(source_values, source_coordinates, target_coordinates,
                        support, dim, degree, support.radii2, rbf);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

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
}

} // end anonymous namespace

// Test cases for centroid to node mapping using MLS
TEST_CASE("meshfields_spr_test")
{

  auto lib = Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  Real cutoffDistance = 0.3;
  cutoffDistance = cutoffDistance * cutoffDistance;

  const auto dim = mesh.dim();

  const auto& target_coordinates = mesh.coords();

  const auto& nfaces = mesh.nfaces();

  const auto& ntargets = mesh.nverts();

  Write<Real> radii2(
    ntargets, cutoffDistance,
    "populate initial square of cutoffDistance to all target points");
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

  Points source_points;

  source_points.coordinates =
    PointsViewType("Number of local source supports", nfaces);
  Kokkos::parallel_for(
    "target points", nfaces, KOKKOS_LAMBDA(int j) {
      source_points.coordinates(j).x = source_coordinates[j * dim];
      source_points.coordinates(j).y = source_coordinates[j * dim + 1];
    });

  Points target_points;

  target_points.coordinates =
    PointsViewType("Number of local source supports", mesh.nverts());
  Kokkos::parallel_for(
    "target points", mesh.nverts(), KOKKOS_LAMBDA(int j) {
      target_points.coordinates(j).x = target_coordinates[j * dim];
      target_points.coordinates(j).y = target_coordinates[j * dim + 1];
    });

  SECTION("test interpolation degree 1, function degree 0")
  {

    int degree = 1;
    LO min_num_supports = 10;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpolation degree 1, function degree 1")
  {

    int degree = 1;
    LO min_num_supports = 10;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpo degree 2 poly degree 0")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 2);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 2);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpolation degree 2 poly degree 1")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpolation degree 2, function degree 2")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpolation degree 3, function degree 2")
  {

    int degree = 3;
    LO min_num_supports = 20;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }

  SECTION("test interpolation degree 3, function degree 3")
  {

    int degree = 3;
    LO min_num_supports = 20;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports, Reals(source_values),
         Reals(exact_target_values), Reals(source_coordinates),
         Reals(target_coordinates));
  }
}
