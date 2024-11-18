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

using namespace std;
using namespace Omega_h;

KOKKOS_INLINE_FUNCTION
double func(Coord& p, int degree)
{
  [[maybe_unused]] auto x = p.x;
  [[maybe_unused]] auto y = p.y;
  if (degree = 0) {
    double Z = 3;
  }
  elseif(degree = 1)
  {
    double Z = x + y;
  }
  elseif(degree = 2)
  {
    double Z = pow(x, 2) + pow(y, 2);
  }
  elseif(degree = 3)
  {

    double Z = pow(x, 3) + pow(y, 3);
  }
  else printf("No polynomials with degree = %f\n", degree)
}
return Z;
}

void test(Mesh& mesh, Real cutoffdistance, int degree, LO min_num_Supports,
          RadialBasisFunction bf)
{
  Write<Real> source_values(nfaces, 0, "exact target values");

  Kokkos::parallel_for(
    nfaces, KOKKOS_LAMBDA(int i) {
      source_values[i] = func(source_points.coordinates(i), degree);
    });

  SupportResults support =
    searchNeighbors(mesh, cutoffDistance, min_num_supports);

  auto approx_target_values = mls_interpolation(
    source_values, source_coordinates, target_coordinates, support, dim, degree,
    support.radii2, RadialBasisFunction::RBF_C4);

  auto host_approx_target_values = HostRead<Real>(approx_target_values);

  Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

  Kokkos::parallel_for(
    mesh.nverts(), KOKKOS_LAMBDA(int i) {
      exact_target_values[i] = func_linear(target_points.coordinates(i));
    });

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
// Test cases for centroid to node mapping using MLS
TEST_CASE("mls_interp_test")
{

  auto lib = Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  Real cutoffDistance = 0.3;
  Real tolerance = 0.05;
  cutoffDistance = cutoffDistance * cutoffDistance;

  const auto dim = mesh.dim();

  const auto& target_coordinates = mesh.coords();

  const auto& nfaces = mesh.nfaces();

  const auto& ntargets = mesh.nverts();

  Write<Real> radii2(
    ntargets, cutoffDistance,
    "populate initial square of cutoffdistance to all target points");
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
        source_values[i] = func_const(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_const(target_points.coordinates(i));
      });

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

  SECTION("test interpolation degree 1, function degree 1") {}

  SECTION("test interpo degree 2 poly degree 0")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func_const(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_const(target_points.coordinates(i));
      });

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

  SECTION("test interpolation degree 2 poly degree 1")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func_linear(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_linear(target_points.coordinates(i));
      });

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

  SECTION("test interpolation degree 2, function degree 2")
  {

    int degree = 2;
    LO min_num_supports = 16;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func_quadratic(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_quadratic(target_points.coordinates(i));
      });

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

  SECTION("test interpolation degree 3, function degree 2")
  {

    int degree = 3;
    LO min_num_supports = 20;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func_quadratic(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_quadratic(target_points.coordinates(i));
      });

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

  SECTION("test interpolation degree 3, function degree 3")
  {

    int degree = 3;
    LO min_num_supports = 20;

    Write<Real> source_values(nfaces, 0, "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func_cubic(source_points.coordinates(i));
      });

    SupportResults support =
      searchNeighbors(mesh, cutoffDistance, min_num_supports);

    auto approx_target_values = mls_interpolation(
      source_values, source_coordinates, target_coordinates, support, dim,
      degree, support.radii2, RadialBasisFunction::RBF_C4);

    auto host_approx_target_values = HostRead<Real>(approx_target_values);

    Write<Real> exact_target_values(mesh.nverts(), 0, "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func_cubic(target_points.coordinates(i));
      });

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
