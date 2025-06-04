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

void test(Omega_h::Mesh& mesh, Omega_h::Real cutoffDistance, int degree,
          Omega_h::LO min_num_supports, Omega_h::Reals source_values,
          Omega_h::Reals exact_target_values, Omega_h::Reals source_coordinates,
          Omega_h::Reals target_coordinates)
{

  int dim = mesh.dim();
  Omega_h::Real tolerance = 5e-4;

  std::vector<pcms::RadialBasisFunction> rbf_types = {
    pcms::RadialBasisFunction::RBF_GAUSSIAN, pcms::RadialBasisFunction::RBF_C4,
    pcms::RadialBasisFunction::RBF_CONST

  };

  SupportResults support =
    searchNeighbors(mesh, cutoffDistance, min_num_supports);

  for (const auto& rbf : rbf_types) {
    auto approx_target_values =
      mls_interpolation(source_values, source_coordinates, target_coordinates,
                        support, dim, degree, rbf);

    auto host_approx_target_values =
      Omega_h::HostRead<Omega_h::Real>(approx_target_values);

    auto host_exact_target_values =
      Omega_h::HostRead<Omega_h::Real>(exact_target_values);

    int m = exact_target_values.size();
    int n = approx_target_values.size();
    REQUIRE(m == n);

    for (size_t i = 0; i < m; ++i) {
      CAPTURE(i, host_exact_target_values[i], host_approx_target_values[i],
              tolerance);
      CHECK_THAT(
        host_exact_target_values[i],
        Catch::Matchers::WithinAbs(host_approx_target_values[i], tolerance));
    }
  }
}
// Test cases for centroid to node mapping using MLS
TEST_CASE("test_mls_interpolation")
{

  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  Omega_h::Real cutoffDistance = 0.3;
  cutoffDistance = cutoffDistance * cutoffDistance;

  const auto dim = mesh.dim();

  const auto& target_coordinates = mesh.coords();

  const auto& nfaces = mesh.nfaces();

  const auto& ntargets = mesh.nverts();

  Omega_h::Write<Omega_h::Real> source_coordinates(
    dim * nfaces, 0, "stores coordinates of cell centroid of each tri element");

  const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  Kokkos::parallel_for(
    "calculate the centroid in each tri element", nfaces,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      const auto current_el_verts = Omega_h::gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
        Omega_h::gather_vectors<3, 2>(target_coordinates, current_el_verts);
      auto centroid = Omega_h::average(current_el_vert_coords);
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

  SECTION("test interpolation degree 1, function degree 0")
  {
    std::cout << "-------starting test: d0p1------------" << "\n";
    int degree = 1;
    Omega_h::LO min_num_supports = 10;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));

    std::cout << " The test for d=0, p=1, passed " << "\n";
    std::cout << "----------finishing  test: d0p1-------" << "\n";
  }

  SECTION("test interpolation degree 1, function degree 1")
  {

    std::cout << "--------starting test: d1p1---------" << "\n";
    int degree = 1;
    Omega_h::LO min_num_supports = 10;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));

    std::cout << " The test for d=1, p=1, passed " << "\n";
    std::cout << "------------finishing test: d1p1--------" << "\n";
  }

  SECTION("test interpo degree 2 poly degree 0")
  {

    std::cout << "-------------starting test: d0p2------------" << "\n";
    int degree = 2;
    Omega_h::LO min_num_supports = 16;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 2);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 2);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));

    std::cout << " The test for d=0, p=2, passed " << "\n";
    std::cout << "-------------finishing test: d0p2------------" << "\n";
  }

  SECTION("test interpolation degree 2 poly degree 1")
  {

    std::cout << "----------------satrting test: d1p2--------------" << "\n";
    int degree = 2;
    Omega_h::LO min_num_supports = 16;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));

    std::cout << " The test for d=1, p=2, passed " << "\n";
    std::cout << "----------------finishing test: d1p2--------------" << "\n";
  }

  SECTION("test interpolation degree 2, function degree 2")
  {

    std::cout << "-------------starting test: d2p2--------------------" << "\n";
    int degree = 2;
    Omega_h::LO min_num_supports = 16;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));
    std::cout << " The test for d=2, p=2, passed " << "\n";
    std::cout << "-------------finishing test: d2p2--------------------"
              << "\n";
  }

  SECTION("test interpolation degree 3, function degree 2")
  {

    std::cout << "-------------starting test: d2p3--------------------" << "\n";
    int degree = 3;
    Omega_h::LO min_num_supports = 20;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree - 1);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree - 1);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));
    std::cout << " The test for d=2, p=3, passed " << "\n";
    std::cout << "-------------finishing test: d2p3--------------------"
              << "\n";
  }

  SECTION("test interpolation degree 3, function degree 3")
  {

    std::cout << "-------------starting test: d3p3--------------------" << "\n";
    int degree = 3;
    Omega_h::LO min_num_supports = 20;

    Omega_h::Write<Omega_h::Real> source_values(nfaces, 0,
                                                "exact target values");

    Kokkos::parallel_for(
      nfaces, KOKKOS_LAMBDA(int i) {
        source_values[i] = func(source_points.coordinates(i), degree);
      });

    Omega_h::Write<Omega_h::Real> exact_target_values(mesh.nverts(), 0,
                                                      "exact target values");

    Kokkos::parallel_for(
      mesh.nverts(), KOKKOS_LAMBDA(int i) {
        exact_target_values[i] = func(target_points.coordinates(i), degree);
      });

    test(mesh, cutoffDistance, degree, min_num_supports,
         Omega_h::Reals(source_values), Omega_h::Reals(exact_target_values),
         Omega_h::Reals(source_coordinates),
         Omega_h::Reals(target_coordinates));

    std::cout << " The test for d=3, p=3, passed " << "\n";
    std::cout << "-------------finishing test: d3p3--------------------"
              << "\n";
  }
}
