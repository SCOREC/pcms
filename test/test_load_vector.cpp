#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_vtk.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <algorithm>
#include <string>
#include <pcms/interpolator/mesh_intersection/load_vector_integrator.hpp>

TEST_CASE("Load vector computation on intersected regions", "[load_vector]")
{

  Omega_h::Library lib;

  // 2D coordinates (x,y) : 4 vertices
  Omega_h::Reals coords({
    0.0, 0.0, // v0
    1.0, 0.0, // v1
    1.0, 1.0, // v2
    0.0, 1.0  // v3
  });

  // Target Mesh with two traingles
  // Two triangles, CCW
  // T0: (v0,v1,v3) = (0,1,3)
  // T1: (v1,v2,v3) = (1,2,3)
  Omega_h::LOs ev2v_source({0, 1, 3, 1, 2, 3});

  Omega_h::Mesh target_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&target_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_source, coords);

  // Source Mesh with two traingles
  // Two triangles, CCW
  // T0: (v0,v1,v3) = (0,1,2)
  // T1: (v1,v2,v3) = (0,2,3)
  Omega_h::LOs ev2v_target({0, 1, 2, 0, 2, 3});

  Omega_h::Mesh source_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&source_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_target, coords);

  REQUIRE(source_mesh.dim() == 2);
  REQUIRE(target_mesh.dim() == 2);

  int num_tgt_elems = target_mesh.nelems();

  auto intersection = intersectTargets(source_mesh, target_mesh);
  SECTION("check localization routine for coincident cases")
  {
    Kokkos::View<Omega_h::Real* [2]> points("test_points", 3);
    auto points_h = Kokkos::create_mirror_view(points);
    points_h(0, 0) = 0.0;
    points_h(0, 1) = 0.0;
    points_h(1, 0) = 1.0;
    points_h(0, 1) = 0.0;
    points_h(2, 0) = 1.0;
    points_h(2, 1) = 1.0;

    Kokkos::deep_copy(points, points_h);
    pcms::GridPointSearch search_cell(target_mesh, 20, 20);

    auto d_results = search_cell(points);

    auto h_results = Kokkos::create_mirror_view(d_results);

    Kokkos::deep_copy(h_results, d_results);

    REQUIRE(h_results.extent(0) == 3);

    for (int i = 0; i < h_results.extent(0); ++i) {
      auto r = h_results(i);
      std::cout << "tri_id = " << r.tri_id << "\n";
      std::cout << "dim = " << static_cast<int>(r.dimensionality) << "\n";
      // std::cout << "parametric = " << r.parametric_coords << "\n";
      INFO("tri_id : " << r.tri_id);
      REQUIRE(r.tri_id >= 0);
    }
  }

  SECTION("Basic shape and consistency of result")
  {

    // Fill dummy values at each vertex of source
    Omega_h::Write<Omega_h::Real> values(source_mesh.nverts(), 1.0);

    auto load_vector =
      pcms::buildLoadVector(target_mesh, source_mesh, intersection, values);

    auto load_vector_host = Kokkos::create_mirror(load_vector);
    Kokkos::deep_copy(load_vector_host, load_vector);

    REQUIRE(load_vector_host.extent(0) == num_tgt_elems * 3);

    for (int i = 0; i < load_vector_host.extent(0); ++i)
      REQUIRE(load_vector_host(i) >=
              0.0); // Since function is 1.0 and everything is positive
  }

  SECTION("Integration is zero if source values are zero")
  {
    Omega_h::Write<double> zero_field(source_mesh.nverts(), 0.0);

    auto load_vector =
      pcms::buildLoadVector(target_mesh, source_mesh, intersection, zero_field);

    auto load_vector_host = Kokkos::create_mirror_view(load_vector);
    Kokkos::deep_copy(load_vector_host, load_vector);

    for (int i = 0; i < load_vector_host.extent(0); ++i) {
      REQUIRE(load_vector_host(i) == Catch::Approx(0.0));
    }
  }

  SECTION("load vector computed after the intersection of simple target and "
          "source elements")
  {
    // the source elements are triangle1 (0,0), (1,0) & (1,1) and triangle2
    // (0,0), (1,1) & (0,1) the target elements are traingle1 (0,0), (1,0) &
    // (0,1) and triangle2 (1,0), (1,1) & (0,1)

    Omega_h::Write<double> constant_field(source_mesh.nverts(), 2.0);

    auto load_vector = pcms::buildLoadVector(target_mesh, source_mesh,
                                             intersection, constant_field);

    auto load_vector_host = Kokkos::create_mirror_view(load_vector);
    Kokkos::deep_copy(load_vector_host, load_vector);

    double expected_load_vector[6] = {0.333333, 0.333333, 0.333333,
                                      0.333333, 0.333333, 0.333333};
    double tolerance = 1e-6;
    for (int i = 0; i < load_vector_host.extent(0); ++i) {

      CAPTURE(i, expected_load_vector[i], load_vector_host[i], tolerance);
      CHECK_THAT(expected_load_vector[i],
                 Catch::Matchers::WithinAbs(load_vector_host[i], tolerance));
    }
  }
}
