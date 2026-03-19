#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_library.hpp>
#include <array>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <pcms/transfer/load_vector_integrator.hpp>

constexpr double tol = 1e-12;

OMEGA_H_INLINE double linear_field(double x, double y)
{
  return 2.0 * x + 3.0 * y + 1.0;
}

Omega_h::Mesh make_single_triangle_mesh(Omega_h::Library& lib)
{
  // Triangle with vertices:
  // v0 = (0,0), v1 = (1,0), v2 = (0,1)
  Omega_h::Reals coords({0.0, 0.0, 1.0, 0.0, 0.0, 1.0});

  Omega_h::LOs ev2v({0, 1, 2});

  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, 2, ev2v, coords);

  return mesh;
}

Omega_h::Mesh make_unit_square_two_tri_mesh(Omega_h::Library& lib)
{
  // Square:
  // v0 = (0,0), v1 = (1,0), v2 = (1,1), v3 = (0,1)
  //
  // Two CCW triangles:
  // T0 = (0,1,3)
  // T1 = (1,2,3)
  Omega_h::Reals coords({0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0});

  Omega_h::LOs ev2v({0, 1, 3, 1, 2, 3});

  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, 2, ev2v, coords);

  return mesh;
}

Kokkos::View<MeshField::Real* [3]> make_barycentric_view(
  const std::vector<std::array<MeshField::Real, 3>>& vals)
{
  Kokkos::View<MeshField::Real* [3]> view("barycentric_coords", vals.size());
  auto host = Kokkos::create_mirror_view(view);

  for (std::size_t i = 0; i < vals.size(); ++i) {
    host(i, 0) = vals[i][0];
    host(i, 1) = vals[i][1];
    host(i, 2) = vals[i][2];
  }

  Kokkos::deep_copy(view, host);
  return view;
}

Kokkos::View<pcms::Real* [2]> make_points_view(
  const std::vector<std::array<pcms::Real, 2>>& vals)
{
  Kokkos::View<pcms::Real* [2]> view("points", vals.size());
  auto host = Kokkos::create_mirror_view(view);

  for (std::size_t i = 0; i < vals.size(); ++i) {
    host(i, 0) = vals[i][0];
    host(i, 1) = vals[i][1];
  }

  Kokkos::deep_copy(view, host);
  return view;
}

Omega_h::Reals make_linear_nodal_field_values(Omega_h::Mesh& mesh)
{
  const int dim = mesh.dim();
  const int nverts = mesh.nverts();
  auto coords = mesh.coords();

  Omega_h::Write<Omega_h::Real> values(nverts, 0.0);

  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(const int i) {
      const Omega_h::Real x = coords[i * dim + 0];
      const Omega_h::Real y = coords[i * dim + 1];
      values[i] = linear_field(x, y);
    });

  return Omega_h::read(values);
}

TEST_CASE("read_sobol_barycentric_samples_from_file reads header and values",
          "[mc_projection][sobol_io]")
{
  const std::string filename = "tmp_sobol_barycentric_samples.txt";
  std::ofstream out(filename);
  REQUIRE(out.is_open());
  out << "l0 l1 l2\n";
  out << "1.0 0.0 0.0\n";
  out << "0.0 1.0 0.0\n";
  out << "0.2 0.3 0.5\n";
  out.close();

  auto samples = pcms::read_sobol_barycentric_samples_from_file(filename);
  auto host = Kokkos::create_mirror_view(samples);
  Kokkos::deep_copy(host, samples);

  REQUIRE(host.extent(0) == 3);
  REQUIRE(host.extent(1) == 3);

  CHECK(host(0, 0) == Catch::Approx(1.0).margin(tol));
  CHECK(host(0, 1) == Catch::Approx(0.0).margin(tol));
  CHECK(host(0, 2) == Catch::Approx(0.0).margin(tol));

  CHECK(host(1, 0) == Catch::Approx(0.0).margin(tol));
  CHECK(host(1, 1) == Catch::Approx(1.0).margin(tol));
  CHECK(host(1, 2) == Catch::Approx(0.0).margin(tol));

  CHECK(host(2, 0) == Catch::Approx(0.2).margin(tol));
  CHECK(host(2, 1) == Catch::Approx(0.3).margin(tol));
  CHECK(host(2, 2) == Catch::Approx(0.5).margin(tol));

  std::remove(filename.c_str());
}

TEST_CASE(
  "generate_uniform_random_barycentric_coords returns valid barycentric rows",
  "[mc_projection][sampling]")
{
  const int npoints_each_tri = 280;

  auto bary =
    pcms::generate_uniform_random_barycentric_coords(npoints_each_tri);
  auto host = Kokkos::create_mirror_view(bary);
  Kokkos::deep_copy(host, bary);

  REQUIRE(host.extent(0) == npoints_each_tri);
  REQUIRE(host.extent(1) == 3);

  for (int i = 0; i < npoints_each_tri; ++i) {
    const double l0 = host(i, 0);
    const double l1 = host(i, 1);
    const double l2 = host(i, 2);

    CHECK(l0 >= 0.0);
    CHECK(l1 >= 0.0);
    CHECK(l2 >= 0.0);

    CHECK(l0 <= 1.0);
    CHECK(l1 <= 1.0);
    CHECK(l2 <= 1.0);

    CHECK(l0 + l1 + l2 == Catch::Approx(1.0).margin(1e-12));
  }
}

TEST_CASE(
  "global_coords_from_ref_barycentric_coords maps reference samples correctly",
  "[mc_projection][mapping]")
{
  Omega_h::Library lib;
  auto mesh = make_single_triangle_mesh(lib);

  auto ref_barycentric_coords = make_barycentric_view({{1.0, 0.0, 0.0},
                                                       {0.0, 1.0, 0.0},
                                                       {0.0, 0.0, 1.0},
                                                       {0.5, 0.5, 0.0},
                                                       {0.25, 0.25, 0.50}});

  auto global_points = pcms::global_coords_from_ref_barycentric_coords(
    mesh, ref_barycentric_coords);

  auto host_points = Kokkos::create_mirror_view(global_points);
  Kokkos::deep_copy(host_points, global_points);

  REQUIRE(host_points.extent(0) == 5);
  REQUIRE(host_points.extent(1) == 2);

  // (1,0,0) -> v0 = (0,0)
  CHECK(host_points(0, 0) == Catch::Approx(0.0).margin(tol));
  CHECK(host_points(0, 1) == Catch::Approx(0.0).margin(tol));

  // (0,1,0) -> v1 = (1,0)
  CHECK(host_points(1, 0) == Catch::Approx(1.0).margin(tol));
  CHECK(host_points(1, 1) == Catch::Approx(0.0).margin(tol));

  // (0,0,1) -> v2 = (0,1)
  CHECK(host_points(2, 0) == Catch::Approx(0.0).margin(tol));
  CHECK(host_points(2, 1) == Catch::Approx(1.0).margin(tol));

  // (0.5,0.5,0) -> (0.5,0)
  CHECK(host_points(3, 0) == Catch::Approx(0.5).margin(tol));
  CHECK(host_points(3, 1) == Catch::Approx(0.0).margin(tol));

  // (0.25,0.25,0.5) -> (0.25,0.5)
  CHECK(host_points(4, 0) == Catch::Approx(0.25).margin(tol));
  CHECK(host_points(4, 1) == Catch::Approx(0.50).margin(tol));
}

TEST_CASE("localize_points_in_mesh identifies inside and outside points",
          "[mc_projection][localization]")
{
  Omega_h::Library lib;
  auto mesh = make_unit_square_two_tri_mesh(lib);

  auto points = make_points_view({
    {0.2, 0.2}, // inside
    {0.8, 0.8}, // inside
    {1.2, 0.5}  // outside
  });

  auto results = pcms::localize_points_in_mesh(mesh, points);
  auto host_results = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(host_results, results);

  REQUIRE(host_results.extent(0) == 3);

  CHECK(host_results(0).element_id >= 0);
  CHECK(host_results(1).element_id >= 0);
  CHECK(host_results(2).element_id < 0);
}

TEST_CASE(
  "evaluate_field_from_point_localization reproduces linear field exactly",
  "[mc_projection][field_eval]")
{
  Omega_h::Library lib;
  auto mesh = make_single_triangle_mesh(lib);

  // For f(x,y) = 2x + 3y + 1:
  // at (0,0) -> 1
  // at (1,0) -> 3
  // at (0,1) -> 4
  Omega_h::Reals nodal_field_values({1.0, 3.0, 4.0});

  auto points = make_points_view({{0.2, 0.2}, {0.1, 0.7}, {0.3, 0.1}});

  auto results = pcms::localize_points_in_mesh(mesh, points);

  auto field_values = pcms::evaluate_field_from_point_localization(
    mesh, nodal_field_values, results);

  auto host_points = Kokkos::create_mirror_view(points);
  Kokkos::deep_copy(host_points, points);

  auto host_values = Omega_h::HostRead<Omega_h::Real>(field_values);

  for (int i = 0; i < 3; ++i) {
    const double x = host_points(i, 0);
    const double y = host_points(i, 1);
    const double exact = linear_field(x, y);

    CHECK(host_values[i] == Catch::Approx(exact).margin(1e-12));
  }
}

TEST_CASE("evaluate_field_from_point_localization works on multi-element mesh",
          "[mc_projection][field_eval][multi_elem]")
{
  Omega_h::Library lib;
  auto mesh = make_unit_square_two_tri_mesh(lib);

  auto nodal_field_values = make_linear_nodal_field_values(mesh);

  auto points = make_points_view({{0.2, 0.2}, {0.8, 0.8}, {0.75, 0.25}});

  auto results = pcms::localize_points_in_mesh(mesh, points);

  auto field_values = pcms::evaluate_field_from_point_localization(
    mesh, nodal_field_values, results);

  auto host_points = Kokkos::create_mirror_view(points);
  Kokkos::deep_copy(host_points, points);

  auto host_values = Omega_h::HostRead<Omega_h::Real>(field_values);

  for (int i = 0; i < 3; ++i) {
    const double x = host_points(i, 0);
    const double y = host_points(i, 1);
    const double exact = linear_field(x, y);

    CHECK(host_values[i] == Catch::Approx(exact).margin(1e-12));
  }
}

TEST_CASE("sampling mapping localization and field evaluation pipeline works",
          "[mc_projection][integration]")
{
  Omega_h::Library lib;
  auto source_mesh = make_unit_square_two_tri_mesh(lib);
  auto target_mesh = make_single_triangle_mesh(lib);

  auto ref_barycentric_coords = make_barycentric_view(
    {{0.6, 0.2, 0.2}, {0.2, 0.6, 0.2}, {0.2, 0.2, 0.6}, {0.3, 0.3, 0.4}});

  auto sampled_points = pcms::global_coords_from_ref_barycentric_coords(
    target_mesh, ref_barycentric_coords);

  auto results = pcms::localize_points_in_mesh(source_mesh, sampled_points);

  auto nodal_field_values = make_linear_nodal_field_values(source_mesh);

  auto sampled_field_values = pcms::evaluate_field_from_point_localization(
    source_mesh, nodal_field_values, results);

  auto host_points = Kokkos::create_mirror_view(sampled_points);
  Kokkos::deep_copy(host_points, sampled_points);

  auto host_field_values =
    Omega_h::HostRead<Omega_h::Real>(sampled_field_values);

  for (int i = 0; i < 4; ++i) {
    const double x = host_points(i, 0);
    const double y = host_points(i, 1);
    const double exact = linear_field(x, y);

    CHECK(host_field_values[i] == Catch::Approx(exact).margin(1e-12));
  }
}
