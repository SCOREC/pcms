#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/conservative_projection_solver.hpp>
#include <pcms/transfer/load_vector_integrator.hpp>

#include <fstream>

double integrate_linear_field(Omega_h::Mesh& mesh, const Omega_h::Reals& u)
{
  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_verts = mesh.ask_elem_verts();
  const auto nverts = mesh.nverts();

  REQUIRE(static_cast<Omega_h::LO>(u.size()) == nverts);

  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);
  const auto elem_verts_h = Omega_h::HostRead<Omega_h::LO>(elem_verts);
  const auto u_h = Omega_h::HostRead<Omega_h::Real>(u);

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const Omega_h::LO v0 = elem_verts_h[3 * e + 0];
    const Omega_h::LO v1 = elem_verts_h[3 * e + 1];
    const Omega_h::LO v2 = elem_verts_h[3 * e + 2];

    const double area = elem_areas_h[e];
    const double avg = (u_h[v0] + u_h[v1] + u_h[v2]) / 3.0;
    integral += area * avg;
  }
  return integral;
}

void write_barycentric_samples_file(
  const std::string& filename,
  const std::vector<std::array<double, 3>>& bary)
{
  std::ofstream out(filename);
  REQUIRE(out.is_open());
  out << "l0 l1 l2\n";
  for (const auto& s : bary) {
    out << s[0] << " " << s[1] << " " << s[2] << "\n";
  }
  out.close();
}

Omega_h::Reals evaluate_field_at_barycentric_samples(
  Omega_h::Mesh& mesh,
  const std::vector<std::array<double, 3>>& bary,
  const std::function<double(double, double)>& f)
{
  const int npoints_each_tri = static_cast<int>(bary.size());
  const auto coords_h = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  const auto ev2v_h = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());

  Omega_h::Write<Omega_h::Real> vals(mesh.nelems() * npoints_each_tri);

  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const Omega_h::LO v0 = ev2v_h[3 * e + 0];
    const Omega_h::LO v1 = ev2v_h[3 * e + 1];
    const Omega_h::LO v2 = ev2v_h[3 * e + 2];

    const double x0 = coords_h[2 * v0 + 0];
    const double y0 = coords_h[2 * v0 + 1];
    const double x1 = coords_h[2 * v1 + 0];
    const double y1 = coords_h[2 * v1 + 1];
    const double x2 = coords_h[2 * v2 + 0];
    const double y2 = coords_h[2 * v2 + 1];

    const int base = e * npoints_each_tri;
    for (int i = 0; i < npoints_each_tri; ++i) {
      const double l0 = bary[i][0];
      const double l1 = bary[i][1];
      const double l2 = bary[i][2];

      const double x = l0 * x0 + l1 * x1 + l2 * x2;
      const double y = l0 * y0 + l1 * y1 + l2 * y2;

      vals[base + i] = f(x, y);
    }
  }

  return Omega_h::Reals(vals);
}


TEST_CASE("monte carlo projection preserves constant and linear fields",
          "[transfer][mc_projection]")
{
  Omega_h::Library lib;

  Omega_h::Reals coords({
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0
  });

  Omega_h::LOs ev2v_target({0, 1, 3, 1, 2, 3});
  Omega_h::Mesh target_mesh(&lib);
  Omega_h::build_from_elems_and_coords(
    &target_mesh, OMEGA_H_SIMPLEX, 2, ev2v_target, coords);

  // Deterministic barycentric samples for the reference triangle.
  // Replace this with many more Sobol samples in real application
  std::vector<std::array<double, 3>> bary;
  bary.reserve(8);
  bary.push_back({1.0/3.0, 1.0/3.0, 1.0/3.0});
  bary.push_back({0.6, 0.2, 0.2});
  bary.push_back({0.2, 0.6, 0.2});
  bary.push_back({0.2, 0.2, 0.6});
  bary.push_back({0.5, 0.4, 0.1});
  bary.push_back({0.5, 0.1, 0.4});
  bary.push_back({0.1, 0.5, 0.4});
  bary.push_back({0.1, 0.4, 0.5});

  const std::string sobol_filename = "tmp_mc_barycentric_samples.txt";
  write_barycentric_samples_file(sobol_filename, bary);

  const int npoints_each_tri = static_cast<int>(bary.size());

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;

    auto field_values_at_points = evaluate_field_at_barycentric_samples(
      target_mesh, bary, [c](double, double) { return c; });

    auto projected = pcms::solveGalerkinProjectionMC(
      target_mesh, field_values_at_points, npoints_each_tri,
      pcms::SamplingMethod::SOBOL,sobol_filename);

    auto projected_h = Omega_h::HostRead<Omega_h::Real>(projected);

    REQUIRE(static_cast<Omega_h::LO>(projected.size()) == target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(projected_h[i] == Catch::Approx(c).margin(1e-10));
    }

    const double exact_integral = c * 1.0; // unit square area = 1
    const double projected_integral =
      integrate_linear_field(target_mesh, projected);
    REQUIRE(projected_integral == Catch::Approx(exact_integral).margin(1e-10));
  }

  SECTION("linear field is reproduced approximately")
  {
    auto field_values_at_points = evaluate_field_at_barycentric_samples(
      target_mesh, bary, [](double x, double y) { return x + y; });

    auto projected = pcms::solveGalerkinProjectionMC(
      target_mesh, field_values_at_points, npoints_each_tri,
      pcms::SamplingMethod::SOBOL, sobol_filename);

    auto projected_h = Omega_h::HostRead<Omega_h::Real>(projected);
    auto tgt_coords_h = Omega_h::HostRead<Omega_h::Real>(target_mesh.coords());

    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h[2 * i + 0] + tgt_coords_h[2 * i + 1];
      REQUIRE(projected_h[i] == Catch::Approx(expected).margin(1e-6));
    }

    const double exact_integral = 1.0; // integral of x+y over unit square
    const double projected_integral =
      integrate_linear_field(target_mesh, projected);
    REQUIRE(projected_integral == Catch::Approx(exact_integral).margin(1e-6));
  }
}
