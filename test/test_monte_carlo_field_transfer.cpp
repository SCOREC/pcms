#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/conservative_projection_solver.hpp>
#include <pcms/transfer/load_vector_integrator.hpp>

#include <fstream>

Kokkos::View<MeshField::Real* [3]> make_sobol_barycentric_samples()
{
  static constexpr double data[100][3] = {
    {0.77247435, 0.02761234, 0.19991331},
    {0.00066204, 0.93645417, 0.06288378},
    {0.24669859, 0.20114795, 0.55215346},
    {0.37601817, 0.44040299, 0.18357885},
    {0.45721860, 0.21639685, 0.32638456},
    {0.14873985, 0.50203922, 0.34922092},
    {0.08265261, 0.21955372, 0.69779366},
    {0.61711447, 0.30763833, 0.07524720},
    {0.55756878, 0.14595971, 0.29647151},
    {0.13181260, 0.55790360, 0.31028380},
    {0.17533046, 0.04910693, 0.77556261},
    {0.39746965, 0.60204336, 0.00048699},
    {0.32622724, 0.11952704, 0.55424572},
    {0.27416548, 0.62836549, 0.09746903},
    {0.04762209, 0.43946455, 0.51291336},
    {0.67610228, 0.17059590, 0.15330182},
    {0.70636627, 0.13868367, 0.15495006},
    {0.05423990, 0.50286400, 0.44289610},
    {0.26554739, 0.10223092, 0.63222169},
    {0.31251591, 0.56755339, 0.11993070},
    {0.44019781, 0.00429113, 0.55551107},
    {0.20250786, 0.75080234, 0.04668980},
    {0.10680470, 0.31658569, 0.57660961},
    {0.50457450, 0.32955129, 0.16587422},
    {0.59831281, 0.08088912, 0.32079807},
    {0.07370338, 0.70727787, 0.21901875},
    {0.15846319, 0.34444054, 0.49709627},
    {0.47101819, 0.31449160, 0.21449022},
    {0.34034786, 0.19216696, 0.46748519},
    {0.21550916, 0.57121237, 0.21327847},
    {0.02496562, 0.06798091, 0.90705347},
    {0.92252284, 0.06814746, 0.00932970},
    {0.85768371, 0.03637656, 0.10593973},
    {0.01619319, 0.68252604, 0.30128077},
    {0.22649636, 0.07900530, 0.69449834},
    {0.35119355, 0.59511537, 0.05369109},
    {0.48835031, 0.11281048, 0.39883921},
    {0.16631331, 0.65300768, 0.18067901},
    {0.06666197, 0.36178388, 0.57155415},
    {0.57652986, 0.24449693, 0.17897321},
    {0.52216375, 0.01874412, 0.45909213},
    {0.11478599, 0.86792052, 0.01729349},
    {0.19362412, 0.25593162, 0.55044426},
    {0.42505536, 0.36317736, 0.21176728},
    {0.30305820, 0.31329553, 0.38364627},
    {0.25341301, 0.38460776, 0.36197923},
    {0.06382536, 0.14763042, 0.78854422},
    {0.72968702, 0.22878500, 0.04152798},
    {0.65489003, 0.05202033, 0.29308965},
    {0.03867714, 0.80507099, 0.15625187},
    {0.28610086, 0.35098517, 0.36291398},
    {0.33688728, 0.36564553, 0.29746719},
    {0.41257474, 0.21920914, 0.36821612},
    {0.18330823, 0.55984932, 0.25684245},
    {0.12427171, 0.01628781, 0.85944048},
    {0.53763460, 0.44108409, 0.02128131},
    {0.64025214, 0.15451016, 0.20523770},
    {0.09045185, 0.55789331, 0.35165485},
    {0.14044321, 0.18374592, 0.67581086},
    {0.44167163, 0.43248795, 0.12584041},
    {0.36556721, 0.05186860, 0.58256418},
    {0.23481881, 0.68177609, 0.08340511},
    {0.00969458, 0.30812735, 0.68217807},
    {0.80334176, 0.14697275, 0.04968549},
    {0.79815936, 0.07689627, 0.12494437},
    {0.01472633, 0.55979307, 0.42548060},
    {0.23349852, 0.17606379, 0.59043770},
    {0.37341621, 0.49494930, 0.13163449},
    {0.45414551, 0.05244857, 0.49340592},
    {0.13694848, 0.78356864, 0.07948288},
    {0.09804517, 0.23897485, 0.66297999},
    {0.63200480, 0.25752600, 0.11046920},
    {0.54417024, 0.07635918, 0.37947058},
    {0.11645350, 0.75270655, 0.13083994},
    {0.18700282, 0.36095309, 0.45204410},
    {0.40096509, 0.30271633, 0.29631858},
    {0.32946503, 0.21998856, 0.45054641},
    {0.28751929, 0.45375342, 0.25872729},
    {0.03352118, 0.03288756, 0.93359126},
    {0.65789378, 0.33176733, 0.01033890},
    {0.74042350, 0.00658366, 0.25299284},
    {0.05660849, 0.90854533, 0.03484618},
    {0.25725536, 0.27042428, 0.47232036},
    {0.29337009, 0.47961516, 0.22701475},
    {0.41676634, 0.29019891, 0.29303475},
    {0.19477005, 0.45163734, 0.35359261},
    {0.10936446, 0.12581988, 0.76481566},
    {0.52412918, 0.39570142, 0.08016940},
    {0.57432353, 0.12821589, 0.29746058},
    {0.07183352, 0.68862912, 0.23953735},
    {0.16520587, 0.07277389, 0.76202024},
    {0.49782046, 0.45243601, 0.04974353},
    {0.36175953, 0.12969465, 0.50854582},
    {0.22280546, 0.59788723, 0.17930731},
    {0.02311096, 0.42476114, 0.55212791},
    {0.83893443, 0.10049592, 0.06056965},
    {0.91142950, 0.02213425, 0.06643625},
    {0.02989667, 0.78460762, 0.18549570},
    {0.21430992, 0.30905766, 0.47663242},
    {0.34768998, 0.37781691, 0.27449311},
    {0.48400348, 0.14291404, 0.37308247},
    {0.15508466, 0.60091608, 0.24399926}
  };

  Kokkos::View<MeshField::Real* [3]> samples("sobol_bary", 100);
  auto host = Kokkos::create_mirror_view(samples);

  for (int i = 0; i < 100; ++i) {
    host(i, 0) = data[i][0];
    host(i, 1) = data[i][1];
    host(i, 2) = data[i][2];
  }

  Kokkos::deep_copy(samples, host);
  return samples;
}

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


Omega_h::Reals evaluate_field_at_barycentric_samples(
  Omega_h::Mesh& mesh,
  const Kokkos::View<MeshField::Real*[3]>& bary,
  const std::function<double(double, double)>& f)
{
  const int npoints_each_tri = bary.extent(0); 
  const auto bary_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), bary);
  const auto coords_h = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  const auto ev2v_h = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());

  Omega_h::HostWrite<Omega_h::Real> vals_h(mesh.nelems() * npoints_each_tri);

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

      vals_h[base + i] = f(x, y);
    }
  }

  return Omega_h::Reals(vals_h);
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

  auto bary = make_sobol_barycentric_samples();
  const int npoints_each_tri = bary.extent(0);
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
