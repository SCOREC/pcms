#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <Kokkos_Core.hpp>
#include <vector>

TEST_CASE("interpolate omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::OmegaHFieldLayout(mesh, pcms::OmegaHFieldLayoutLocation::Linear, 2);
  const auto nverts = mesh.nents(0);
  auto mesh_coords = mesh.coords();
  auto f = [](double x, double y) { return -0.3 * x + 0.5 * y; };
  Omega_h::Write<double> test_f(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      double x = mesh_coords[2 * i + 0];
      double y = mesh_coords[2 * i + 1];
      test_f[i] = f(x, y);
    });
  mesh.add_tag<double>(0, "test_f", 1, Omega_h::Read(test_f));
  pcms::OmegaHField2 field("test_f", pcms::CoordinateSystem::Cartesian, layout, mesh);
  pcms::OmegaHField2 interpolated("interpolated", pcms::CoordinateSystem::Cartesian, layout, mesh);

  std::vector<double> coords = {
    0.7681, 0.886,
    0.5337, 0.5205,
    0.8088, 0.1513,
    0.13, 0.43,
    0.5484, 0.8263,
    0.006119, 0.8642,
    0.5889, 0.5622,
    0.9268, 0.1749,
    0.2615, 0.1468,
    0.9793, 0.9612,
  };

  pcms::interpolate_field2(field, interpolated);

  std::vector<double> evaluation(coords.size() / 2);
  pcms::Rank1View<double, pcms::HostMemorySpace> eval_view{evaluation.data(),
                                                           evaluation.size()};
  pcms::Rank2View<const double, pcms::HostMemorySpace> coords_view(
    coords.data(), coords.size() / 2, 2);
  pcms::FieldDataView<double, pcms::HostMemorySpace> data_view(
    eval_view, interpolated.GetCoordinateSystem());
  pcms::CoordinateView<pcms::HostMemorySpace> coordinate_view{
    interpolated.GetCoordinateSystem(), coords_view};

  auto locale = interpolated.GetLocalizationHint(coordinate_view);
  interpolated.Evaluate(locale, data_view);

  for (int i = 0; i < coords.size() / 2; ++i) {
    double x = coords[2 * i + 0];
    double y = coords[2 * i + 1];

    CAPTURE(evaluation[i], x, y, f(x, y));
    REQUIRE(std::abs(evaluation[i] - f(x, y)) < 1e-12);
  }
}
