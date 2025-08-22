#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <Kokkos_Core.hpp>
#include <vector>

TEST_CASE("interpolate linear 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::OmegaHFieldLayout(mesh, {1, 0, 0, 0}, 1, pcms::CoordinateSystem::Cartesian);
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
  pcms::OmegaHField2 field("", layout, mesh);
  pcms::OmegaHField2 interpolated("", layout, mesh);
  pcms::Rank1View<const double, pcms::HostMemorySpace> array_view{
    std::data(test_f), std::size(test_f)};
  pcms::FieldDataView<const double, pcms::HostMemorySpace> field_data_view{
    array_view, field.GetCoordinateSystem()};
  field.SetDOFHolderData(field_data_view);

  pcms::interpolate_field2(field, interpolated);
  auto interpolated_dof = interpolated.GetDOFHolderData();
  auto original_dof = field.GetDOFHolderData();
  REQUIRE(interpolated_dof.GetCoordinateSystem() == original_dof.GetCoordinateSystem());
  REQUIRE( interpolated_dof.Size() == original_dof.Size());
  // assumes that GetDOFHolderData will return a host view
  for (int i = 0; i < interpolated_dof.Size(); ++i) {
    REQUIRE_THAT(interpolated_dof.GetValues()[i],
                 Catch::Matchers::WithinRel(original_dof.GetValues()[i], 0.001) ||
                 Catch::Matchers::WithinAbs(original_dof.GetValues()[i], 1E-10)
                 );
  }
}

TEST_CASE("interpolate quadratic 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::OmegaHFieldLayout(mesh, {1, 1, 0, 0}, 1, pcms::CoordinateSystem::Cartesian);
  const auto nverts = mesh.nents(0);
  const auto nedges = mesh.nents(1);
  auto mesh_coords = mesh.coords();
  auto edge_verts = mesh.ask_verts_of(1);
  auto f = [](double x, double y) { return -0.3 * x + 0.5 * y; };
  Omega_h::Write<double> test_f(nverts + nedges);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      double x = mesh_coords[2 * i + 0];
      double y = mesh_coords[2 * i + 1];
      test_f[i] = f(x, y);
    });
  Omega_h::parallel_for(
    nedges, OMEGA_H_LAMBDA(int i) {
      auto endpoints = Omega_h::gather_verts<2>(edge_verts, i);
      double x0 = mesh_coords[2 * endpoints[0] + 0];
      double y0 = mesh_coords[2 * endpoints[0] + 1];
      double x1 = mesh_coords[2 * endpoints[1] + 0];
      double y1 = mesh_coords[2 * endpoints[1] + 1];
      double cx = (x0 + x1) / 2;
      double cy = (y0 + y1) / 2;
      test_f[nverts + i] = f(cx, cy);
    });

  pcms::OmegaHField2 field("", layout, mesh);
  pcms::OmegaHField2 interpolated("", layout, mesh);
  pcms::Rank1View<const double, pcms::HostMemorySpace> array_view{
    std::data(test_f), std::size(test_f)
  };
  field.SetDOFHolderData({array_view, field.GetCoordinateSystem()});

  // interpolate the field from one mesh to another mesh with the same coordinates
  pcms::interpolate_field2(field, interpolated);

  auto interpolated_dof = interpolated.GetDOFHolderData();
  auto original_dof = field.GetDOFHolderData();
  REQUIRE(interpolated_dof.GetCoordinateSystem() == original_dof.GetCoordinateSystem());
  REQUIRE( interpolated_dof.Size() == original_dof.Size());
  // assumes that GetDOFHolderData will return a host view
  for (int i = 0; i < interpolated_dof.Size(); ++i) {
    REQUIRE_THAT(interpolated_dof.GetValues()[i],
                 Catch::Matchers::WithinRel(original_dof.GetValues()[i], 0.001) ||
                 Catch::Matchers::WithinAbs(original_dof.GetValues()[i], 1E-10)
                 );
  }
}
