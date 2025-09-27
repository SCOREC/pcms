#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/create_field.h"
#include <Kokkos_Core.hpp>
#include <vector>

TEST_CASE("evaluate linear 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  const auto nverts = mesh.nents(0);
  auto mesh_coords = mesh.coords();
  auto f = [](double x, double y) { return std::sin(20 * x * y) / 2 + 0.5; };
  Omega_h::Write<double> test_f(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      double x = mesh_coords[2 * i + 0];
      double y = mesh_coords[2 * i + 1];
      test_f[i] = f(x, y);
    });
  auto field = layout->CreateField();
  field->SetDOFHolderData(pcms::make_const_array_view(test_f));

  std::vector<double> coords = {
    0.7681, 0.886,  0.5337, 0.5205,   0.8088, 0.1513, 0.13,
    0.43,   0.5484, 0.8263, 0.006119, 0.8642, 0.5889, 0.5622,
    0.9268, 0.1749, 0.2615, 0.1468,   0.9793, 0.9612,
  };

  std::vector<double> evaluation(coords.size() / 2);
  pcms::Rank1View<double, pcms::HostMemorySpace> eval_view{evaluation.data(),
                                                           evaluation.size()};
  pcms::Rank2View<const double, pcms::HostMemorySpace> coords_view(
    coords.data(), coords.size() / 2, 2);
  pcms::FieldDataView<double, pcms::HostMemorySpace> data_view(
    eval_view, field->GetCoordinateSystem());
  pcms::CoordinateView<pcms::HostMemorySpace> coordinate_view{
    field->GetCoordinateSystem(), coords_view};

  auto locale = field->GetLocalizationHint(coordinate_view);
  field->Evaluate(locale, data_view);

  for (int i = 0; i < coords.size() / 2; ++i) {
    double x = coords[2 * i + 0];
    double y = coords[2 * i + 1];

    double test_value = evaluation[i];
    double reference_value = f(x, y);
    double percent_error =
      100 * std::abs(test_value - reference_value) / reference_value;

    REQUIRE(percent_error < 1.0);
  }
}

TEST_CASE("evaluate quadratic 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto layout =
    pcms::CreateLagrangeLayout(mesh, 2, 1, pcms::CoordinateSystem::Cartesian);
  const auto nverts = mesh.nents(0);
  const auto nedges = mesh.nents(1);
  auto mesh_coords = mesh.coords();
  auto edge_verts = mesh.ask_verts_of(1);
  auto f = [](double x, double y) { return std::sin(20 * x * y) / 2 + 0.5; };
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

  auto field = layout->CreateField();
  field->SetDOFHolderData(pcms::make_const_array_view(test_f));

  std::vector<double> coords = {
    0.7681, 0.886,  0.5337, 0.5205,   0.8088, 0.1513, 0.13,
    0.43,   0.5484, 0.8263, 0.006119, 0.8642, 0.5889, 0.5622,
    0.9268, 0.1749, 0.2615, 0.1468,   0.9793, 0.9612,
  };

  std::vector<double> evaluation(coords.size() / 2);
  pcms::Rank1View<double, pcms::HostMemorySpace> eval_view{evaluation.data(),
                                                           evaluation.size()};
  pcms::Rank2View<const double, pcms::HostMemorySpace> coords_view(
    coords.data(), coords.size() / 2, 2);
  pcms::FieldDataView<double, pcms::HostMemorySpace> data_view(
    eval_view, field->GetCoordinateSystem());
  pcms::CoordinateView<pcms::HostMemorySpace> coordinate_view{
    field->GetCoordinateSystem(), coords_view};

  auto locale = field->GetLocalizationHint(coordinate_view);
  field->Evaluate(locale, data_view);

  for (int i = 0; i < coords.size() / 2; ++i) {
    double x = coords[2 * i + 0];
    double y = coords[2 * i + 1];

    double test_value = evaluation[i];
    double reference_value = f(x, y);
    double percent_error =
      100 * std::abs(test_value - reference_value) / reference_value;

    REQUIRE(percent_error < 1.0);
  }
}
