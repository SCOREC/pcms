#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/create_field.h"
#include <Kokkos_Core.hpp>
#include <vector>

using pcms::Real;

TEST_CASE("omega_h_field2 out of bounds FILL mode")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  // Create a 1x1 box mesh (coords from 0 to 1)
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);
  auto layout =
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  const auto nverts = mesh.nents(0);
  auto mesh_coords = mesh.coords();

  // Set up a simple linear field
  auto f = KOKKOS_LAMBDA(Real x, Real y)
  {
    return x + y;
  };
  Omega_h::Write<Real> test_f(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      Real x = mesh_coords[2 * i + 0];
      Real y = mesh_coords[2 * i + 1];
      test_f[i] = f(x, y);
    });
  Omega_h::HostWrite<Real> test_f_host(test_f);
  auto field = layout->CreateField();
  field->SetDOFHolderData(pcms::make_const_array_view(test_f_host));

  // Set FILL mode with fill value of -999.0
  Real fill_value = -999.0;
  field->SetOutOfBoundsMode(pcms::OutOfBoundsMode::FILL, fill_value);

  // Test points - mix of inside and outside
  std::vector<Real> coords = {
    0.5,  0.5,  // inside - should evaluate normally
    1.5,  0.5,  // outside (x > 1) - should return fill_value
    0.5,  -0.1, // outside (y < 0) - should return fill_value
    0.3,  0.7,  // inside - should evaluate normally
    -0.1, 0.5,  // outside (x < 0) - should return fill_value
  };

  std::vector<Real> evaluation(coords.size() / 2);
  pcms::Rank1View<Real, pcms::HostMemorySpace> eval_view{evaluation.data(),
                                                         evaluation.size()};
  pcms::Rank2View<const Real, pcms::HostMemorySpace> coords_view(
    coords.data(), coords.size() / 2, 2);
  pcms::FieldDataView<Real, pcms::HostMemorySpace> data_view(
    eval_view, field->GetCoordinateSystem());
  pcms::CoordinateView<pcms::HostMemorySpace> coordinate_view{
    field->GetCoordinateSystem(), coords_view};

  auto locale = field->GetLocalizationHint(coordinate_view);
  field->Evaluate(locale, data_view);

  // Check results
  // Point 0: (0.5, 0.5) - inside, should be close to f(0.5, 0.5) = 1.0
  REQUIRE(std::abs(evaluation[0] - 1.0) < 0.1);

  // Point 1: (1.5, 0.5) - outside, should be fill_value
  REQUIRE(evaluation[1] == fill_value);

  // Point 2: (0.5, -0.1) - outside, should be fill_value
  REQUIRE(evaluation[2] == fill_value);

  // Point 3: (0.3, 0.7) - inside, should be close to f(0.3, 0.7) = 1.0
  REQUIRE(std::abs(evaluation[3] - 1.0) < 0.1);

  // Point 4: (-0.1, 0.5) - outside, should be fill_value
  REQUIRE(evaluation[4] == fill_value);
}