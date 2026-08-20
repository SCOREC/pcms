#include <pcms/utility/eqdsk.h>
#include <pcms/utility/arrays.h>
#include <pcms/field/eqdsk_field.h>
#include <pcms/field/function_space/spline.h>
#include <pcms/field/coordinate_view.hpp>
#include <pcms/field/evaluation_request.h>
#include <Kokkos_Core.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "field_test_utils.h"
#include <iostream>
#include <cmath>
#include "pcms/field/coordinate_systems/cylindrical.hpp"

using Catch::Approx;
using pcms::CoordinateSystem;
using pcms::EQDSKData;

TEST_CASE("EQDSKData manual constructor")
{
  auto lib = Omega_h::Library{};

  SECTION("Basic grid construction")
  {
    // Create a simple test grid
    // divisions stores cell counts: 4 R cells, 3 Z cells → 5x4 grid points
    pcms::Uniform2DGrid test_grid;
    test_grid.divisions[0] = 4;     // 4 R cells → 5 R points
    test_grid.divisions[1] = 3;     // 3 Z cells → 4 Z points
    test_grid.bot_left[0] = 1.0;    // R_min
    test_grid.bot_left[1] = -1.0;   // Z_min
    test_grid.edge_length[0] = 2.0; // R_max - R_min = 3.0 - 1.0
    test_grid.edge_length[1] = 2.0; // Z_max - Z_min = 1.0 - (-1.0)

    // Create test PSIZR data (5x4 = 20 values)
    std::vector<pcms::Real> test_psirz(20);
    for (int i = 0; i < 20; ++i) {
      test_psirz[i] = static_cast<pcms::Real>(i * 0.1);
    }

    // Create test psi_grid and I_psi data
    const int test_npsi = 5;
    std::vector<pcms::Real> test_psi_grid(test_npsi);
    std::vector<pcms::Real> test_I_psi(test_npsi);
    for (int i = 0; i < test_npsi; ++i) {
      test_psi_grid[i] = static_cast<pcms::Real>(i * 0.2);
      test_I_psi[i] = static_cast<pcms::Real>(1.0 + i * 0.5);
    }

    // Construct EQDSKData manually
    EQDSKData eqdsk_data(test_grid, test_psirz, test_psi_grid, test_I_psi);

    // Verify grid properties
    REQUIRE(eqdsk_data.GetNW() == 5);
    REQUIRE(eqdsk_data.GetNH() == 4);
    REQUIRE(eqdsk_data.GetRMin() == Approx(1.0));
    REQUIRE(eqdsk_data.GetRMax() == Approx(3.0));
    REQUIRE(eqdsk_data.GetZMin() == Approx(-1.0));
    REQUIRE(eqdsk_data.GetZMax() == Approx(1.0));
    REQUIRE(eqdsk_data.GetRDIM() == Approx(2.0));
    REQUIRE(eqdsk_data.GetZDIM() == Approx(2.0));

    // Verify grid spacing
    REQUIRE(eqdsk_data.GetDeltaR() == Approx(0.5));
    REQUIRE(eqdsk_data.GetDeltaZ() == Approx(2.0 / 3.0));

    // Verify PSIZR data
    REQUIRE(eqdsk_data.PSIZR.extent(0) == 20);
    auto psizr_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), eqdsk_data.PSIZR);
    for (int i = 0; i < 20; ++i) {
      REQUIRE(psizr_host(i) == Approx(test_psirz[i]));
    }

    // Verify psi_grid and I_psi data
    REQUIRE(eqdsk_data.GetNPsi() == test_npsi);
    auto psi_grid_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), eqdsk_data.psi_grid);
    auto I_psi_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), eqdsk_data.I_psi);

    for (int i = 0; i < test_npsi; ++i) {
      REQUIRE(psi_grid_host(i) == Approx(test_psi_grid[i]));
      REQUIRE(I_psi_host(i) == Approx(test_I_psi[i]));
    }
  }

  SECTION("Grid accessor methods")
  {
    // Create a test grid (4 R cells, 3 Z cells → 5x4 grid points)
    pcms::Uniform2DGrid test_grid;
    test_grid.divisions[0] = 4;
    test_grid.divisions[1] = 3;
    test_grid.bot_left[0] = 1.0;
    test_grid.bot_left[1] = -1.0;
    test_grid.edge_length[0] = 2.0;
    test_grid.edge_length[1] = 2.0;

    std::vector<pcms::Real> test_psirz(20, 0.0);
    std::vector<pcms::Real> test_psi_grid(5, 0.0);
    std::vector<pcms::Real> test_I_psi(5, 0.0);

    EQDSKData eqdsk_data(test_grid, test_psirz, test_psi_grid, test_I_psi);

    // Verify index calculations
    REQUIRE(eqdsk_data.GetPsiIndex(0, 0) == 0);
    REQUIRE(eqdsk_data.GetPsiIndex(4, 0) == 4);
    REQUIRE(eqdsk_data.GetPsiIndex(0, 3) == 15);
    REQUIRE(eqdsk_data.GetPsiIndex(4, 3) == 19);

    // Verify coordinate calculations
    REQUIRE(eqdsk_data.GetR(0) == Approx(1.0));
    REQUIRE(eqdsk_data.GetR(4) == Approx(3.0));
    REQUIRE(eqdsk_data.GetZ(0) == Approx(-1.0));
    REQUIRE(eqdsk_data.GetZ(3) == Approx(1.0));
  }

  SECTION("Input validation")
  {
    pcms::Uniform2DGrid test_grid;
    test_grid.divisions[0] = 4;
    test_grid.divisions[1] = 3;
    test_grid.bot_left[0] = 1.0;
    test_grid.bot_left[1] = -1.0;
    test_grid.edge_length[0] = 2.0;
    test_grid.edge_length[1] = 2.0;

    std::vector<pcms::Real> test_psi_grid(5, 0.0);
    std::vector<pcms::Real> test_I_psi(5, 0.0);

    // Wrong PSIZR size (should be 20)
    std::vector<pcms::Real> wrong_psirz(15, 0.0);
    REQUIRE_THROWS_AS(
      EQDSKData(test_grid, wrong_psirz, test_psi_grid, test_I_psi),
      std::runtime_error);

    // Mismatched psi_grid and I_psi sizes
    std::vector<pcms::Real> correct_psirz(20, 0.0);
    std::vector<pcms::Real> wrong_I_psi(3, 0.0);
    REQUIRE_THROWS_AS(
      EQDSKData(test_grid, correct_psirz, test_psi_grid, wrong_I_psi),
      std::runtime_error);
  }
}

TEST_CASE("EQDSKData with SplineFunctionSpace")
{
  auto lib = Omega_h::Library{};

  // Create a test grid (9 R cells, 7 Z cells → 10x8 grid points)
  pcms::Uniform2DGrid test_grid;
  test_grid.divisions[0] = 9; // 9 R cells → 10 R points
  test_grid.divisions[1] = 7; // 7 Z cells → 8 Z points
  test_grid.bot_left[0] = 1.0;
  test_grid.bot_left[1] = -1.0;
  test_grid.edge_length[0] = 2.0;
  test_grid.edge_length[1] = 2.0;

  const int NW = test_grid.divisions[0] + 1; // 10
  const int NH = test_grid.divisions[1] + 1; // 8

  // Create test PSIZR data with a simple quadratic pattern
  std::vector<pcms::Real> test_psirz(static_cast<size_t>(NW) * NH);
  for (int j = 0; j < NH; ++j) {
    for (int i = 0; i < NW; ++i) {
      double r = 1.0 + i * 2.0 / (NW - 1);
      double z = -1.0 + j * 2.0 / (NH - 1);
      test_psirz[i + j * NW] = static_cast<pcms::Real>(r * r + z * z);
    }
  }

  // Create test psi_grid and I_psi data
  const int test_npsi = 10;
  std::vector<pcms::Real> test_psi_grid(test_npsi);
  std::vector<pcms::Real> test_I_psi(test_npsi);
  for (int i = 0; i < test_npsi; ++i) {
    test_psi_grid[i] = static_cast<pcms::Real>(i * 0.5);
    test_I_psi[i] = static_cast<pcms::Real>(1.0 + i * 0.3);
  }

  EQDSKData eqdsk_data(test_grid, test_psirz, test_psi_grid, test_I_psi);

  SECTION("Spline function space creation and evaluation")
  {
    // auto spline_space = pcms::SplineFunctionSpace::FromEQDSK(eqdsk_data);

    // auto psi_field = spline_space->CreateFunction<pcms::Real>();

    auto [spline_space, psi_field] = pcms::MakeEQDSKField(eqdsk_data);
    auto layout = spline_space->GetLayout();
    const int layout_size = layout->OwnedSize();
    const int psizr_size = static_cast<int>(eqdsk_data.PSIZR.extent(0));

    REQUIRE(layout_size == psizr_size);

    // psi_field.GetData().SetDOFHolderData(
    //   pcms::Rank1View<const pcms::Real, pcms::DeviceMemorySpace>(
    //     eqdsk_data.PSIZR.data(), eqdsk_data.PSIZR.extent(0)));

    // Evaluate at test points
    const pcms::Real R_mid =
      (eqdsk_data.GetRMin() + eqdsk_data.GetRMax()) / 2.0;
    const pcms::Real Z_mid =
      (eqdsk_data.GetZMin() + eqdsk_data.GetZMax()) / 2.0;

    std::vector<pcms::Real> eval_coords = {
      R_mid,
      Z_mid, // center
      eqdsk_data.GetRMin() + 0.25 * eqdsk_data.GetRDIM(),
      eqdsk_data.GetZMin() + 0.25 * eqdsk_data.GetZDIM(), // lower left
      eqdsk_data.GetRMin() + 0.75 * eqdsk_data.GetRDIM(),
      eqdsk_data.GetZMin() + 0.75 * eqdsk_data.GetZDIM() // upper right
    };

    const int num_eval_points = 3;

    auto device_coords = pcms::test::CreateDeviceCoordinateView(
      eval_coords, pcms::csys::CylindricalRZ::Create());
    auto eval_request =
      pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view);

    auto evaluator =
      spline_space->CreatePointEvaluator<pcms::Real>(eval_request);

    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> eval_results(
      "eval_results", num_eval_points, 1);
    evaluator->Evaluate(psi_field, pcms::MakeRank2View(eval_results));

    auto results_host = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), eval_results);

    // Verify that all results are finite
    for (int i = 0; i < num_eval_points; ++i) {
      REQUIRE(std::isfinite(results_host(i, 0)));
    }
  }
}
