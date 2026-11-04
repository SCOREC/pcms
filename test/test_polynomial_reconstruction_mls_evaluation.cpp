#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "field_test_utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"

using pcms::CoordinateSystem;
using pcms::CoordinateView;
using pcms::DeviceMemorySpace;
using pcms::HostMemorySpace;
using pcms::Rank1View;
using pcms::Rank2View;
using pcms::Real;

namespace
{

// Build a regular NxN grid of source points over [0,1]^2.
std::vector<Real> MakeGrid2D(int N)
{
  std::vector<Real> pts;
  pts.reserve(static_cast<size_t>(N) * N * 2);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) {
      pts.push_back(static_cast<Real>(i) / (N - 1));
      pts.push_back(static_cast<Real>(j) / (N - 1));
    }
  }
  return pts;
}

std::vector<Real> MakeGrid3D(int N)
{
  std::vector<Real> pts;
  pts.reserve(static_cast<size_t>(N) * N * N * 3);
  for (int k = 0; k < N; ++k) {
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        pts.push_back(static_cast<Real>(i) / (N - 1));
        pts.push_back(static_cast<Real>(j) / (N - 1));
        pts.push_back(static_cast<Real>(k) / (N - 1));
      }
    }
  }
  return pts;
}

// Linear MLS on supports tuned for a 7x7 source grid on [0,1]^2.
constexpr pcms::MLSOptions Linear{.degree = 1};
const pcms::NeighbourhoodOptions GridSupports{.radius = 0.35,
                                              .min_supports = 6};
const pcms::NeighbourhoodOptions GridSupports3D{.radius = 0.8,
                                                .min_supports = 10};

const char* BasisName(pcms::RadialBasisFunction basis)
{
  switch (basis) {
    case pcms::RadialBasisFunction::RBF_GAUSSIAN: return "RBF_GAUSSIAN";
    case pcms::RadialBasisFunction::RBF_C4: return "RBF_C4";
    case pcms::RadialBasisFunction::RBF_CONST: return "RBF_CONST";
    case pcms::RadialBasisFunction::RBF_MULTIQUADRIC: return "RBF_MULTIQUADRIC";
    default: return "UNEXPECTED_BASIS";
  }
}

template <typename Func>
void CheckPolynomialReproduction(unsigned degree,
                                 pcms::RadialBasisFunction basis, Func func,
                                 double abs_tol)
{
  CAPTURE(BasisName(basis));
  CAPTURE(degree);

  auto src = MakeGrid2D(9);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 81, 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::Cartesian::Deferred(),
                                            coords_view),
      pcms::NeighbourhoodOptions{.radius = 0.45,
                                 .min_supports = degree < 2u ? 8 : 16}),
    pcms::MLSOptions{.degree = degree, .basis = basis});
  auto field = fs->CreateFunction<Real>();
  pcms::test::SetField(field, func);

  auto pts = pcms::test::StandardEvalCoords2D();
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));
  pcms::test::CheckEvaluation(*evaluator, field, pts, func, abs_tol);
}

} // namespace

// ============================================================================
// Polynomial reproduction sweep
// ============================================================================

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: reproduces "
          "representative polynomials "
          "across degree and basis options")
{
  auto bases = std::array{pcms::RadialBasisFunction::RBF_GAUSSIAN,
                          pcms::RadialBasisFunction::RBF_C4,
                          pcms::RadialBasisFunction::RBF_CONST,
                          pcms::RadialBasisFunction::RBF_MULTIQUADRIC};

  for (auto basis : bases) {
    CheckPolynomialReproduction(
      0, basis, OMEGA_H_LAMBDA(Real, Real) { return 3.14; }, 5e-3);
    CheckPolynomialReproduction(
      1, basis,
      OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
      5e-3);
    CheckPolynomialReproduction(
      2, basis,
      OMEGA_H_LAMBDA(Real x, Real y) { return x * x + x * y + 2.0 * y * y; },
      5e-3);
  }
}

// ============================================================================
// PointEvaluator reused across two FieldData objects
// ============================================================================

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: same PointEvaluator "
          "reused for two "
          "FieldData objects")
{
  auto src = MakeGrid2D(7);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 49, 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::Cartesian::Deferred(),
                                            coords_view),
      GridSupports),
    Linear);
  auto field_a = fs->CreateFunction<Real>();
  auto field_b = fs->CreateFunction<Real>();

  pcms::test::SetField(
    field_a,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  const Real cval = 7.0;
  pcms::test::SetField(field_b, OMEGA_H_LAMBDA(Real, Real) { return cval; });

  auto pts = pcms::test::StandardEvalCoords2D();
  // Create the PointEvaluator once and reuse it for both fields.
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));

  pcms::test::CheckEvaluation(
    *evaluator, field_a, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
    5e-3);
  pcms::test::CheckEvaluation(
    *evaluator, field_b, pts, OMEGA_H_LAMBDA(Real, Real) { return cval; },
    5e-3);
}

// ============================================================================
// Non-scalar output view must throw
// ============================================================================

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: Evaluate throws for "
          "num_components != 1")
{
  auto src = MakeGrid2D(5);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 25, 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view)));
  auto field = fs->CreateFunction<Real>();

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));

  // Two-component output — must throw
  Kokkos::View<Real**, DeviceMemorySpace> out_device("out", n, 2);
  REQUIRE_THROWS(
    evaluator->Evaluate(field, pcms::test::TagLike(field, out_device)));
}

// ============================================================================
// Default MLSOptions (no explicit options) — smoke test
// ============================================================================

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: default MLSOptions — "
          "smoke evaluation")
{
  auto src = MakeGrid2D(7);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 49, 2);

  // Use default options — no third argument
  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view)));
  auto field = fs->CreateFunction<Real>();
  pcms::test::SetField(field, OMEGA_H_LAMBDA(Real, Real) { return Real(1.0); });

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));

  Kokkos::View<Real**, DeviceMemorySpace> out_device("out", n, 1);
  // Just verify it runs without error and returns finite values
  REQUIRE_NOTHROW(
    evaluator->Evaluate(field, pcms::test::TagLike(field, out_device)));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);
  for (int i = 0; i < n; ++i)
    REQUIRE(std::isfinite(out_host(i, 0)));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: CreatePointEvaluator "
          "rejects a query system that is not the space's")
{
  auto src = MakeGrid2D(5);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 25, 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::CylindricalRZ::Create(),
                                            coords_view),
      GridSupports),
    Linear);

  auto pts = pcms::test::StandardEvalCoords2D();
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  // The hoisted domain-coordinate system identity check fires before any MLS
  // logic.
  REQUIRE_THROWS(fs->CreatePointEvaluator<Real>(
    fs->GetGeometry()->Bind(device_coords.coordinate_view)));
}

// ensure that MLS reconstruction works with CylindricalRZ which has identity
// metric
TEST_CASE("PolynomialReconstructionFunctionSpace MLS: accepts the (R, Z) "
          "poloidal section, whose metric is the identity")
{
  auto src = MakeGrid2D(7);
  Rank2View<Real, HostMemorySpace> coords_view(
    src.data(), static_cast<pcms::LO>(src.size() / 2), 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::CylindricalRZ::Create(),
                                            coords_view),
      GridSupports),
    Linear);

  // A linear field in (R, Z), which a degree-1 MLS must reproduce.
  auto linear_rz = [](Real r, Real z) { return 2.0 + 3.0 * r - 1.5 * z; };
  auto field = fs->CreateFunction<Real>();
  std::vector<Real> dof_values;
  dof_values.reserve(src.size() / 2);
  for (size_t i = 0; i < src.size() / 2; ++i) {
    dof_values.push_back(linear_rz(src[2 * i], src[2 * i + 1]));
  }
  field.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
      pcms::ValueBasis{},
      Rank2View<const Real, HostMemorySpace>(
        dof_values.data(), static_cast<pcms::LO>(dof_values.size()), 1)));

  std::vector<Real> pts{0.3, 0.4, 0.6, 0.7};
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));

  Kokkos::View<Real**, DeviceMemorySpace> out_device("out", 2, 1);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  REQUIRE(out_host(0, 0) == Catch::Approx(linear_rz(0.3, 0.4)).margin(5e-3));
  REQUIRE(out_host(1, 0) == Catch::Approx(linear_rz(0.6, 0.7)).margin(5e-3));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: rejects coordinates "
          "whose metric is not the identity")
{
  auto src = MakeGrid3D(3);
  Rank2View<Real, HostMemorySpace> coords_view(
    src.data(), static_cast<pcms::LO>(src.size() / 3), 3);

  // h_theta = r, so summing squared coordinate differences is not a distance.
  // The geometry checks its sites' system once, when it is created.
  REQUIRE_THROWS_WITH(
    pcms::CreatePolynomialReconstructionSpace(
      pcms::CreatePointCloudGeometry(
        pcms::CoordinateView<HostMemorySpace>(
          pcms::csys::CylindricalRThetaZ::Create(), coords_view),
        GridSupports3D),
      Linear),
    Catch::Matchers::ContainsSubstring("whose metric is the identity"));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: spaces of two degrees "
          "on one geometry share a bind")
{
  auto src = MakeGrid2D(7);
  Rank2View<Real, HostMemorySpace> coords_view(
    src.data(), static_cast<pcms::LO>(src.size() / 2), 2);
  auto first = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::Cartesian::Create(2),
                                            coords_view),
      GridSupports),
    Linear);
  auto geometry = std::dynamic_pointer_cast<const pcms::PointCloudGeometry>(
    first->GetGeometry());
  REQUIRE(geometry != nullptr);
  auto linear = pcms::CreatePolynomialReconstructionSpace(geometry, Linear);
  auto quadratic = pcms::CreatePolynomialReconstructionSpace(
    geometry, pcms::MLSOptions{.degree = 2});
  REQUIRE(linear->GetGeometry() == quadratic->GetGeometry());

  // a quadratic field: degree 2 reproduces it, degree 1 does not
  auto quad = [](Real x, Real y) {
    return 1.0 + x - 2.0 * y + 3.0 * x * x - x * y + 0.5 * y * y;
  };
  std::vector<Real> dof_values;
  for (size_t i = 0; i < src.size() / 2; ++i) {
    dof_values.push_back(quad(src[2 * i], src[2 * i + 1]));
  }
  const pcms::ValueView<const Real, HostMemorySpace> dofs(
    pcms::ValueBasis{},
    Rank2View<const Real, HostMemorySpace>(
      dof_values.data(), static_cast<pcms::LO>(dof_values.size()), 1));
  auto f_linear = linear->CreateFunction<Real>();
  f_linear.SetDOFHolderDataHost(dofs);
  auto f_quadratic = quadratic->CreateFunction<Real>();
  f_quadratic.SetDOFHolderDataHost(dofs);

  std::vector<Real> pts{0.31, 0.42, 0.55, 0.68, 0.47, 0.23};
  const auto bound = pcms::test::Bind(*geometry, pts);
  const auto quadratic_values =
    pcms::test::EvaluateToHost(*quadratic->CreatePointEvaluator<Real>(bound),
                               f_quadratic, pcms::ValueBasis{});
  const auto linear_values = pcms::test::EvaluateToHost(
    *linear->CreatePointEvaluator<Real>(bound), f_linear, pcms::ValueBasis{});
  // the first space accepts it too
  REQUIRE_NOTHROW(first->CreatePointEvaluator<Real>(bound));
  Real linear_error = 0;
  for (size_t i = 0; i < 3; ++i) {
    const Real exact = quad(pts[2 * i], pts[2 * i + 1]);
    CAPTURE(i, exact, quadratic_values(i, 0), linear_values(i, 0));
    REQUIRE(quadratic_values(i, 0) == Catch::Approx(exact).margin(1e-8));
    linear_error =
      std::max(linear_error, std::abs(linear_values(i, 0) - exact));
  }
  REQUIRE(linear_error > 1e-4);
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: radius option is "
          "interpreted as a physical cutoff")
{
  std::vector<Real> src{0.0, 0.0, 0.6, 0.0};
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 2, 2);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::Cartesian::Deferred(),
                                            coords_view),
      pcms::NeighbourhoodOptions{
        .radius = 0.5, .min_supports = 1, .adapt_radius = false}),
    pcms::MLSOptions{.degree = 0,
                     .basis = pcms::RadialBasisFunction::RBF_CONST});
  auto field = fs->CreateFunction<Real>();
  std::vector<Real> dof_values{1.0, 5.0};
  field.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
      pcms::ValueBasis{},
      Rank2View<const Real, HostMemorySpace>(
        dof_values.data(), static_cast<pcms::LO>(dof_values.size()), 1)));

  std::vector<Real> pts{0.0, 0.0};
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));
  Kokkos::View<Real**, DeviceMemorySpace> out_device("out", 1, 1);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  REQUIRE(out_host(0, 0) == Catch::Approx(1.0).margin(1e-8));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: 3D point clouds preserve "
          "z coordinates in "
          "layout and evaluation")
{
  auto src = MakeGrid3D(3);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 27, 3);

  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      pcms::CoordinateView<HostMemorySpace>(pcms::csys::Cartesian::Deferred(),
                                            coords_view),
      GridSupports3D),
    Linear);
  auto layout_coords_device =
    fs->GetLayout()->GetDOFHolderCoordinates().GetValues();
  auto layout_coords = pcms::test::CopyCoordinatesToHost(layout_coords_device);
  REQUIRE(static_cast<int>(layout_coords.extent(1)) == 3);
  for (int i = 0; i < 27; ++i) {
    REQUIRE(layout_coords(i, 0) == Catch::Approx(src[3 * i + 0]));
    REQUIRE(layout_coords(i, 1) == Catch::Approx(src[3 * i + 1]));
    REQUIRE(layout_coords(i, 2) == Catch::Approx(src[3 * i + 2]));
  }

  auto field = fs->CreateFunction<Real>();
  pcms::test::SetField(
    field,
    KOKKOS_LAMBDA(Real x, Real y, Real z) { return x + 2.0 * y + 3.0 * z; });

  std::vector<Real> pts{0.5, 0.5, 0.25, 0.5, 0.5, 0.75};
  auto evaluator =
    fs->CreatePointEvaluator<Real>(pcms::test::Bind(*fs->GetGeometry(), pts));
  Kokkos::View<Real**, DeviceMemorySpace> out_device("out", 2, 1);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  // This 3D case is a small, low-resolution support cloud with MLS weights
  // built from a finite-radius neighborhood rather than exact nodal lookup, so
  // it is checked with a slightly looser tolerance than the denser 2D tests.
  REQUIRE(out_host(0, 0) == Catch::Approx(2.25).margin(1e-2));
  REQUIRE(out_host(1, 0) == Catch::Approx(3.75).margin(1e-2));
  REQUIRE(out_host(1, 0) - out_host(0, 0) == Catch::Approx(1.5).margin(1e-2));
}
