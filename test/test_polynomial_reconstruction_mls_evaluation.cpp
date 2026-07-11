#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/field_metadata.h"
#include "pcms/field/coordinate_system.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "field_test_utils.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

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

// MLSOptions tuned for a 7x7 source grid on [0,1]^2.
pcms::MLSOptions DefaultTestOptions()
{
  pcms::MLSOptions opts;
  opts.radius = 0.35;
  opts.min_req_supports = 6;
  opts.degree = 1;
  opts.adapt_radius = true;
  return opts;
}

pcms::MLSOptions DefaultTestOptions3D()
{
  pcms::MLSOptions opts;
  opts.radius = 0.8;
  opts.min_req_supports = 10;
  opts.degree = 1;
  opts.adapt_radius = true;
  return opts;
}

// Interior query points for a unit box — same as StandardEvalCoords2D.
std::vector<Real> QueryPoints()
{
  return pcms::test::StandardEvalCoords2D();
}

pcms::MLSOptions SweepTestOptions(unsigned degree,
                                  pcms::RadialBasisFunction basis)
{
  pcms::MLSOptions opts;
  opts.radius = 0.45;
  opts.min_req_supports = (degree < 2u) ? 8u : 16u;
  opts.degree = degree;
  opts.adapt_radius = true;
  opts.basis = basis;
  return opts;
}

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

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian, SweepTestOptions(degree, basis));
  auto field = fs.CreateField<Real>(pcms::FieldMetadata{});
  pcms::test::SetField(field.GetData(), *fs.GetLayout(), func);

  auto pts = QueryPoints();
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
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
      1, basis, OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; }, 5e-3);
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

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian, DefaultTestOptions());
  auto field_a = fs.CreateField<Real>(pcms::FieldMetadata{});
  auto field_b = fs.CreateField<Real>(pcms::FieldMetadata{});

  pcms::test::SetField(
    field_a.GetData(), *fs.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
  const Real cval = 7.0;
  pcms::test::SetField(
    field_b.GetData(), *fs.GetLayout(),
    OMEGA_H_LAMBDA(Real, Real) { return cval; });

  auto pts = QueryPoints();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<Real*, DeviceMemorySpace> out_a_device("out_a", n);
  Kokkos::View<Real*, DeviceMemorySpace> out_b_device("out_b", n);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<Real, DeviceMemorySpace, LayoutPolicy> view_a(out_a_device.data(),
                                                          n, 1);
  Rank2View<Real, DeviceMemorySpace, LayoutPolicy> view_b(out_b_device.data(),
                                                          n, 1);

  evaluator->Evaluate(field_a, view_a);
  evaluator->Evaluate(field_b, view_b);

  auto out_a_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_a_device);
  auto out_b_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_b_device);

  for (int i = 0; i < n; ++i) {
    Real x = pts[2 * static_cast<size_t>(i)],
         y = pts[2 * static_cast<size_t>(i) + 1];
    REQUIRE(out_a_host(i) ==
            Catch::Approx(pcms::test::linear_f(x, y)).margin(5e-3));
    REQUIRE(out_b_host(i) == Catch::Approx(cval).margin(5e-3));
  }
}

// ============================================================================
// Non-scalar output view must throw
// ============================================================================

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: Evaluate throws for "
          "num_components != 1")
{
  auto src = MakeGrid2D(5);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 25, 2);

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);
  auto field = fs.CreateField<Real>(pcms::FieldMetadata{});

  auto pts = QueryPoints();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  // Two-component output — must throw
  Kokkos::View<Real*, DeviceMemorySpace> out_device("out",
                                                    static_cast<size_t>(n) * 2);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  auto out_view =
    Rank2View<Real, DeviceMemorySpace, LayoutPolicy>(out_device.data(), n, 2);
  REQUIRE_THROWS(evaluator->Evaluate(field, out_view));
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
  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);
  auto field = fs.CreateField<Real>(pcms::FieldMetadata{});
  pcms::test::SetField(
    field.GetData(), *fs.GetLayout(),
    OMEGA_H_LAMBDA(Real, Real) { return Real(1.0); });

  auto pts = QueryPoints();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<Real*, DeviceMemorySpace> out_device("out", n);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  auto out_view =
    Rank2View<Real, DeviceMemorySpace, LayoutPolicy>(out_device.data(), n, 1);
  // Just verify it runs without error and returns finite values
  REQUIRE_NOTHROW(evaluator->Evaluate(field, out_view));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);
  for (int i = 0; i < n; ++i)
    REQUIRE(std::isfinite(out_host(i)));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: CreatePointEvaluator "
          "rejects coordinate "
          "system mismatch")
{
  auto src = MakeGrid2D(5);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 25, 2);

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cylindrical, DefaultTestOptions());

  auto pts = QueryPoints();
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  REQUIRE_THROWS(fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view)));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: CreatePointEvaluator "
          "rejects non-Cartesian "
          "point-cloud coordinates")
{
  auto src = MakeGrid2D(5);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 25, 2);

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cylindrical, DefaultTestOptions());

  auto pts = QueryPoints();
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cylindrical);
  REQUIRE_THROWS(fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view)));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: radius option is "
          "interpreted as a physical cutoff")
{
  std::vector<Real> src{0.0, 0.0, 0.6, 0.0};
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 2, 2);

  pcms::MLSOptions opts;
  opts.radius = 0.5;
  opts.min_req_supports = 1;
  opts.degree = 0;
  opts.adapt_radius = false;
  opts.basis = pcms::RadialBasisFunction::RBF_CONST;

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian, opts);
  auto field = fs.CreateField<Real>(pcms::FieldMetadata{});
  std::vector<Real> dof_values{1.0, 5.0};
  field.GetData().SetDOFHolderDataHost(Rank2View<const Real, HostMemorySpace>(
    dof_values.data(), static_cast<pcms::LO>(dof_values.size()), 1));

  std::vector<Real> pts{0.0, 0.0};
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  Kokkos::View<Real*, DeviceMemorySpace> out_device("out", 1);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<Real, DeviceMemorySpace, LayoutPolicy> out_view(out_device.data(),
                                                            1, 1);
  evaluator->Evaluate(field, out_view);
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  REQUIRE(out_host(0) == Catch::Approx(1.0).margin(1e-8));
}

TEST_CASE("PolynomialReconstructionFunctionSpace MLS: 3D point clouds preserve "
          "z coordinates in "
          "layout and evaluation")
{
  auto src = MakeGrid3D(3);
  Rank2View<Real, HostMemorySpace> coords_view(src.data(), 27, 3);

  auto fs = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian, DefaultTestOptions3D());
  auto layout_coords_device =
    fs.GetLayout()->GetDOFHolderCoordinates().GetValues();
  auto layout_coords =
    pcms::test::CopyCoordinatesToHost(layout_coords_device, 27, 3);
  REQUIRE(static_cast<int>(layout_coords.extent(1)) == 3);
  for (int i = 0; i < 27; ++i) {
    REQUIRE(layout_coords(i, 0) == Catch::Approx(src[3 * i + 0]));
    REQUIRE(layout_coords(i, 1) == Catch::Approx(src[3 * i + 1]));
    REQUIRE(layout_coords(i, 2) == Catch::Approx(src[3 * i + 2]));
  }

  auto field = fs.CreateField<Real>(pcms::FieldMetadata{});
  std::vector<Real> dof_values(27);
  for (int i = 0; i < 27; ++i) {
    const Real x = src[3 * i + 0];
    const Real y = src[3 * i + 1];
    const Real z = src[3 * i + 2];
    dof_values[i] = x + 2.0 * y + 3.0 * z;
  }
  field.GetData().SetDOFHolderDataHost(Rank2View<const Real, HostMemorySpace>(
    dof_values.data(), static_cast<pcms::LO>(dof_values.size()), 1));

  std::vector<Real> pts{0.5, 0.5, 0.25, 0.5, 0.5, 0.75};
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, CoordinateSystem::Cartesian, 3);
  auto evaluator = fs.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  Kokkos::View<Real*, DeviceMemorySpace> out_device("out", 2);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<Real, DeviceMemorySpace, LayoutPolicy> out_view(out_device.data(),
                                                            2, 1);
  evaluator->Evaluate(field, out_view);
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  // This 3D case is a small, low-resolution support cloud with MLS weights
  // built from a finite-radius neighborhood rather than exact nodal lookup, so
  // it is checked with a slightly looser tolerance than the denser 2D tests.
  REQUIRE(out_host(0) == Catch::Approx(2.25).margin(1e-2));
  REQUIRE(out_host(1) == Catch::Approx(3.75).margin(1e-2));
  REQUIRE(out_host(1) - out_host(0) == Catch::Approx(1.5).margin(1e-2));
}
