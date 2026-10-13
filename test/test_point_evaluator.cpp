#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
#include <stdexcept>

#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/function_space/spline.h"
#include "pcms/field/field_data.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/value_view.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/utility/arrays.h"
#include "field_test_utils.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"

using pcms::CoordinateSystem;
using pcms::Real;

namespace
{

// Evaluate checks the output's shape and tag before evaluating, the same way
// for every evaluator. `field` must be scalar.
void CheckOutputValidation(const pcms::PointEvaluator<Real>& evaluator,
                           const pcms::Field<Real>& field, int num_points)
{
  using OutView = Kokkos::View<Real**, pcms::DeviceMemorySpace>;
  REQUIRE(evaluator.NumPoints() == num_points);
  REQUIRE(field.GetLayout().GetNumComponents() == 1);
  auto scalar = [](OutView v) {
    return pcms::ValueView<Real, pcms::DeviceMemorySpace>(
      pcms::ValueBasis{}, pcms::MakeRank2View(v));
  };
  const auto n = static_cast<size_t>(num_points);

  REQUIRE_NOTHROW(evaluator.Evaluate(field, scalar(OutView("ok", n, 1))));
  REQUIRE_THROWS_AS(
    evaluator.Evaluate(field, scalar(OutView("short", n - 1, 1))),
    std::invalid_argument);
  REQUIRE_THROWS_AS(
    evaluator.Evaluate(field, scalar(OutView("long", n + 1, 1))),
    std::invalid_argument);
  REQUIRE_THROWS_AS(evaluator.Evaluate(field, scalar(OutView("wide", n, 2))),
                    std::invalid_argument);
  // One component is also a valid 1D vector, so only the tag is wrong.
  OutView out("tagged", n, 1);
  REQUIRE_THROWS_AS(
    evaluator.Evaluate(field,
                       pcms::ValueView<Real, pcms::DeviceMemorySpace>(
                         pcms::ValueBasis{pcms::csys::Cartesian::Create(1),
                                          pcms::ComponentScaling::Physical,
                                          pcms::values::Vector},
                         pcms::MakeRank2View(out))),
    pcms::pcms_error);
}

// Two points inside [0,1]^2 and two outside.
const std::vector<Real> kMixedPoints = {0.25, 0.25, -0.5, 0.5,
                                        0.7,  0.3,  1.5,  0.5};

std::unique_ptr<pcms::PointEvaluator<Real>> CreateFillEvaluator(
  const pcms::FunctionSpace& space)
{
  auto coords = pcms::test::CreateDeviceCoordinateView(
    kMixedPoints, pcms::csys::Cartesian::Deferred());
  return space.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(
      coords.coordinate_view,
      pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -1.0}));
}

} // namespace

TEST_CASE("PointEvaluator::Evaluate checks the output for every evaluator")
{
  auto lib = Omega_h::Library{};
  const int n = static_cast<int>(kMixedPoints.size()) / 2;
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {10, 10};

  SECTION("Omega_h Lagrange")
  {
    auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10,
                                   10, 0, false);
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    CheckOutputValidation(*CreateFillEvaluator(*space),
                          space->CreateFunction<Real>(), n);
  }
#ifdef PCMS_ENABLE_MESHFIELDS
  SECTION("MeshFields Lagrange")
  {
    auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10,
                                   10, 0, false);
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
      pcms::LagrangeFunctionSpace::Backend::MeshFields);
    CheckOutputValidation(*CreateFillEvaluator(*space),
                          space->CreateFunction<Real>(), n);
  }
#endif
  SECTION("uniform grid Lagrange")
  {
    auto space = pcms::LagrangeFunctionSpace::FromUniformGrid(
      grid, 1, pcms::csys::Cartesian::Deferred(), 1);
    CheckOutputValidation(*CreateFillEvaluator(*space),
                          space->CreateFunction<Real>(), n);
  }
  SECTION("uniform grid spline")
  {
    auto space = pcms::SplineFunctionSpace::FromUniformGrid(
      grid, pcms::csys::Cartesian::Deferred());
    CheckOutputValidation(*CreateFillEvaluator(*space),
                          space->CreateFunction<Real>(), n);
  }
  SECTION("MLS point cloud")
  {
    std::vector<Real> sources;
    for (int j = 0; j < 5; ++j) {
      for (int i = 0; i < 5; ++i) {
        sources.push_back(0.25 * i);
        sources.push_back(0.25 * j);
      }
    }
    auto space = pcms::PolynomialReconstructionFunctionSpace::Create(
      pcms::CoordinateView<pcms::HostMemorySpace>(
        pcms::csys::Cartesian::Deferred(),
        pcms::Rank2View<const Real, pcms::HostMemorySpace>(sources.data(), 25,
                                                           2)));
    const auto pts = pcms::test::StandardEvalCoords2D();
    auto coords = pcms::test::CreateDeviceCoordinateView(
      pts, pcms::csys::Cartesian::Deferred());
    auto evaluator = space->CreatePointEvaluator<Real>(
      pcms::EvaluationRequest::FromCoordinates(coords.coordinate_view));
    CheckOutputValidation(*evaluator, space->CreateFunction<Real>(),
                          static_cast<int>(pts.size()) / 2);
  }
}

// ============================================================================
// OmegaH order-1 — basic evaluation via new API
// ============================================================================

TEST_CASE("PointEvaluator: OmegaH order-1 linear evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  pcms::test::CheckEvaluation(
    *evaluator, field_data, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
}

// ============================================================================
// OmegaH order-1 — repeated evaluation: same PointEvaluator, two FieldDatas
// ============================================================================

TEST_CASE("PointEvaluator: same evaluator reused for two FieldData objects")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 50, 50, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  auto field_a = factory->CreateFunction<Real>();
  auto field_b = factory->CreateFunction<Real>();

  // field_a: linear_f;  field_b: constant 42
  pcms::test::SetField(
    field_a.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  pcms::test::SetField(
    field_b.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real, Real) { return Real(42); });

  auto pts = pcms::test::StandardEvalCoords2D();
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());

  // Create the PointEvaluator once and reuse it for both fields.
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  pcms::test::CheckEvaluation(
    *evaluator, field_a, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  pcms::test::CheckEvaluation(
    *evaluator, field_b, pts, OMEGA_H_LAMBDA(Real, Real) { return Real(42); });
}

// ============================================================================
// OmegaH order-1 — OutOfBoundsPolicy::FILL
// ============================================================================

TEST_CASE("PointEvaluator: OmegaH order-1 out-of-bounds fill")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 20, 20, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  // Points clearly outside [0,1]^2
  const auto outside_pts = pcms::test::StandardOutsideCoords2D();
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    outside_pts, pcms::csys::Cartesian::Deferred());
  pcms::OutOfBoundsPolicy policy{pcms::OutOfBoundsMode::FILL, -999.0};
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view,
                                             policy));
  pcms::test::CheckFillMode(*evaluator, field_data, -999.0, outside_pts);
}

// ============================================================================
// UniformGrid — basic evaluation via new API
// ============================================================================

TEST_CASE("PointEvaluator: UniformGrid order-1 linear evaluation")
{
  // 2D grid: [0,1]^2 with 10x10 divisions
  const int N = 10;
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {N, N};

  auto factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::csys::Cartesian::Deferred(), 1);

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  pcms::test::CheckEvaluation(
    *evaluator, field_data, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
    1e-8);
}

TEST_CASE("PointEvaluator: SplineFunctionSpace uniform-grid evaluation")
{
  const int N = 10;
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {N, N};

  auto factory = pcms::SplineFunctionSpace::FromUniformGrid(
    grid, pcms::csys::Cartesian::Deferred());

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  pcms::test::CheckEvaluation(
    *evaluator, field_data, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
    1e-8);
}

TEST_CASE("PointEvaluator: uniform grids report the points they fill")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {10, 10};
  const std::vector<Real> pts = {0.25, 0.25, -0.5, 0.5, 0.7, 0.3, 1.5, 0.5};
  const std::vector<bool> is_inside = {true, false, true, false};
  auto linear = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return pcms::test::linear_f(x, y);
  };

  SECTION("Lagrange")
  {
    auto factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
      grid, 1, pcms::csys::Cartesian::Deferred(), 1);
    auto field = factory->CreateFunction<Real>();
    pcms::test::SetField(field.GetData(), *factory->GetLayout(), linear);
    pcms::test::CheckEvaluationWithFill(factory, field, pts, is_inside, linear,
                                        -999.0, 1e-8);
  }
  SECTION("spline")
  {
    auto factory = pcms::SplineFunctionSpace::FromUniformGrid(
      grid, pcms::csys::Cartesian::Deferred());
    auto field = factory->CreateFunction<Real>();
    pcms::test::SetField(field.GetData(), *factory->GetLayout(), linear);
    pcms::test::CheckEvaluationWithFill(factory, field, pts, is_inside, linear,
                                        -999.0, 1e-8);
  }
  SECTION("no fill without the FILL policy")
  {
    auto factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
      grid, 1, pcms::csys::Cartesian::Deferred(), 1);
    auto coords = pcms::test::CreateDeviceCoordinateView(
      pts, pcms::csys::Cartesian::Deferred());
    auto evaluator = factory->CreatePointEvaluator<Real>(
      pcms::EvaluationRequest::FromCoordinates(
        coords.coordinate_view,
        pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::NEAREST_BOUNDARY}));
    REQUIRE(evaluator->FilledPoints().extent(0) == 0);
  }
}

// ============================================================================
// FieldLayout — metadata interface
// ============================================================================

TEST_CASE("FieldLayout: metadata queries")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto layout = factory->GetLayout();

  auto coords = layout->GetDOFHolderCoordinates();
  REQUIRE(pcms::SameCoordinateSystem(coords.GetCoordinateSystem(),
                                     pcms::csys::Cartesian::Create(2)));
  REQUIRE(coords.GetValues().extent(0) > 0);
  REQUIRE(coords.GetValues().extent(1) == 2);
}

// ============================================================================
// FieldData / FieldLayout metadata queries
// ============================================================================

TEST_CASE("FieldData: layout metadata queries")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field_data = factory->CreateFunction<Real>();

  auto coords = factory->GetLayout()->GetDOFHolderCoordinates();
  REQUIRE(pcms::SameCoordinateSystem(coords.GetCoordinateSystem(),
                                     pcms::csys::Cartesian::Create(2)));
  REQUIRE(coords.GetValues().extent(0) > 0);
  REQUIRE(coords.GetValues().extent(1) == 2);
}

// ============================================================================
// CreateFieldData / SimpleFieldData round-trip
// ============================================================================

TEST_CASE("SimpleFieldData: set and get DOF holder data round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field_data = factory->CreateFunction<Real>();

  auto& layout = *factory->GetLayout();
  int n = layout.GetNumOwnedDofHolder();
  REQUIRE(n > 0);

  // Write sequential values
  std::vector<Real> data_in(n);
  for (int i = 0; i < n; ++i)
    data_in[i] = static_cast<Real>(i) * 0.5;

  field_data.GetData().SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(data_in.data(), n, 1));

  auto data_out =
    pcms::FlattenToRank1View(field_data.GetData().GetDOFHolderDataHost());
  REQUIRE(static_cast<int>(data_out.size()) == n);
  for (int i = 0; i < n; ++i) {
    REQUIRE(data_out[i] == Catch::Approx(data_in[i]));
  }
}

// ============================================================================
// MeshFields FieldEvaluatorFactory metadata (only when MeshFields is enabled)
// ============================================================================

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("FieldLayout: MeshFields metadata queries")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto layout = factory->GetLayout();
  auto coords = layout->GetDOFHolderCoordinates();
  REQUIRE(pcms::SameCoordinateSystem(coords.GetCoordinateSystem(),
                                     pcms::csys::Cartesian::Create(2)));
  REQUIRE(coords.GetValues().extent(0) > 0);
  REQUIRE(coords.GetValues().extent(1) == 2);
}

TEST_CASE("PointEvaluator: MeshFields order-1 linear evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  pcms::test::CheckEvaluation(
    *evaluator, field_data, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
}

TEST_CASE("PointEvaluator: MeshFields out-of-bounds fill")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 20, 20, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto field_data = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_data.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  const auto outside_pts = pcms::test::StandardOutsideCoords2D();
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    outside_pts, pcms::csys::Cartesian::Deferred());
  pcms::OutOfBoundsPolicy policy{pcms::OutOfBoundsMode::FILL, -999.0};
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view,
                                             policy));
  pcms::test::CheckFillMode(*evaluator, field_data, -999.0, outside_pts);
}

TEST_CASE(
  "PointEvaluator: MeshFields same evaluator reused for two FieldData objects")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 50, 50, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto field_a = factory->CreateFunction<Real>();
  auto field_b = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field_a.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  pcms::test::SetField(
    field_b.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real, Real) { return Real(42); });

  auto pts = pcms::test::StandardEvalCoords2D();
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());

  // Create the PointEvaluator once and reuse it for both fields.
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  pcms::test::CheckEvaluation(
    *evaluator, field_a, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  pcms::test::CheckEvaluation(
    *evaluator, field_b, pts, OMEGA_H_LAMBDA(Real, Real) { return Real(42); });
}

TEST_CASE("PointEvaluator: MeshFields order-1 multi-component (2) evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 50, 50, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 2, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto field = factory->CreateFunction<Real>();
  auto layout = factory->GetLayout();
  REQUIRE(layout->GetNumComponents() == 2);
  int num_dof = layout->GetNumOwnedDofHolder();

  // Component 0: f(x,y)=x+y, Component 1: f(x,y)=2x-y
  auto dof_coords_mdspan = layout->GetDOFHolderCoordinates().GetValues();
  Kokkos::View<Real**, pcms::DeviceMemorySpace> dof_dev("dof_dev", num_dof, 2);
  pcms::ConvertMismatchLayoutView2D(dof_dev, dof_coords_mdspan);
  auto coords_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), dof_dev);

  std::vector<Real> host_data(static_cast<size_t>(num_dof * 2));
  for (int i = 0; i < num_dof; ++i) {
    Real x = coords_host(i, 0);
    Real y = coords_host(i, 1);
    host_data[static_cast<size_t>(i) * 2 + 0] = x + y;
    host_data[static_cast<size_t>(i) * 2 + 1] = 2.0 * x - y;
  }
  field.GetData().SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(host_data.data(),
                                                       num_dof, 2));

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<Real**, pcms::DeviceMemorySpace> out("out", n, 2);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);

  for (int i = 0; i < n; ++i) {
    Real x = pts[2 * static_cast<size_t>(i)];
    Real y = pts[2 * static_cast<size_t>(i) + 1];
    INFO("Point " << i << " (" << x << ", " << y << ")");
    REQUIRE(out_host(i, 0) == Catch::Approx(x + y).margin(1e-8));
    REQUIRE(out_host(i, 1) == Catch::Approx(2.0 * x - y).margin(1e-8));
  }
}

TEST_CASE("PointEvaluator: MeshFields order-1 multi-component (3) evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 50, 50, 0, false);

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 3, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  auto field = factory->CreateFunction<Real>();
  auto layout = factory->GetLayout();
  REQUIRE(layout->GetNumComponents() == 3);
  int num_dof = layout->GetNumOwnedDofHolder();

  // Component 0: x, Component 1: y, Component 2: x*y
  auto dof_coords_mdspan = layout->GetDOFHolderCoordinates().GetValues();
  Kokkos::View<Real**, pcms::DeviceMemorySpace> dof_dev("dof_dev", num_dof, 2);
  pcms::ConvertMismatchLayoutView2D(dof_dev, dof_coords_mdspan);
  auto coords_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), dof_dev);

  std::vector<Real> host_data(static_cast<size_t>(num_dof * 3));
  for (int i = 0; i < num_dof; ++i) {
    Real x = coords_host(i, 0);
    Real y = coords_host(i, 1);
    host_data[static_cast<size_t>(i) * 3 + 0] = x;
    host_data[static_cast<size_t>(i) * 3 + 1] = y;
    host_data[static_cast<size_t>(i) * 3 + 2] = x * y;
  }
  field.GetData().SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(host_data.data(),
                                                       num_dof, 3));

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::Cartesian::Deferred());
  auto evaluator = factory->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<Real**, pcms::DeviceMemorySpace> out("out", n, 3);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);

  for (int i = 0; i < n; ++i) {
    Real x = pts[2 * static_cast<size_t>(i)];
    Real y = pts[2 * static_cast<size_t>(i) + 1];
    INFO("Point " << i << " (" << x << ", " << y << ")");
    REQUIRE(out_host(i, 0) == Catch::Approx(x).margin(1e-8));
    REQUIRE(out_host(i, 1) == Catch::Approx(y).margin(1e-8));
    REQUIRE(out_host(i, 2) == Catch::Approx(x * y).margin(1e-8));
  }
}
#endif // PCMS_ENABLE_MESHFIELDS
