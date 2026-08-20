#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/basis_transformation.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/transfer/interpolator.h>
#include <pcms/transfer/transformed_transfer_operator.hpp>
#include "field_test_utils.h"
#include <cmath>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
namespace values = pcms::values;
using pcms::ComponentScaling;
using pcms::HostMemorySpace;
using pcms::LagrangeFunctionSpace;
using pcms::Real;
using pcms::TransformedTransferOperator;
using pcms::ValueBasis;
using pcms::ValueView;
using pcms::Variance;
using pcms::VarianceSignature;

namespace
{

constexpr double kTol = 1e-10;

// A 3D simplex mesh whose coordinates are interpreted as (r, theta, z).
Omega_h::Mesh BuildCylindricalMesh(Omega_h::Library& lib)
{
  return Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1, 8, 8,
                            6, false);
}

std::shared_ptr<LagrangeFunctionSpace> BuildCylindricalSpace(
  Omega_h::Mesh& mesh, int num_components)
{
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::CylindricalRThetaZ::Create(), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

} // namespace

TEST_CASE("TransformedTransferOperator: same-system basis-only rotation")
{
  auto lib = Omega_h::Library{};
  // One coordinate system, two meshes: source stores CARTESIAN components on a
  // cylindrical-system space (borrowed-basis storage); target declares native
  // cylindrical components.
  auto src_mesh = BuildCylindricalMesh(lib);
  auto src_space = BuildCylindricalSpace(src_mesh, 3);
  auto tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.9, 1.5,
                                     1.0, 5, 5, 5, false);
  auto tgt_space = BuildCylindricalSpace(tgt_mesh, 3);

  auto src = src_space->CreateFunction<Real>("b", values::Vector,
                                             pcms::csys::Cartesian::Create(3));
  auto tgt = tgt_space->CreateFunction<Real>("b", values::Vector,
                                             ComponentScaling::Physical);
  // Constant Cartesian x-unit vector; in cylindrical components at (r,θ,z)
  // it must become (cos θ, -sin θ, 0).
  pcms::test::SetFieldComponents(src, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  TransformedTransferOperator<Real> op(
    std::in_place_type<pcms::Interpolator<Real>>, *src_space, *tgt_space,
    std::make_shared<pcms::CartesianToCylindricalBasis>());
  op.Apply(src, tgt);

  const auto tgt_coords = tgt_space->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(tgt_coords.GetValues().extent(0));
  const auto coords_host =
    pcms::test::CopyCoordinatesToHost(tgt_coords.GetValues());
  auto result = tgt.GetDOFHolderDataHost();
  for (int i = 0; i < n; ++i) {
    const Real theta = coords_host(i, 1);
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(std::cos(theta), kTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(-std::sin(theta), kTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(0.0, kTol));
  }
}

TEST_CASE("TransformedTransferOperator: matching bases take the inner "
          "committed path")
{
  auto lib = Omega_h::Library{};
  auto src_mesh = BuildCylindricalMesh(lib);
  auto src_space = BuildCylindricalSpace(src_mesh, 1);
  auto tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.9, 1.5,
                                     1.0, 5, 5, 5, false);
  auto tgt_space = BuildCylindricalSpace(tgt_mesh, 1);

  auto src = src_space->CreateFunction<Real>("f");
  auto tgt = tgt_space->CreateFunction<Real>("f");
  pcms::test::SetField(
    src, OMEGA_H_LAMBDA(Real r, Real theta, Real z) {
      return 2.0 * r + 3.0 * theta - z;
    });

  // Scalars never need a value transformation; the wrapper delegates straight
  // to the inner committed Apply.
  TransformedTransferOperator<Real> op(
    std::in_place_type<pcms::Interpolator<Real>>, *src_space, *tgt_space,
    std::make_shared<pcms::CartesianToCylindricalBasis>());
  op.Apply(src, tgt);

  const auto tgt_coords = tgt_space->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(tgt_coords.GetValues().extent(0));
  const auto coords_host =
    pcms::test::CopyCoordinatesToHost(tgt_coords.GetValues());
  auto result = tgt.GetDOFHolderDataHost();
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0),
                 WithinAbs(2.0 * coords_host(i, 0) + 3.0 * coords_host(i, 1) -
                             coords_host(i, 2),
                           kTol));
  }
}

TEST_CASE("TransformedTransferOperator: construction and gate errors")
{
  auto lib = Omega_h::Library{};
  auto src_mesh = BuildCylindricalMesh(lib);
  auto src_space = BuildCylindricalSpace(src_mesh, 3);

  SECTION("a null transformation is rejected")
  {
    REQUIRE_THROWS_WITH(TransformedTransferOperator<Real>(
                          std::in_place_type<pcms::Interpolator<Real>>,
                          *src_space, *src_space,
                          std::shared_ptr<const pcms::BasisTransformation>{}),
                        ContainsSubstring("must not be null"));
  }

  SECTION("a transformation that does not bridge the bases throws at Apply")
  {
    auto tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.9, 1.5,
                                       1.0, 4, 4, 4, false);
    auto tgt_space = BuildCylindricalSpace(tgt_mesh, 3);
    // The source stores native cylindrical components, but the transformation
    // starts from borrowed Cartesian ones.
    auto src = src_space->CreateFunction<Real>("b", values::Vector,
                                               ComponentScaling::Physical);
    auto tgt = tgt_space->CreateFunction<Real>(
      "b", values::Vector, pcms::csys::Cartesian::Create(3));
    TransformedTransferOperator<Real> op(
      std::in_place_type<pcms::Interpolator<Real>>, *src_space, *tgt_space,
      std::make_shared<pcms::CartesianToCylindricalBasis>());
    REQUIRE_THROWS_WITH(op.Apply(src, tgt),
                        ContainsSubstring("does not bridge"));
  }
}

// The supported cross-system route while no operator composes it: map the
// query points once, build the native evaluator on the mapped points, and
// rotate the evaluated components with a transformation bound to the same
// points.
TEST_CASE("manual composition: map + native evaluator + basis transformation")
{
  auto lib = Omega_h::Library{};
  auto mesh = BuildCylindricalMesh(lib);
  auto space = BuildCylindricalSpace(mesh, 3);
  auto b = space->CreateFunction<Real>("b", values::Vector,
                                       ComponentScaling::Physical);
  // Constant radial unit field in cylindrical components.
  pcms::test::SetFieldComponents(b, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  const std::vector<Real> cyl_pts = {0.7, 0.3, 0.4, 1.0, 0.9, 0.9};
  std::vector<Real> cart_pts(cyl_pts.size());
  for (size_t i = 0; i < cyl_pts.size(); i += 3) {
    cart_pts[i] = cyl_pts[i] * std::cos(cyl_pts[i + 1]);
    cart_pts[i + 1] = cyl_pts[i] * std::sin(cyl_pts[i + 1]);
    cart_pts[i + 2] = cyl_pts[i + 2];
  }
  const int n = static_cast<int>(cyl_pts.size()) / 3;

  auto query = pcms::test::CreateDeviceCoordinateView(
    cart_pts, pcms::csys::Cartesian::Create(3), 3);
  const auto mapped = pcms::CartesianToCylindrical{}.Map(query.coordinate_view);
  auto evaluator = space->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(mapped.View()));
  const auto law = pcms::CylindricalToCartesianBasis{}.Bind(mapped.View());

  Kokkos::View<Real**, DeviceMemorySpace> native("native", n, 3);
  evaluator->Evaluate(b, pcms::MakeRank2View(native));
  Kokkos::View<Real**, DeviceMemorySpace> rotated("rotated", n, 3);
  law->Apply(
    ValueView<const Real, DeviceMemorySpace>(
      ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                 ComponentScaling::Physical,
                 VarianceSignature{Variance::Contravariant}},
      pcms::MakeConstRank2View(native)),
    ValueView<Real, DeviceMemorySpace>(
      ValueBasis{pcms::csys::Cartesian::Create(3), ComponentScaling::Physical,
                 VarianceSignature{Variance::Contravariant}},
      pcms::MakeRank2View(rotated)));

  Kokkos::View<Real**, HostMemorySpace> host("host", n, 3);
  pcms::DeepCopyMismatchLayouts(host, rotated);
  for (int i = 0; i < n; ++i) {
    const Real theta = cyl_pts[3 * i + 1];
    CAPTURE(i);
    REQUIRE_THAT(host(i, 0), WithinAbs(std::cos(theta), kTol));
    REQUIRE_THAT(host(i, 1), WithinAbs(std::sin(theta), kTol));
  }
}

TEST_CASE("tagged writes gate on the field's declaration")
{
  auto lib = Omega_h::Library{};
  auto mesh = BuildCylindricalMesh(lib);
  auto space = BuildCylindricalSpace(mesh, 3);
  auto b = space->CreateFunction<Real>("b", values::Vector,
                                       ComponentScaling::Physical);
  const int n = static_cast<int>(
    space->GetLayout()->GetDOFHolderCoordinates().GetValues().extent(0));
  std::vector<Real> data(static_cast<size_t>(n) * 3, 1.0);
  const pcms::Rank2View<const Real, HostMemorySpace> raw(data.data(), n, 3);

  SECTION("a matching tag is accepted")
  {
    b.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
      ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                 ComponentScaling::Physical,
                 VarianceSignature{Variance::Contravariant}},
      raw));
    REQUIRE_THAT(b.GetDOFHolderDataHost()(0, 0),
                 WithinAbs(1.0, 1e-14)); // tagged getter forwards indexing
  }

  SECTION("a mismatched basis claim is rejected")
  {
    REQUIRE_THROWS_WITH(
      b.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
        ValueBasis{pcms::csys::Cartesian::Create(3), ComponentScaling::Physical,
                   VarianceSignature{Variance::Contravariant}},
        raw)),
      ContainsSubstring("does not match this field's declaration"));
  }

  SECTION("a mismatched value type is rejected")
  {
    REQUIRE_THROWS_WITH(
      b.SetDOFHolderDataHost(
        ValueView<const Real, HostMemorySpace>(ValueBasis{}, raw)),
      ContainsSubstring("does not match this field's declaration"));
  }

  SECTION("the unchecked path inherits the declaration unverified")
  {
    REQUIRE_NOTHROW(b.SetDOFHolderDataUncheckedHost(raw));
  }
}
