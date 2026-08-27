#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/basis_transformation.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/transfer/interpolator.h>
#include <pcms/transfer/transformed_transfer_operator.hpp>
#include "field_test_utils.h"
#include <memory>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

using Catch::Matchers::ContainsSubstring;
using pcms::DeviceMemorySpace;
using pcms::Field;
using pcms::LagrangeFunctionSpace;
using pcms::Real;
using pcms::TransferKey;
using pcms::TransformedTransferOperator;
using pcms::ValueView;
namespace values = pcms::values;
using pcms::ComponentScaling;

namespace
{

// The keyed Apply is passkey-protected; only a TransferOperator can mint the
// key, so the tests borrow one through a derived type.
class KeyMaker : public pcms::TransferOperator<Real>
{
public:
  static TransferKey Key() { return MakeTransferKey(); }

  void Apply(const Field<Real>&, Field<Real>&) const override {}
  void Apply(TransferKey, const Field<Real>&,
             ValueView<Real, DeviceMemorySpace>) const override
  {
  }
};

// A 3D simplex mesh whose coordinates are interpreted as (r, theta, z).
Omega_h::Mesh BuildCylindricalMesh(Omega_h::Library& lib, int divisions)
{
  return Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1,
                            divisions, divisions, divisions, false);
}

std::shared_ptr<LagrangeFunctionSpace> BuildCylindricalSpace(
  Omega_h::Mesh& mesh, int num_components)
{
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::CylindricalRThetaZ::Create(), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

int NumDOFHolders(const LagrangeFunctionSpace& space)
{
  return static_cast<int>(
    space.GetLayout()->GetDOFHolderCoordinates().GetValues().extent(0));
}

} // namespace

TEST_CASE("PointEvaluator::Evaluate gates on the source field's stored basis")
{
  auto lib = Omega_h::Library{};
  auto mesh = BuildCylindricalMesh(lib, 6);
  auto space = BuildCylindricalSpace(mesh, 3);
  // Borrowed-basis storage: Cartesian components on a cylindrical-system space.
  auto cartesian = space->CreateFunction<Real>(
    "b", values::Vector, pcms::csys::Cartesian::Create(3));
  auto cylindrical = space->CreateFunction<Real>("b_native", values::Vector,
                                                 ComponentScaling::Physical);
  pcms::test::SetFieldComponents(cartesian, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  const std::vector<Real> pts = {0.7, 0.3, 0.4, 1.0, 0.9, 0.9};
  auto query = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
  auto evaluator = space->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(query.coordinate_view));
  const int n = static_cast<int>(pts.size()) / 3;
  Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 3);

  SECTION("a matching tag is accepted")
  {
    REQUIRE_NOTHROW(evaluator->Evaluate(
      cartesian,
      ValueView<Real, DeviceMemorySpace>(cartesian.GetData().GetValueBasis(),
                                         pcms::MakeRank2View(out))));
  }

  SECTION("a mismatched basis claim is rejected")
  {
    REQUIRE_THROWS_WITH(
      evaluator->Evaluate(cartesian, ValueView<Real, DeviceMemorySpace>(
                                       cylindrical.GetData().GetValueBasis(),
                                       pcms::MakeRank2View(out))),
      ContainsSubstring("does not match the basis this call writes"));
  }
}

TEST_CASE("Interpolator's keyed Apply gates on the source's stored basis")
{
  auto lib = Omega_h::Library{};
  auto src_mesh = BuildCylindricalMesh(lib, 6);
  auto src_space = BuildCylindricalSpace(src_mesh, 3);
  auto tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.9, 1.5,
                                     1.0, 4, 4, 4, false);
  auto tgt_space = BuildCylindricalSpace(tgt_mesh, 3);

  auto src = src_space->CreateFunction<Real>("b", values::Vector,
                                             ComponentScaling::Physical);
  auto foreign = src_space->CreateFunction<Real>(
    "b_cart", values::Vector, pcms::csys::Cartesian::Create(3));
  pcms::test::SetFieldComponents(src, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  pcms::Interpolator<Real> interp(*src_space, *tgt_space);
  Kokkos::View<Real**, DeviceMemorySpace> out("out", NumDOFHolders(*tgt_space),
                                              3);

  SECTION("a matching tag is accepted")
  {
    REQUIRE_NOTHROW(
      interp.Apply(KeyMaker::Key(), src,
                   ValueView<Real, DeviceMemorySpace>(
                     src.GetData().GetValueBasis(), pcms::MakeRank2View(out))));
  }

  SECTION("a mismatched basis claim is rejected")
  {
    REQUIRE_THROWS_WITH(
      interp.Apply(
        KeyMaker::Key(), src,
        ValueView<Real, DeviceMemorySpace>(foreign.GetData().GetValueBasis(),
                                           pcms::MakeRank2View(out))),
      ContainsSubstring("does not match the basis this call writes"));
  }
}

TEST_CASE("TransformedTransferOperator's keyed Apply writes the "
          "transformation's target basis")
{
  auto lib = Omega_h::Library{};
  auto src_mesh = BuildCylindricalMesh(lib, 6);
  auto src_space = BuildCylindricalSpace(src_mesh, 3);
  auto tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.9, 1.5,
                                     1.0, 4, 4, 4, false);
  auto tgt_space = BuildCylindricalSpace(tgt_mesh, 3);

  // Source stores borrowed Cartesian components; the transformation carries
  // them to the native cylindrical components the target declares.
  auto src = src_space->CreateFunction<Real>("b", values::Vector,
                                             pcms::csys::Cartesian::Create(3));
  auto tgt = tgt_space->CreateFunction<Real>("b", values::Vector,
                                             ComponentScaling::Physical);
  pcms::test::SetFieldComponents(src, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  TransformedTransferOperator<Real> op(
    std::in_place_type<pcms::Interpolator<Real>>, *src_space, *tgt_space,
    std::make_shared<pcms::CartesianToCylindricalBasis>());
  Kokkos::View<Real**, DeviceMemorySpace> out("out", NumDOFHolders(*tgt_space),
                                              3);

  SECTION("the transformation's target basis is accepted")
  {
    REQUIRE_NOTHROW(
      op.Apply(KeyMaker::Key(), src,
               ValueView<Real, DeviceMemorySpace>(tgt.GetData().GetValueBasis(),
                                                  pcms::MakeRank2View(out))));
  }

  SECTION("the source's stored basis is rejected")
  {
    REQUIRE_THROWS_WITH(
      op.Apply(KeyMaker::Key(), src,
               ValueView<Real, DeviceMemorySpace>(src.GetData().GetValueBasis(),
                                                  pcms::MakeRank2View(out))),
      ContainsSubstring("does not match the basis this call writes"));
  }
}
