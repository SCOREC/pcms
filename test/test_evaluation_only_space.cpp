#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space.h>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/layout/empty.h>
#include <pcms/field/transformed_source.hpp>
#include <pcms/transfer/interpolator.h>
#include "field_test_utils.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::LO;
using pcms::Real;
using pcms::ValueBasis;

// An evaluation-only source, the pattern external codes (NIMROD, Fusion-IO)
// follow: a FunctionSpace on an EmptyFieldLayout whose field data names a
// quantity but provides no DOFAccess. Here the quantity is analytic:
// component c of f at (R, theta, Z) is 2R + 3 theta - Z + c.
namespace
{

KOKKOS_INLINE_FUNCTION Real MockValue(Real r, Real theta, Real z, int c)
{
  return 2.0 * r + 3.0 * theta - z + c;
}

class MockData final : public pcms::FieldData<Real>
{
public:
  pcms::FieldValueType GetValueType() const override
  {
    return pcms::FieldValueType::Scalar;
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

private:
  ValueBasis basis_;
};

class MockEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  explicit MockEvaluator(const pcms::EvaluationRequest& request)
    : points_("mock_points", request.coords.GetValues().extent(0),
              request.coords.GetValues().extent(1))
  {
    pcms::CopyDeviceRank2ViewToRank2View(pcms::MakeRank2View(points_),
                                         request.coords.GetValues());
  }

  void Evaluate(const pcms::Field<Real>& field,
                pcms::ValueView<Real, DeviceMemorySpace> values) const override
  {
    CheckEvaluateWriteTag("MockEvaluator", field, values);
    if (dynamic_cast<const MockData*>(&field.GetData()) == nullptr) {
      throw pcms::pcms_error("MockEvaluator: not a field of this source");
    }
    const auto out = values.GetValues();
    const auto points = points_;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
        0, static_cast<LO>(points.extent(0))),
      KOKKOS_LAMBDA(const LO i) {
        for (size_t c = 0; c < out.extent(1); ++c) {
          out(i, c) =
            MockValue(points(i, 0), points(i, 1), points(i, 2), int(c));
        }
      });
  }

  Kokkos::View<const LO*, DeviceMemorySpace> FilledPoints() const override
  {
    return {};
  }

private:
  Kokkos::View<Real**, DeviceMemorySpace> points_;
};

class MockSpace final : public pcms::FunctionSpace
{
public:
  static std::shared_ptr<MockSpace> Create(int num_components)
  {
    return std::shared_ptr<MockSpace>(new MockSpace(num_components));
  }

  std::shared_ptr<const pcms::FieldLayout> GetLayout() const noexcept override
  {
    return layout_;
  }

protected:
  pcms::FieldVariant CreateFieldImpl(pcms::Type, ValueBasis) const override
  {
    throw pcms::pcms_error("MockSpace: pass the data naming the quantity");
  }

  pcms::FieldVariant CreateFieldImpl(pcms::FieldDataVariant data) const override
  {
    auto& real = std::get<std::unique_ptr<pcms::FieldData<Real>>>(data);
    if (dynamic_cast<const MockData*>(real.get()) == nullptr) {
      throw pcms::pcms_error("MockSpace: requires MockData");
    }
    return WrapField<Real>(layout_, std::move(real));
  }

  pcms::PointEvaluatorVariant CreatePointEvaluatorImpl(
    pcms::Type, const pcms::EvaluationRequest& request) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<MockEvaluator>(request));
  }

private:
  explicit MockSpace(int num_components)
    : layout_(std::make_shared<pcms::EmptyFieldLayout>(
        pcms::csys::CylindricalRThetaZ::Create(), num_components))
  {
  }

  std::shared_ptr<const pcms::EmptyFieldLayout> layout_;
};

std::shared_ptr<pcms::LagrangeFunctionSpace> MakePoloidalTarget(
  Omega_h::Mesh& mesh, int num_components)
{
  return pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::CylindricalRZ::Create(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
}

} // namespace

TEST_CASE("Evaluation-only space interpolates onto a poloidal plane")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildUnitSquare(lib, 1);
  const int num_components = GENERATE(1, 2);
  auto source = MockSpace::Create(num_components);
  auto target = MakePoloidalTarget(mesh, num_components);
  const Real theta0 = -0.3;
  pcms::TransformedSource on_plane(
    source, std::make_shared<pcms::PoloidalPlaneEmbedding>(theta0));
  pcms::Interpolator<Real> interpolator(on_plane, *target);

  auto f = source->CreateFunction<Real>("f", std::make_unique<MockData>());
  auto g = target->CreateFunction<Real>("g");
  interpolator.Apply(f, g);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    target->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto result = g.GetDOFHolderDataHost().GetValues();
  for (size_t i = 0; i < coords.extent(0); ++i) {
    for (int c = 0; c < num_components; ++c) {
      CAPTURE(i, c);
      REQUIRE_THAT(result(i, c),
                   WithinAbs(MockValue(coords(i, 0), theta0, coords(i, 1), c),
                             pcms::test::ExactTol));
    }
  }

  // The embedding maps points only, so vector values cannot pass through it.
  auto evaluator = on_plane.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromFunctionSpace(*target));
  REQUIRE_THROWS_WITH(
    evaluator->OutputBasis(ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                                      pcms::ComponentScaling::Physical,
                                      pcms::values::Vector}),
    ContainsSubstring("cannot re-express"));
}

TEST_CASE("Functions of an evaluation-only space hold no values")
{
  auto source = MockSpace::Create(1);
  auto f = source->CreateFunction<Real>("f", std::make_unique<MockData>());
  REQUIRE(f.GetLayout().GetNumOwnedDofHolder() == 0);
  REQUIRE_THROWS_WITH(f.GetDOFHolderData(), ContainsSubstring("no DOFAccess"));
  REQUIRE_THROWS_WITH(source->CreateFunction<Real>("default"),
                      ContainsSubstring("naming the quantity"));
}
