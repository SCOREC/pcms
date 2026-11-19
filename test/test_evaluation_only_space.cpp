#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_library.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/field/layout/empty.hpp>
#include <pcms/geometry/geometry.hpp>
#include <pcms/transfer/interpolator.hpp>
#include "field_test_utils.hpp"

namespace pcms
{
using field::EmptyFieldLayout;
using omega_h::CreateGeometry;
} // namespace pcms

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::LO;
using pcms::Real;
using pcms::ValueBasis;

// An evaluation-only source, the pattern external codes (NIMROD, Fusion-IO)
// follow: a FunctionSpace on an EmptyFieldLayout whose field data names a
// quantity but provides no DOFAccess, created on a geometry of its own. Here
// the quantity is analytic: component c of f at (R, theta, Z) is
// 2R + 3 theta - Z + c, and the geometry is all of space.
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

// The mock evaluates every query point, so it needs no state of its own.
class MockBound final : public pcms::BoundPoints
{
public:
  MockBound(std::shared_ptr<const pcms::Geometry> geometry,
            const pcms::BindRequest& request)
    : BoundPoints(std::move(geometry), request)
  {
  }

protected:
  // Reset has marked every point found
  void RebindImpl() override {}
};

// All of (R, theta, Z) space: every point is found, without a search.
class AnalyticGeometry final : public pcms::Geometry
{
public:
  AnalyticGeometry() : system_(pcms::csys::CreateCylindricalRThetaZ()) {}

  const std::shared_ptr<const pcms::CoordinateSystem>& GetCoordinateSystem()
    const noexcept override
  {
    return system_;
  }

protected:
  std::shared_ptr<pcms::BoundPoints> BindImpl(
    const pcms::BindRequest& request) const override
  {
    return std::make_shared<MockBound>(shared_from_this(), request);
  }

private:
  std::shared_ptr<const pcms::CoordinateSystem> system_;
};

class MockEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  MockEvaluator(std::shared_ptr<const pcms::BoundPoints> bound,
                std::shared_ptr<const pcms::FieldLayout> layout)
    : PointEvaluator(std::move(bound), pcms::DifferentialOperator::Value,
                     std::move(layout))
  {
  }

  void EvaluateImpl(
    const pcms::Field<Real>& field,
    pcms::ValueView<Real, DeviceMemorySpace> values) const override
  {
    if (dynamic_cast<const MockData*>(&field.GetData()) == nullptr) {
      throw pcms::pcms_error("MockEvaluator: not a field of this source");
    }
    const auto out = values.GetValues();
    const auto points = Bound().LocatedPointValues();
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
};

// The evaluator source of the mock: it hands out MockEvaluators on the
// analytic geometry.
class MockEvaluatorFactory final : public pcms::PointEvaluatorFactory
{
public:
  MockEvaluatorFactory(std::shared_ptr<const pcms::Geometry> geometry,
                       std::shared_ptr<const pcms::FieldLayout> layout)
    : geometry_(std::move(geometry)), layout_(std::move(layout))
  {
  }

  std::shared_ptr<const pcms::Geometry> GetGeometry() const noexcept override
  {
    return geometry_;
  }

  bool Supports(pcms::DifferentialOperator op) const noexcept override
  {
    return op == pcms::DifferentialOperator::Value;
  }

protected:
  pcms::PointEvaluatorVariant CreatePointEvaluatorImpl(
    pcms::Type, std::shared_ptr<const pcms::BoundPoints> bound,
    pcms::DifferentialOperator /*op*/) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<MockEvaluator>(std::move(bound), layout_));
  }

private:
  std::shared_ptr<const pcms::Geometry> geometry_;
  std::shared_ptr<const pcms::FieldLayout> layout_;
};

// A space on an empty layout whose fields are made from data naming the
// quantity (MockData), evaluated by the mock on the analytic geometry.
std::shared_ptr<const pcms::FunctionSpace> MakeMockSpace(int num_components)
{
  const auto layout = std::make_shared<pcms::EmptyFieldLayout>(
    pcms::csys::CreateCylindricalRThetaZ(), num_components);
  const auto geometry = std::make_shared<AnalyticGeometry>();
  return std::make_shared<pcms::FunctionSpace>(
    layout,
    [](pcms::Type, ValueBasis) -> pcms::FieldDataVariant {
      throw pcms::pcms_error("MockSpace: pass the data naming the quantity");
    },
    std::make_shared<MockEvaluatorFactory>(geometry, layout), geometry);
}

std::shared_ptr<const pcms::FunctionSpace> MakePoloidalTarget(
  Omega_h::Mesh& mesh, int num_components)
{
  return pcms::CreateLagrangeSpace(
    pcms::CreateGeometry(mesh, pcms::csys::CreateCylindricalRZ()), 1,
    num_components,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
}

} // namespace

TEST_CASE("Evaluation-only space interpolates onto a poloidal plane")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildUnitSquare(lib, 1);
  const int num_components = GENERATE(1, 2);
  auto source = MakeMockSpace(num_components);
  auto target = MakePoloidalTarget(mesh, num_components);
  const Real theta0 = -0.3;
  const auto target_layout = target->GetLayout();
  const auto on_plane = source->GetGeometry()->Bind(
    target_layout->GetDOFHolderCoordinates(), pcms::OutOfBoundsPolicy{},
    pcms::CreatePoloidalPlaneEmbedding(theta0));
  const auto interpolator =
    pcms::CreateInterpolator<Real>(source, target, on_plane);

  auto f = source->CreateFunction<Real>("f", std::make_unique<MockData>());
  auto g = target->CreateFunction<Real>("g");
  interpolator->Apply(f, g);

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
  auto evaluator = source->CreatePointEvaluator<Real>(on_plane);
  REQUIRE_THROWS_WITH(
    evaluator->OutputBasis(ValueBasis{pcms::csys::CreateCylindricalRThetaZ(),
                                      pcms::ComponentScaling::Physical,
                                      pcms::values::Vector}),
    ContainsSubstring("cannot re-express"));
}

TEST_CASE("Functions of an evaluation-only space hold no values")
{
  auto source = MakeMockSpace(1);
  auto f = source->CreateFunction<Real>("f", std::make_unique<MockData>());
  REQUIRE(f.GetLayout().GetNumDofHolder() == 0);
  REQUIRE_THROWS_WITH(f.GetDOFHolderData(), ContainsSubstring("no DOFAccess"));
  REQUIRE_THROWS_WITH(source->CreateFunction<Real>("default"),
                      ContainsSubstring("naming the quantity"));
}

TEST_CASE("FunctionSpace requires the geometry of its evaluator source")
{
  const auto layout = std::make_shared<pcms::EmptyFieldLayout>(
    pcms::csys::CreateCylindricalRThetaZ(), 1);
  const auto create_data = [](pcms::Type,
                              ValueBasis) -> pcms::FieldDataVariant {
    throw pcms::pcms_error("no data");
  };
  const auto geometry = std::make_shared<AnalyticGeometry>();
  const auto evaluator =
    std::make_shared<MockEvaluatorFactory>(geometry, layout);
  REQUIRE_THROWS_WITH(
    pcms::FunctionSpace(layout, create_data, evaluator, nullptr),
    Catch::Matchers::ContainsSubstring("the geometry must not be null"));
  REQUIRE_THROWS_WITH(
    pcms::FunctionSpace(layout, create_data, nullptr, nullptr),
    Catch::Matchers::ContainsSubstring("the geometry must not be null"));
  REQUIRE_THROWS_WITH(
    pcms::FunctionSpace(layout, create_data, evaluator,
                        std::make_shared<AnalyticGeometry>()),
    Catch::Matchers::ContainsSubstring("the evaluator source's geometry"));
}
