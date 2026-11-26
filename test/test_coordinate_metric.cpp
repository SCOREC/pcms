#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/function_space.hpp>
#include <pcms/field/layout/empty.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include "field_test_utils.hpp"
#include <numbers>
#include <optional>

using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::LO;
using pcms::Real;

namespace
{

// The metric and volume factor of `system` at the points, on the host.
struct BoundOnHost
{
  std::vector<Real> volume;
  std::vector<Real> metric; // [point][i][j]
  bool metric_broadcast;
};

BoundOnHost BindOnHost(const pcms::CoordinateSystem& system,
                       const std::vector<Real>& points)
{
  const int dim = system.Dimension();
  const auto n = static_cast<LO>(points.size() / dim);
  const auto coords = pcms::test::CreateDeviceRank2View(points, dim);
  const auto bound = system.Bind(pcms::MakeConstRank2View(coords));
  Kokkos::View<Real*, DeviceMemorySpace> volume("volume", n);
  Kokkos::View<Real***, DeviceMemorySpace> metric("metric", n, dim, dim);
  const auto v = bound->VolumeFactor();
  const auto g = bound->Metric();
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO p) {
      volume(p) = v(p);
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          metric(p, i, j) = g(p, i, j);
        }
      }
    });
  const auto v_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), volume);
  const auto g_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), metric);
  BoundOnHost out{{}, {}, g.stride(0) == 0};
  for (LO p = 0; p < n; ++p) {
    out.volume.push_back(v_host(p));
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        out.metric.push_back(g_host(p, i, j));
      }
    }
  }
  return out;
}

} // namespace

TEST_CASE("Bound coordinate systems give the metric and the volume factor",
          "[coordinate_system][metric]")
{
  SECTION("Cartesian: identity, allocated once")
  {
    const auto bound =
      BindOnHost(*pcms::csys::CreateCartesian(2), {0.1, 0.2, 3.0, -4.0});
    CHECK(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{1.0, 1.0});
    CHECK(bound.metric ==
          std::vector<Real>{1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0});
  }
  SECTION("cylindrical (R, theta, Z): diag(1, R^2, 1), volume factor R")
  {
    const auto system = pcms::csys::CreateCylindricalRThetaZ();
    const auto bound = BindOnHost(*system, {2.0, 0.3, 1.0, 3.0, -1.0, 0.0});
    CHECK_FALSE(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{2.0, 3.0});
    CHECK(bound.metric[4] == 4.0);
    CHECK(bound.metric[9 + 4] == 9.0);
    CHECK(bound.metric[1] == 0.0);
    CHECK(system->Period(1) == std::optional<Real>(2.0 * std::numbers::pi));
    CHECK_FALSE(system->Period(0).has_value());
  }
  SECTION("cylindrical (R, Z): identity metric, volume factor R")
  {
    const auto bound =
      BindOnHost(*pcms::csys::CreateCylindricalRZ(), {2.0, 0.3, 3.0, 1.0});
    CHECK(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{2.0, 3.0});
    CHECK(bound.metric ==
          std::vector<Real>{1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0});
  }
}

TEST_CASE("Parameter systems are labels without a metric",
          "[coordinate_system][metric]")
{
  const auto system = pcms::csys::CreateParameterSystem(
    "r_theta_zeta", 3, {std::nullopt, 2.0 * std::numbers::pi, 1.5});
  CHECK(system->Kind() == "r_theta_zeta");
  CHECK(system->Dimension() == 3);
  CHECK_FALSE(system->HasMetric());
  CHECK_FALSE(system->IsOrthonormal());
  CHECK_FALSE(system->Period(0).has_value());
  CHECK(system->Period(1) == std::optional<Real>(2.0 * std::numbers::pi));
  CHECK(system->Period(2) == std::optional<Real>(1.5));
  CHECK(pcms::csys::CreateCartesian(2)->HasMetric());
  const auto coords = pcms::test::CreateDeviceRank2View({0.1, 0.2, 0.3}, 3);
  CHECK_THROWS_AS(system->Bind(pcms::MakeConstRank2View(coords)),
                  pcms::pcms_error);
  // equal by identity
  CHECK(*system == *system);
  CHECK_FALSE(*system == *pcms::csys::CreateParameterSystem(
                           "r_theta_zeta", 3,
                           {std::nullopt, 2.0 * std::numbers::pi, 1.5}));
  CHECK_FALSE(*pcms::csys::CreateParameterSystem("x", 2) ==
              *pcms::csys::CreateCartesian(2));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 0));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 2, {1.0}));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 1, {-1.0}));
}

namespace
{

class Data final : public pcms::FieldData<Real>
{
public:
  pcms::FieldValueType GetValueType() const override
  {
    return pcms::FieldValueType::Scalar;
  }
  const pcms::ValueBasis& GetValueBasis() const override { return basis_; }

private:
  pcms::ValueBasis basis_;
};

class BoundAll final : public pcms::BoundPoints
{
public:
  BoundAll(std::shared_ptr<const pcms::Geometry> geometry,
           const pcms::BindRequest& request)
    : BoundPoints(std::move(geometry), request)
  {
  }

protected:
  void RebindImpl() override {}
};

// All of (R, theta, Z) space.
class Everywhere final : public pcms::Geometry
{
public:
  Everywhere() : system_(pcms::csys::CreateCylindricalRThetaZ()) {}
  const std::shared_ptr<const pcms::CoordinateSystem>& GetCoordinateSystem()
    const noexcept override
  {
    return system_;
  }

protected:
  std::shared_ptr<pcms::BoundPoints> BindImpl(
    const pcms::BindRequest& request) const override
  {
    return std::make_shared<BoundAll>(shared_from_this(), request);
  }

private:
  std::shared_ptr<const pcms::CoordinateSystem> system_;
};

// The partials of f = 2R + 3 theta - Z along the coordinates: (2, 3, -1).
class PartialsEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  PartialsEvaluator(std::shared_ptr<const pcms::BoundPoints> bound,
                    std::shared_ptr<const pcms::FieldLayout> layout)
    : PointEvaluator(std::move(bound), pcms::DifferentialOperator::Gradient,
                     std::move(layout))
  {
  }

  void EvaluateImpl(
    const pcms::Field<Real>&,
    pcms::ValueView<Real, DeviceMemorySpace> values) const override
  {
    const auto out = values.GetValues();
    pcms::detail::ForEach(
      "partials", static_cast<LO>(out.extent(0)), KOKKOS_LAMBDA(const LO p) {
        out(p, 0) = 2.0;
        out(p, 1) = 3.0;
        out(p, 2) = -1.0;
      });
  }
};

class PartialsFactory final : public pcms::PointEvaluatorFactory
{
public:
  PartialsFactory(std::shared_ptr<const pcms::Geometry> geometry,
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
    return op == pcms::DifferentialOperator::Gradient;
  }

protected:
  pcms::PointEvaluatorVariant CreatePointEvaluatorImpl(
    pcms::Type, std::shared_ptr<const pcms::BoundPoints> bound,
    pcms::DifferentialOperator) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<PartialsEvaluator>(std::move(bound), layout_));
  }

private:
  std::shared_ptr<const pcms::Geometry> geometry_;
  std::shared_ptr<const pcms::FieldLayout> layout_;
};

class PartialsSpace final : public pcms::FunctionSpace
{
public:
  using pcms::FunctionSpace::FunctionSpace;
};

} // namespace

TEST_CASE("Gradients in cylindrical coordinates are natural components; "
          "ToPhysical scales them",
          "[coordinate_system][metric][gradient]")
{
  const auto cyl = pcms::csys::CreateCylindricalRThetaZ();
  const auto layout = std::make_shared<pcms::field::EmptyFieldLayout>(cyl, 1);
  const auto geometry = std::make_shared<Everywhere>();
  const auto space = std::make_shared<PartialsSpace>(
    layout, std::make_shared<PartialsFactory>(geometry, layout), geometry);
  const auto f = space->CreateFunction<Real>("f", std::make_unique<Data>());
  const auto bound =
    pcms::test::Bind(*geometry, {2.0, 0.3, 1.0, 4.0, -1.0, 0.0});
  const auto evaluator = space->CreatePointEvaluator<Real>(
    bound, pcms::DifferentialOperator::Gradient);
  const auto natural = evaluator->OutputBasis(f.GetData().GetValueBasis());
  CHECK(natural.scaling == pcms::ComponentScaling::Natural);
  Kokkos::View<Real**, DeviceMemorySpace> values("gradient", 2, 3);
  evaluator->Evaluate(f, pcms::ValueView<Real, DeviceMemorySpace>(
                           natural, pcms::MakeRank2View(values)));
  const auto to_physical = pcms::detail::ToPhysical(cyl, bound->Metric());
  to_physical.Apply(
    pcms::ValueView<Real, DeviceMemorySpace>(natural,
                                             pcms::MakeRank2View(values)),
    pcms::ValueView<Real, DeviceMemorySpace>(to_physical.OutputBasis(natural),
                                             pcms::MakeRank2View(values)));
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), values);
  // (dR, dtheta / R, dZ)
  CHECK_THAT(host(0, 0), WithinAbs(2.0, pcms::test::ExactTol));
  CHECK_THAT(host(0, 1), WithinAbs(3.0 / 2.0, pcms::test::ExactTol));
  CHECK_THAT(host(0, 2), WithinAbs(-1.0, pcms::test::ExactTol));
  CHECK_THAT(host(1, 1), WithinAbs(3.0 / 4.0, pcms::test::ExactTol));
}
