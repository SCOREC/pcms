#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "embedding_test_utils.hpp"
#include <numbers>
#include <random>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using pcms::DeviceMemorySpace;
using pcms::Real;
using pcms::test::Annulus;

namespace
{

// the analytic Jacobian J(j, i) = d x_j / d xi_i of the annulus
std::array<std::array<Real, 2>, 2> AnnulusJacobian(Real xi1, Real xi2)
{
  const Real rho = Annulus::Rho(xi1);
  const Real theta = Annulus::Theta(xi2);
  return {{{0.5 * std::cos(theta), -2.0 * rho * std::sin(theta)},
           {0.5 * std::sin(theta), 2.0 * rho * std::cos(theta)}}};
}

const std::vector<Real> Params = {0.31, 0.42, 0.77, 0.13, 0.5, 0.95, 0.05, 0.6};
// the position is a spline of 32 x 32 cells: a smooth map is reproduced to
// roughly its fourth derivative times h^4
constexpr Real SplineTol = 2e-6;

struct Fixture
{
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakeSplinePosition(pcms::test::MakeUnitSquareGrid(32),
                                   Annulus::Position);
  std::shared_ptr<const pcms::CoordinateMap> map =
    pcms::CreateCoordinateMap(position);
  pcms::test::DeviceCoordinates xi = pcms::test::CreateDeviceCoordinateView(
    Params, map->GetSourceCoordinateSystem(), 2);
};

} // namespace

TEST_CASE("The coordinate map of a position maps parameter points and gives "
          "its Jacobian",
          "[embedding]")
{
  const Fixture f;
  CHECK(f.map->GetTargetCoordinateSystem()->Kind() == "cylindrical_rz");
  CHECK(f.map->GetSourceCoordinateSystem()->Kind() == "xi");
  const auto bound = f.map->Bind(f.xi.coordinate_view);
  const auto mapped = pcms::test::CopyCoordinatesToHost(bound->Mapped());
  const auto jacobian_view = bound->Jacobian();
  Kokkos::View<Real***, pcms::HostMemorySpace> jacobian("j", 4, 2, 2);
  Kokkos::View<Real***, DeviceMemorySpace> device("j", 4, 2, 2);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 4),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          device(p, i, j) = jacobian_view(p, i, j);
        }
      }
    });
  Kokkos::deep_copy(jacobian, device);
  for (int p = 0; p < 4; ++p) {
    const Real xi1 = Params[std::size_t{2} * p];
    const Real xi2 = Params[2 * p + 1];
    const auto x = Annulus::Position(xi1, xi2);
    const auto j = AnnulusJacobian(xi1, xi2);
    CAPTURE(p);
    CHECK_THAT(mapped(p, 0), WithinAbs(x[0], SplineTol));
    CHECK_THAT(mapped(p, 1), WithinAbs(x[1], SplineTol));
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        CHECK_THAT(jacobian(p, a, b), WithinAbs(j[a][b], 100 * SplineTol));
      }
    }
  }
}

TEST_CASE("Coordinates pass through value transforms", "[embedding]")
{
  const Fixture f;
  const auto bound = f.map->Bind(f.xi.coordinate_view);
  const auto coordinates =
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCartesian(3));
  CHECK(bound->ValueTransforms(coordinates).empty());
  const pcms::ValueBasis vector{pcms::csys::CreateCylindricalRZ(),
                                pcms::ComponentScaling::Physical,
                                pcms::values::Vector};
  for (const auto* transform : bound->ValueTransforms(vector)) {
    CHECK(transform->OutputBasis(coordinates) == coordinates);
  }
  // coordinates are not scalars, and coordinates of different coordinate
  // systems differ
  CHECK_FALSE(pcms::detail::SameValueBasis(coordinates, pcms::ValueBasis{}));
  CHECK_FALSE(pcms::detail::SameValueBasis(
    coordinates,
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCylindricalRThetaZ())));
  CHECK(pcms::detail::SameValueBasis(
    f.position->GetData().GetValueBasis(),
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCylindricalRZ())));
}

TEST_CASE("A function of several components defaults to independent "
          "components",
          "[embedding]")
{
  const auto space =
    pcms::CreateLagrangeSpace(pcms::test::MakeUnitSquareGrid(2), 1, 2);
  CHECK(space->CreateFunction<Real>().GetData().GetValueBasis().kind ==
        pcms::ValueKind::Array);
  CHECK_THROWS_WITH(space->CreateFunction<Real>("x", pcms::ValueBasis{}),
                    ContainsSubstring("a scalar has one component"));
}

TEST_CASE("The coordinate map of a position refuses what it cannot do",
          "[embedding]")
{
  const auto geometry = pcms::test::MakeUnitSquareGrid(4);
  CHECK_THROWS_WITH(pcms::CreateCoordinateMap(nullptr),
                    ContainsSubstring("must not be null"));
  // values that are not coordinates
  const auto space = pcms::CreateSplineSpace(geometry, 2);
  const auto scalars = std::make_shared<pcms::Function<Real>>(
    space->CreateFunction<Real>("x", {.kind = pcms::ValueKind::Array}));
  CHECK_THROWS_WITH(pcms::CreateCoordinateMap(scalars),
                    ContainsSubstring("must be coordinates"));
  // two components are not the coordinates of a 3D coordinate system
  CHECK_THROWS_WITH(
    space->CreateFunction<Real>(
      "x", pcms::ValueBasis::Coordinates(pcms::csys::CreateCartesian(3))),
    ContainsSubstring("require 3 components"));
}

namespace
{

// An embedded annulus sector, and a P1 space on its parameter grid put on it.
struct EmbeddedFixture
{
  std::shared_ptr<const pcms::Geometry> grid =
    pcms::test::MakeUnitSquareGrid(32);
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakeSplinePosition(grid, Annulus::Position);
  std::shared_ptr<const pcms::EmbeddedGeometry> embedded =
    pcms::CreateEmbeddedGeometry(grid, position);
  std::shared_ptr<const pcms::FunctionSpace> space =
    pcms::CreateEmbeddedSpace(embedded, pcms::CreateLagrangeSpace(grid, 1, 1));
};

// f = xi1 + 2 xi2 at the points of `bound`
std::vector<Real> EvaluateLinearInParameters(
  const EmbeddedFixture& f,
  const std::shared_ptr<const pcms::BoundPoints>& bound)
{
  auto field = f.space->CreateFunction<Real>("f");
  pcms::test::SetField(
    field, KOKKOS_LAMBDA(Real xi1, Real xi2) { return xi1 + 2 * xi2; });
  const auto evaluator = f.space->CreatePointEvaluator<Real>(bound);
  Kokkos::View<Real**, DeviceMemorySpace> out("out", bound->NumPoints(), 1);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out));
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
  return std::vector<Real>(host.data(), host.data() + host.extent(0));
}

// The geometric factors at the points of `bound`, on the host: J, g
// ([point][4] each) and the volume factor.
struct HostFactors
{
  Kokkos::View<Real**, pcms::DeviceMemorySpace>::host_mirror_type jacobian;
  Kokkos::View<Real**, pcms::DeviceMemorySpace>::host_mirror_type metric;
  Kokkos::View<Real*, pcms::DeviceMemorySpace>::host_mirror_type volume;
};

HostFactors FactorsOnHost(const pcms::BoundPoints& bound)
{
  const auto n = bound.NumPoints();
  Kokkos::View<Real**, DeviceMemorySpace> jacobian("j", n, 4);
  Kokkos::View<Real**, DeviceMemorySpace> metric("g", n, 4);
  Kokkos::View<Real*, DeviceMemorySpace> volume("v", n);
  const auto j = bound.Jacobian();
  const auto g = bound.Metric();
  const auto v = bound.VolumeFactor();
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
          jacobian(p, 2 * a + b) = j(p, a, b);
          metric(p, 2 * a + b) = g(p, a, b);
        }
      }
      volume(p) = v(p);
    });
  return {
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), jacobian),
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), metric),
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), volume)};
}

} // namespace

TEST_CASE("An embedded geometry locates physical points by inverting the "
          "position",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  CHECK(f.embedded->GetCoordinateSystem()->Kind() == "cylindrical_rz");
  CHECK(f.embedded->TopologicalDimension() == 2);
  CHECK(*f.embedded->GetParameterSystem() == *f.grid->GetCoordinateSystem());
  CHECK(f.embedded->GetPosition() == f.position);
  CHECK(f.embedded->GetMetric() == nullptr);

  std::vector<Real> physical;
  std::vector<Real> expected;
  for (int i = 0; i < 4; ++i) {
    const auto x = Annulus::Position(Params[std::size_t{2} * i],
                                     Params[std::size_t{2} * i + 1]);
    physical.insert(physical.end(), x.begin(), x.end());
    expected.push_back(Params[std::size_t{2} * i] +
                       2 * Params[std::size_t{2} * i + 1]);
  }
  const auto bound = pcms::test::Bind(*f.embedded, physical);
  CHECK(bound->NumPoints() == 4);
  // the parameters are recovered to the error of the position's spline
  const auto values = EvaluateLinearInParameters(f, bound);
  for (int i = 0; i < 4; ++i) {
    CAPTURE(i);
    CHECK_THAT(values[i], WithinAbs(expected[i], 1e-5));
  }
}

TEST_CASE("The geometric factors of an embedded geometry are the position's "
          "Jacobian and the induced metric",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  const auto factors = FactorsOnHost(*bound);
  for (int p = 0; p < 4; ++p) {
    const Real xi1 = Params[std::size_t{2} * p];
    const Real xi2 = Params[std::size_t{2} * p + 1];
    const Real rho = Annulus::Rho(xi1);
    const auto x = Annulus::Position(xi1, xi2);
    const auto j = AnnulusJacobian(xi1, xi2);
    CAPTURE(p);
    for (int k = 0; k < 4; ++k) {
      CHECK_THAT(factors.jacobian(p, k),
                 WithinAbs(j[k / 2][k % 2], 100 * SplineTol));
    }
    // g = diag(0.25, 4 rho^2)
    CHECK_THAT(factors.metric(p, 0), WithinAbs(0.25, 100 * SplineTol));
    CHECK_THAT(factors.metric(p, 3), WithinAbs(4 * rho * rho, 100 * SplineTol));
    CHECK_THAT(factors.metric(p, 1), WithinAbs(0.0, 100 * SplineTol));
    // sqrt(det g) = rho, times the volume factor of (R, Z), R
    CHECK_THAT(factors.volume(p), WithinAbs(rho * x[0], 100 * SplineTol));
  }
}

TEST_CASE("Parameter-system vectors and covectors of an embedded space are "
          "refused",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  const auto space = pcms::CreateEmbeddedSpace(
    f.embedded, pcms::CreateLagrangeSpace(f.grid, 1, 2));
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  const auto evaluator = space->CreatePointEvaluator<Real>(bound);
  for (const auto variance : {pcms::values::Vector, pcms::values::Covector}) {
    const auto field = space->CreateFunction<Real>(
      "e1", {.system = f.embedded->GetParameterSystem(),
             .scaling = pcms::ComponentScaling::Natural,
             .variance = variance});
    CHECK_THROWS_WITH(evaluator->OutputBasis(field.GetData().GetValueBasis()),
                      ContainsSubstring("not yet supported"));
  }
}

namespace
{

// The evaluation at point p is the DOFs of holder p: a linear map with a
// transpose, on any geometry with as many DOF holders as points.
class SampleEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  SampleEvaluator(std::shared_ptr<const pcms::BoundPoints> bound,
                  std::shared_ptr<const pcms::FieldLayout> layout)
    : PointEvaluator(std::move(bound), pcms::DifferentialOperator::Value,
                     std::move(layout))
  {
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }

  void EvaluateImpl(
    const pcms::Field<Real>& field,
    pcms::ValueView<Real, DeviceMemorySpace> values) const override
  {
    const auto dofs = field.GetDOFHolderData().GetValues();
    const auto out = values.GetValues();
    const int nc = static_cast<int>(out.extent(1));
    pcms::detail::ForEach(
      "sample", NumPoints(), KOKKOS_LAMBDA(const pcms::LO p) {
        for (int c = 0; c < nc; ++c) {
          out(p, c) = dofs(p, c);
        }
      });
  }

  void EvaluateTransposeImpl(
    pcms::Rank2View<const Real, DeviceMemorySpace> values,
    pcms::Field<Real>& functional) const override
  {
    const auto current = functional.GetDOFHolderData();
    const auto old = current.GetValues();
    Kokkos::View<Real**, DeviceMemorySpace> sum("sum", old.extent(0),
                                                old.extent(1));
    const int nc = static_cast<int>(old.extent(1));
    const pcms::LO np = NumPoints();
    pcms::detail::ForEach(
      "sample_transpose", static_cast<pcms::LO>(old.extent(0)),
      KOKKOS_LAMBDA(const pcms::LO h) {
        for (int c = 0; c < nc; ++c) {
          sum(h, c) = old(h, c) + (h < np ? values(h, c) : 0.0);
        }
      });
    functional.SetDOFHolderData(pcms::ValueView<const Real, DeviceMemorySpace>(
      current.GetBasis(), pcms::MakeConstRank2View(sum)));
  }
};

class SampleFactory final : public pcms::PointEvaluatorFactory
{
public:
  explicit SampleFactory(std::shared_ptr<const pcms::FunctionSpace> space)
    : space_(std::move(space))
  {
  }
  std::shared_ptr<const pcms::Geometry> GetGeometry() const noexcept override
  {
    return space_->GetGeometry();
  }
  bool Supports(pcms::DifferentialOperator op) const noexcept override
  {
    return op == pcms::DifferentialOperator::Value;
  }

protected:
  pcms::PointEvaluatorVariant CreatePointEvaluatorImpl(
    pcms::Type, std::shared_ptr<const pcms::BoundPoints> bound,
    pcms::DifferentialOperator) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<SampleEvaluator>(std::move(bound), space_->GetLayout()));
  }

private:
  std::shared_ptr<const pcms::FunctionSpace> space_;
};

// Random values, on the device.
Kokkos::View<Real**, DeviceMemorySpace> RandomDeviceValues(std::size_t rows,
                                                           std::size_t cols,
                                                           std::mt19937& rng)
{
  std::uniform_real_distribution<Real> uniform(-1.0, 1.0);
  Kokkos::View<Real**, DeviceMemorySpace> values("random", rows, cols);
  auto host = Kokkos::create_mirror_view(values);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t c = 0; c < cols; ++c) {
      host(i, c) = uniform(rng);
    }
  }
  Kokkos::deep_copy(values, host);
  return values;
}

// The plain pairing of two host views of the same shape.
template <typename A, typename B>
Real Pairing(const A& a, const B& b)
{
  Real sum = 0;
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    for (std::size_t c = 0; c < a.extent(1); ++c) {
      sum += a(i, c) * b(i, c);
    }
  }
  return sum;
}

} // namespace

TEST_CASE("EvaluateTranspose through an embedded space is the transpose of "
          "Evaluate",
          "[embedding][embedded_geometry][adjoint]")
{
  const EmbeddedFixture f;
  const auto grid_space = pcms::CreateLagrangeSpace(f.grid, 1, 2);
  const auto sampled = std::make_shared<pcms::FunctionSpace>(
    grid_space->GetLayout(), std::make_shared<SampleFactory>(grid_space),
    f.grid);
  const auto space = pcms::CreateEmbeddedSpace(f.embedded, sampled);
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  const auto evaluator = space->CreatePointEvaluator<Real>(bound);
  REQUIRE(evaluator->SupportsTranspose());
  std::mt19937 rng(23);
  // a scalar field of two components
  const pcms::ValueBasis basis_in{.kind = pcms::ValueKind::Array};
  auto x = space->CreateFunction<Real>("x", basis_in);
  auto functional = space->CreateFunction<Real>("functional", basis_in);
  const auto basis = evaluator->OutputBasis(x.GetData().GetValueBasis());
  const auto dofs =
    RandomDeviceValues(grid_space->GetLayout()->GetNumDofHolder(), 2, rng);
  x.SetDOFHolderData(pcms::ValueView<const Real, DeviceMemorySpace>(
    x.GetData().GetValueBasis(), pcms::MakeConstRank2View(dofs)));
  const auto y = RandomDeviceValues(4, 2, rng);
  Kokkos::View<Real**, DeviceMemorySpace> ex("Ex", 4, 2);
  evaluator->Evaluate(x, pcms::ValueView<Real, DeviceMemorySpace>(
                           basis, pcms::MakeRank2View(ex)));
  evaluator->EvaluateTranspose(pcms::ValueView<const Real, DeviceMemorySpace>(
                                 basis, pcms::MakeConstRank2View(y)),
                               functional);
  const auto host = [](const Kokkos::View<Real**, DeviceMemorySpace>& v) {
    return Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), v);
  };
  const Real lhs = Pairing(host(ex), host(y));
  const Real rhs =
    Pairing(host(dofs), functional.GetDOFHolderDataHost().GetValues());
  CAPTURE(lhs, rhs);
  CHECK_THAT(lhs, WithinAbs(rhs, 1e-12 * (1 + std::abs(rhs))));
}

TEST_CASE("An embedded geometry binds points of its parameter system without "
          "a search",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    Params, f.embedded->GetParameterSystem(), 2);
  const auto bound = f.embedded->Bind(coords.coordinate_view);
  // the field is evaluated at exactly the parameters given
  const auto values = EvaluateLinearInParameters(f, bound);
  for (int i = 0; i < 4; ++i) {
    CHECK_THAT(values[i], WithinAbs(Params[std::size_t{2} * i] +
                                      2 * Params[std::size_t{2} * i + 1],
                                    pcms::test::ExactTol));
  }
  // outside the parameter grid
  const auto outside = pcms::test::CreateDeviceCoordinateView(
    {1.5, 0.5}, f.embedded->GetParameterSystem(), 2);
  CHECK_THROWS(f.embedded->Bind(outside.coordinate_view));
  const auto missing = f.embedded->Bind(
    outside.coordinate_view,
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
  CHECK(missing->MissingPoints().extent(0) == 1);
  CHECK(EvaluateLinearInParameters(f, missing)[0] == -7.0);
}

TEST_CASE("The located points of a bind in the parameter system are their "
          "images in the coordinate system",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  const auto located = bound->LocatedPoints();
  CHECK(located.GetCoordinateSystem() == f.embedded->GetCoordinateSystem());
  const auto host = pcms::test::CopyCoordinatesToHost(located.GetValues());
  REQUIRE(host.extent(1) == 2);
  for (int p = 0; p < 4; ++p) {
    CAPTURE(p);
    const auto x = Annulus::Position(Params[std::size_t{2} * p],
                                     Params[std::size_t{2} * p + 1]);
    CHECK_THAT(host(p, 0), WithinAbs(x[0], SplineTol));
    CHECK_THAT(host(p, 1), WithinAbs(x[1], SplineTol));
  }
}

TEST_CASE("A source on another geometry projects onto an embedded space",
          "[embedding][embedded_geometry][transfer]")
{
  const EmbeddedFixture f;
  pcms::UniformGrid<2> box;
  box.bot_left = {1.3, 0.0};
  box.edge_length = {1.8, 1.1};
  box.divisions = {36, 22};
  const auto plane = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(box, f.embedded->GetCoordinateSystem()),
    1, 1);
  auto source = plane->CreateFunction<Real>("f");
  pcms::test::SetFieldComponents(
    source, [](Real r, Real z, Real* out) { out[0] = r + 2 * z; });
  auto target = f.space->CreateFunction<Real>("f");
  pcms::ProjectionOptions options;
  options.quadrature = pcms::MapQuadrature(
    f.embedded,
    pcms::CreateGaussQuadrature(f.embedded->GetParameterGeometry(), 2));
  pcms::CreateProjection(plane, f.space, options)->Apply(source, target);
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  const auto values =
    pcms::test::EvaluateToHost(*f.space->CreatePointEvaluator<Real>(bound),
                               target, target.GetData().GetValueBasis());
  for (int p = 0; p < 4; ++p) {
    CAPTURE(p);
    const auto x = Annulus::Position(Params[std::size_t{2} * p],
                                     Params[std::size_t{2} * p + 1]);
    CHECK_THAT(values(p, 0), WithinAbs(x[0] + 2 * x[1], 5e-3));
  }
}

TEST_CASE("An embedded geometry applies the out-of-bounds policy to points "
          "outside its image",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  // inside the sector, the hole of the annulus, and outside
  const std::vector<Real> points = {
    Annulus::Rc + 0.75, 0.3, Annulus::Rc, 0.0, 5.0, 5.0};
  CHECK_THROWS_AS(pcms::test::Bind(*f.embedded, points), pcms::pcms_error);
  CHECK_THROWS_AS(pcms::test::Bind(*f.embedded, points,
                                   pcms::OutOfBoundsPolicy{
                                     pcms::OutOfBoundsMode::NEAREST_BOUNDARY}),
                  pcms::pcms_error);
  const auto bound = pcms::test::Bind(
    *f.embedded, points,
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
  CHECK(bound->FoundPoints().extent(0) == 1);
  CHECK(bound->MissingPoints().extent(0) == 2);
  const auto values = EvaluateLinearInParameters(f, bound);
  CHECK(values[1] == -7.0);
  CHECK(values[2] == -7.0);
  // far beyond the bins of the search
  const auto far = pcms::test::Bind(
    *f.embedded, {1e12, 0.0, 0.0, 1e300},
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
  CHECK(far->MissingPoints().extent(0) == 2);
}

TEST_CASE("The search of an embedded geometry sampled in batches binds as in "
          "one batch",
          "[embedding][embedded_geometry]")
{
  const std::vector<Real> points = {Annulus::Rc + 0.75,
                                    0.3,
                                    Annulus::Rc - 0.2,
                                    0.6,
                                    Annulus::Rc,
                                    0.0,
                                    5.0,
                                    5.0};
  const pcms::OutOfBoundsPolicy fill{pcms::OutOfBoundsMode::FILL, -7.0};
  const EmbeddedFixture whole;
  const auto expected = EvaluateLinearInParameters(
    whole, pcms::test::Bind(*whole.embedded, points, fill));
  const auto batch = pcms::detail::EmbeddingSampleBatch;
  pcms::detail::EmbeddingSampleBatch = 7;
  const EmbeddedFixture batched;
  pcms::detail::EmbeddingSampleBatch = batch;
  const auto bound = pcms::test::Bind(*batched.embedded, points, fill);
  CHECK(bound->MissingPoints().extent(0) == 2);
  CHECK(EvaluateLinearInParameters(batched, bound) == expected);
}

TEST_CASE("The quadrature of an embedded geometry integrates in the physical "
          "measure",
          "[embedding][embedded_geometry][quadrature]")
{
  const EmbeddedFixture f;
  const auto quadrature = pcms::MapQuadrature(
    f.embedded,
    pcms::CreateGaussQuadrature(f.embedded->GetParameterGeometry(), 4));
  CHECK(quadrature.points->NumPoints() == 32 * 32 * 9);
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(), quadrature.weights);
  Real sum = 0;
  for (size_t i = 0; i < weights.extent(0); ++i) {
    sum += weights(i);
  }
  // the area of the sector times the R of the volume of revolution
  CHECK_THAT(sum, WithinAbs(Annulus::MomentOfR(), 1e-5));
  CHECK_THROWS(pcms::CreateMonteCarloQuadrature(f.embedded, 3));
  // the points are bound on the embedded geometry
  CHECK(quadrature.points->GetGeometry() == f.embedded);
  CHECK_THROWS(pcms::MapQuadrature(f.embedded, quadrature));
}

TEST_CASE("An embedded geometry survives a singular point of the position",
          "[embedding][embedded_geometry]")
{
  // a disc of radius 1 around (2, 0) with the center a pole: rho = xi1,
  // theta = 2 pi xi2, J has the determinant 2 pi rho, zero on the edge xi1 = 0
  const auto position = [](Real xi1, Real xi2) {
    const Real theta = 2 * std::numbers::pi * xi2;
    return std::array<Real, 2>{2.0 + xi1 * std::cos(theta),
                               xi1 * std::sin(theta)};
  };
  const auto grid = pcms::test::MakeUnitSquareGrid(32);
  const auto embedded = pcms::CreateEmbeddedGeometry(
    grid, pcms::test::MakeSplinePosition(grid, position));
  // the pole and a regular point; the policy fills what is not found
  const auto bound = pcms::test::Bind(
    *embedded, {2.0, 0.0, 2.0, 0.5},
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
  CHECK(bound->NumPoints() == 2);
  // the regular point, xi = (0.5, 0.25), is found
  const auto found = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(), bound->FoundPoints());
  bool regular_found = false;
  for (size_t i = 0; i < found.extent(0); ++i) {
    regular_found = regular_found || found(i) == 1;
  }
  CHECK(regular_found);
}

TEST_CASE("Embedded geometries refuse what they cannot do",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  CHECK_THROWS_WITH(pcms::CreateEmbeddedGeometry(nullptr, f.position),
                    ContainsSubstring("a Function on the parameter geometry"));
  // a position on another geometry
  CHECK_THROWS_WITH(
    pcms::CreateEmbeddedGeometry(pcms::test::MakeUnitSquareGrid(4), f.position),
    ContainsSubstring("a Function on the parameter geometry"));
  // the parameter system is the position's coordinate system
  const auto plane =
    pcms::test::MakeUnitSquareGrid(4, pcms::csys::CreateCylindricalRZ());
  CHECK_THROWS_WITH(
    pcms::CreateEmbeddedGeometry(
      plane, pcms::test::MakeSplinePosition(plane, Annulus::Position)),
    ContainsSubstring("must differ"));
  CHECK_THROWS_WITH(pcms::CreateEmbeddedSpace(
                      f.embedded, pcms::CreateLagrangeSpace(
                                    pcms::test::MakeUnitSquareGrid(4), 1, 1)),
                    ContainsSubstring("embedded geometry of the space"));
  const auto bound = pcms::test::Bind(*f.embedded, {Annulus::Rc + 0.75, 0.3});
  CHECK_FALSE(f.space->Supports(pcms::DifferentialOperator::Gradient));
  // Rebind locates new points
  const auto again = f.embedded->Rebind(
    bound, pcms::test::CreateDeviceCoordinateView(
             {Annulus::Rc + 0.8, 0.2}, f.embedded->GetCoordinateSystem(), 2)
             .coordinate_view);
  CHECK(again->NumPoints() == 1);
}

namespace
{

// The affine map x = A xi + b into Cartesian (x, y): its induced metric is
// A^T A, constant, so a P1 metric Function represents it exactly.
constexpr Real A[2][2] = {{2.0, 0.5}, {0.3, 1.5}};

std::array<Real, 2> Affine(Real xi1, Real xi2)
{
  return {1.0 + A[0][0] * xi1 + A[0][1] * xi2,
          -2.0 + A[1][0] * xi1 + A[1][1] * xi2};
}

// A^T A, scaled
std::array<Real, 4> InducedMetric(Real scale)
{
  std::array<Real, 4> g{};
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      g[2 * a + b] = scale * (A[0][a] * A[0][b] + A[1][a] * A[1][b]);
    }
  }
  return g;
}

// A constant metric Function on the grid: `g` (covariant) or its inverse
// (contravariant).
std::shared_ptr<const pcms::Function<Real>> ConstantMetric(
  const std::shared_ptr<const pcms::Geometry>& grid, std::array<Real, 4> g,
  bool contravariant)
{
  const auto space = pcms::CreateLagrangeSpace(grid, 1, 4);
  const auto variance =
    contravariant ? pcms::VarianceSignature{pcms::Variance::Contravariant,
                                            pcms::Variance::Contravariant}
                  : pcms::VarianceSignature{pcms::Variance::Covariant,
                                            pcms::Variance::Covariant};
  if (contravariant) {
    const Real det = g[0] * g[3] - g[1] * g[2];
    g = {g[3] / det, -g[1] / det, -g[2] / det, g[0] / det};
  }
  auto metric =
    std::make_shared<pcms::Function<Real>>(space->CreateFunction<Real>(
      "g", {.scaling = pcms::ComponentScaling::Natural, .variance = variance}));
  pcms::test::SetFieldComponents(*metric, [g](Real, Real, Real* out) {
    for (int k = 0; k < 4; ++k) {
      out[k] = g[k];
    }
  });
  return metric;
}

std::vector<Real> QuadratureWeights(
  const std::shared_ptr<const pcms::EmbeddedGeometry>& geometry)
{
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(),
    pcms::MapQuadrature(geometry, pcms::CreateGaussQuadrature(
                                    geometry->GetParameterGeometry(), 2))
      .weights);
  return std::vector<Real>(weights.data(), weights.data() + weights.size());
}

} // namespace

TEST_CASE("An embedded geometry uses a given metric",
          "[embedding][embedded_geometry][metric]")
{
  const auto grid = pcms::test::MakeUnitSquareGrid(4);
  const auto position = pcms::test::MakeSplinePosition(
    grid, Affine, pcms::csys::CreateCartesian(2));
  const auto induced = pcms::CreateEmbeddedGeometry(grid, position);
  const auto induced_weights = QuadratureWeights(induced);
  Real area = 0;
  for (const Real w : induced_weights) {
    area += w;
  }
  // |det A|
  CHECK_THAT(area, WithinAbs(2.85, 1e-12));

  for (const bool contravariant : {false, true}) {
    CAPTURE(contravariant);
    // the induced metric, given, reproduces the induced weights
    const auto metric = ConstantMetric(grid, InducedMetric(1.0), contravariant);
    const auto same =
      pcms::CreateEmbeddedGeometry(grid, position, {.metric = metric});
    CHECK(same->GetMetric() == metric);
    const auto same_weights = QuadratureWeights(same);
    REQUIRE(same_weights.size() == induced_weights.size());
    for (std::size_t i = 0; i < same_weights.size(); ++i) {
      CHECK_THAT(same_weights[i], WithinAbs(induced_weights[i], 1e-12));
    }
    // another metric is the one used: four times the metric, twice the
    // lengths, four times the areas
    const auto other = pcms::CreateEmbeddedGeometry(
      grid, position,
      {.metric = ConstantMetric(grid, InducedMetric(4.0), contravariant)});
    const auto weights = QuadratureWeights(other);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK_THAT(weights[i], WithinAbs(4 * induced_weights[i], 1e-12));
    }
    const auto bound = other->Bind(pcms::test::CreateDeviceCoordinateView(
                                     {0.3, 0.6}, other->GetParameterSystem(), 2)
                                     .coordinate_view);
    const auto factors = FactorsOnHost(*bound);
    const auto expected = InducedMetric(4.0);
    for (int k = 0; k < 4; ++k) {
      CHECK_THAT(factors.metric(0, k), WithinAbs(expected[k], 1e-12));
    }
  }
  SECTION("a metric that is not one of the parameter system is refused")
  {
    const auto space = pcms::CreateLagrangeSpace(grid, 1, 4);
    const auto mixed =
      std::make_shared<pcms::Function<Real>>(space->CreateFunction<Real>(
        "g", {.scaling = pcms::ComponentScaling::Natural,
              .variance = pcms::VarianceSignature{pcms::Variance::Contravariant,
                                                  pcms::Variance::Covariant}}));
    CHECK_THROWS_WITH(
      pcms::CreateEmbeddedGeometry(grid, position, {.metric = mixed}),
      ContainsSubstring("rank 2 in the parameter system"));
    const auto other = ConstantMetric(pcms::test::MakeUnitSquareGrid(2),
                                      InducedMetric(1.0), false);
    CHECK_THROWS_WITH(
      pcms::CreateEmbeddedGeometry(grid, position, {.metric = other}),
      ContainsSubstring("rank 2 in the parameter system"));
  }
}

TEST_CASE("An embedded geometry uses a given volume factor",
          "[embedding][embedded_geometry][metric]")
{
  const auto grid = pcms::test::MakeUnitSquareGrid(4);
  const auto position = pcms::test::MakeSplinePosition(
    grid, Affine, pcms::csys::CreateCartesian(2));
  const auto scalars = pcms::CreateLagrangeSpace(grid, 1, 1);
  // the integral of 1 on the embedded geometry with `options`
  const auto area = [&](pcms::EmbeddingOptions options) {
    const auto space = pcms::CreateEmbeddedSpace(
      pcms::CreateEmbeddedGeometry(grid, position, std::move(options)),
      scalars);
    auto one = space->CreateFunction<Real>("one");
    pcms::test::SetFieldComponents(one,
                                   [](Real, Real, Real* out) { out[0] = 1; });
    return pcms::L2Integral(space).Apply(one);
  };
  auto volume = std::make_shared<pcms::Function<Real>>(
    scalars->CreateFunction<Real>("volume"));
  pcms::test::SetFieldComponents(*volume,
                                 [](Real, Real, Real* out) { out[0] = 3.5; });
  for (const auto& metric : {std::shared_ptr<const pcms::Function<Real>>{},
                             ConstantMetric(grid, InducedMetric(4.0), true)}) {
    CAPTURE(metric != nullptr);
    CHECK(pcms::CreateEmbeddedGeometry(grid, position,
                                       {.metric = metric, .volume = volume})
            ->GetVolume() == volume);
    // 3.5 times the unit parameter area
    CHECK_THAT(area({.metric = metric, .volume = volume}),
               WithinAbs(3.5, 1e-12));
  }
  // a volume per cell: the mean of 1 + x + 2 y at the 16 cell centers
  auto cells = std::make_shared<pcms::Function<Real>>(
    pcms::CreateLagrangeSpace(grid, 0, 1)->CreateFunction<Real>("volume"));
  pcms::test::SetFieldComponents(
    *cells, [](Real x, Real y, Real* out) { out[0] = 1 + x + 2 * y; });
  CHECK_THAT(area({.volume = cells}), WithinAbs(2.5, 1e-12));
  const auto vectors = pcms::CreateLagrangeSpace(grid, 1, 2);
  const auto vector =
    std::make_shared<pcms::Function<Real>>(vectors->CreateFunction<Real>(
      "volume", {.scaling = pcms::ComponentScaling::Natural,
                 .variance = pcms::values::Vector}));
  CHECK_THROWS_WITH(
    pcms::CreateEmbeddedGeometry(grid, position, {.volume = vector}),
    ContainsSubstring("the volume factor must be a scalar Function"));
}

TEST_CASE("Physical components on a sheared parameter system scale by the "
          "metric diagonal",
          "[embedding][embedded_geometry][metric]")
{
  using pcms::ComponentScaling;
  const auto grid = pcms::test::MakeUnitSquareGrid(4);
  const auto embedded = pcms::CreateEmbeddedGeometry(
    grid, pcms::test::MakeSplinePosition(grid, Affine,
                                         pcms::csys::CreateCartesian(2)));
  const auto xi = embedded->GetParameterSystem();
  const auto bound = embedded->Bind(
    pcms::test::CreateDeviceCoordinateView({0.3, 0.6}, xi, 2).coordinate_view);
  // sqrt(g_ii) and sqrt(g^ii) of g = A^T A
  const auto g = InducedMetric(1.0);
  const Real det = g[0] * g[3] - g[1] * g[2];
  const std::array<Real, 2> vector_factors{std::sqrt(g[0]), std::sqrt(g[3])};
  const std::array<Real, 2> covector_factors{std::sqrt(g[3] / det),
                                             std::sqrt(g[0] / det)};
  const auto to_physical = pcms::detail::ToPhysical(xi, bound->Metric());
  const auto to_natural = pcms::detail::ToNatural(xi, bound->Metric());
  for (const bool covariant : {false, true}) {
    CAPTURE(covariant);
    const auto variance =
      covariant ? pcms::values::Covector : pcms::values::Vector;
    const auto& factors = covariant ? covector_factors : vector_factors;
    const pcms::ValueBasis natural{xi, ComponentScaling::Natural, variance};
    const pcms::ValueBasis physical{xi, ComponentScaling::Physical, variance};
    REQUIRE(to_physical.OutputBasis(natural) == physical);

    // natural (1, 1) to physical and back
    Kokkos::View<Real**, DeviceMemorySpace> values("values", 1, 2);
    Kokkos::deep_copy(values, 1.0);
    to_physical.Apply(pcms::ValueView<Real, DeviceMemorySpace>(
                        natural, pcms::MakeRank2View(values)),
                      pcms::ValueView<Real, DeviceMemorySpace>(
                        physical, pcms::MakeRank2View(values)));
    auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), values);
    CHECK_THAT(host(0, 0), WithinAbs(factors[0], 1e-12));
    CHECK_THAT(host(0, 1), WithinAbs(factors[1], 1e-12));
    to_natural.Apply(pcms::ValueView<Real, DeviceMemorySpace>(
                       physical, pcms::MakeRank2View(values)),
                     pcms::ValueView<Real, DeviceMemorySpace>(
                       natural, pcms::MakeRank2View(values)));
    Kokkos::deep_copy(host, values);
    CHECK_THAT(host(0, 0), WithinAbs(1.0, 1e-12));
    CHECK_THAT(host(0, 1), WithinAbs(1.0, 1e-12));
  }
}

TEST_CASE("Coordinates on an identified axis interpolate through the nearest "
          "image",
          "[embedding]")
{
  // u = (x, y + x / 2) on a grid periodic in y, with u2 stored modulo 1
  const std::vector<pcms::Identification> periodic = {{1, 0.0, 1.0}};
  const auto grid = pcms::test::MakeUnitSquareGrid(
    8, pcms::csys::CreateParameterSystem("xi", 2, periodic));
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(grid, 1, 2),
    [](const Real* xi, Real* u) {
      u[0] = xi[0];
      u[1] = std::fmod(xi[1] + 0.5 * xi[0], 1.0);
    },
    pcms::csys::CreateParameterSystem("u", 2, periodic));
  // across the last cell in y, and inside cells where u2 wraps
  const std::vector<Real> points = {0.3, 0.95, 0.6, 0.5, 0.95, 0.99};
  const auto xi = pcms::test::CreateDeviceCoordinateView(
    points, grid->GetCoordinateSystem(), 2);
  const auto bound =
    pcms::CreateCoordinateMap(position)->Bind(xi.coordinate_view);
  const auto mapped = pcms::test::CopyCoordinatesToHost(bound->Mapped());
  const auto j = bound->Jacobian();
  Kokkos::View<Real**, DeviceMemorySpace> jacobian("j", 3, 4);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 3),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int k = 0; k < 4; ++k) {
        jacobian(p, k) = j(p, k / 2, k % 2);
      }
    });
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), jacobian);
  for (int p = 0; p < 3; ++p) {
    CAPTURE(p);
    const Real x = points[std::size_t{2} * p];
    const Real y = points[std::size_t{2} * p + 1];
    CHECK_THAT(mapped(p, 0), WithinAbs(x, 1e-12));
    CHECK_THAT(std::remainder(mapped(p, 1) - y - 0.5 * x, 1.0),
               WithinAbs(0.0, 1e-12));
    const std::array<Real, 4> expected = {1.0, 0.0, 0.5, 1.0};
    for (int k = 0; k < 4; ++k) {
      CHECK_THAT(host(p, k), WithinAbs(expected[k], 1e-12));
    }
  }
  // the evaluation through the nearest image has no transpose
  const auto& space = position->GetSpace();
  const auto evaluator =
    space.CreatePointEvaluator<Real>(grid->Bind(xi.coordinate_view));
  const auto& basis = position->GetData().GetValueBasis();
  auto functional = space.CreateFunction<Real>("functional", basis);
  const Kokkos::View<Real**, DeviceMemorySpace> loads("loads", 3, 2);
  CHECK_THROWS_WITH(
    evaluator->EvaluateTranspose(
      pcms::ValueView<const Real, DeviceMemorySpace>(
        evaluator->OutputBasis(basis), pcms::MakeConstRank2View(loads)),
      functional),
    ContainsSubstring("nearest image"));
}

namespace
{

struct MappedPoints
{
  Kokkos::View<const Real**, pcms::HostMemorySpace> x;
  // jacobian(p, a * dim + b) = d x_a / d xi_b
  Kokkos::View<Real**, pcms::HostMemorySpace> jacobian;
};

// The image of `points` of `grid` under `position`, and its Jacobian
MappedPoints MapPoints(
  const std::shared_ptr<const pcms::Function<Real>>& position,
  const std::shared_ptr<const pcms::Geometry>& grid,
  const std::vector<Real>& points, int dim)
{
  const auto xi = pcms::test::CreateDeviceCoordinateView(
    points, grid->GetCoordinateSystem(), dim);
  const auto bound =
    pcms::CreateCoordinateMap(position)->Bind(xi.coordinate_view);
  const auto j = bound->Jacobian();
  const auto n = static_cast<pcms::LO>(j.extent(0));
  Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace> jacobian(
    "j", n, static_cast<std::size_t>(dim) * dim);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int k = 0; k < dim * dim; ++k) {
        jacobian(p, k) = j(p, k / dim, k % dim);
      }
    });
  return {
    pcms::test::CopyCoordinatesToHost(bound->Mapped()),
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), jacobian)};
}

// `count` random points with coordinate d in [lower[d], 1)
std::vector<Real> RandomPoints(int count, const std::vector<Real>& lower)
{
  std::mt19937 rng(7);
  std::vector<Real> points;
  for (int p = 0; p < count; ++p) {
    for (const Real l : lower) {
      points.push_back(std::uniform_real_distribution<Real>(l, 1.0)(rng));
    }
  }
  return points;
}

} // namespace

TEST_CASE("Coordinates across a sheared seam interpolate through the image "
          "nearest their neighbour along the seam axis",
          "[embedding]")
{
  // phi = P z + zShift(x, y), z wrapped; zShift jumps 3 P / 2 between y nodes
  constexpr int N = 6;
  constexpr Real P = 2 * std::numbers::pi / 5;
  pcms::UniformGrid<3> cells;
  cells.bot_left = {0.0, 0.0, 0.0};
  cells.edge_length = {1.0, 1.0, 1.0};
  cells.divisions = {N, N, N};
  const auto grid = pcms::uniform_grid::CreateGeometry(
    cells, pcms::csys::CreateParameterSystem("xi", 3, {{2, 0.0, 1.0}}));
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(grid, 1, 3),
    [](const Real* xi, Real* u) {
      u[0] = xi[0];
      u[1] = xi[1];
      u[2] = P * xi[2] + 0.3 * xi[0] + 1.5 * P * N * xi[1];
    },
    pcms::csys::CreateParameterSystem("u", 3, {{2, 0.0, P}}));
  auto points = RandomPoints(20, {0.0, 0.0, 0.0});
  const auto seam = RandomPoints(20, {0.0, 0.0, (N - 1.0) / N});
  points.insert(points.end(), seam.begin(), seam.end());
  const auto mapped = MapPoints(position, grid, points, 3);
  for (std::size_t p = 0; p < points.size() / 3; ++p) {
    CAPTURE(p);
    const Real* xi = &points[3 * p];
    CHECK_THAT(mapped.x(p, 2),
               WithinAbs(P * xi[2] + 0.3 * xi[0] + 1.5 * P * N * xi[1], 1e-9));
    CHECK_THAT(mapped.jacobian(p, 6), WithinAbs(0.3, 1e-9));
    CHECK_THAT(mapped.jacobian(p, 7), WithinAbs(1.5 * P * N, 1e-9));
    CHECK_THAT(mapped.jacobian(p, 8), WithinAbs(P, 1e-9));
  }
}

TEST_CASE("Longitudes across the seam of a periodic grid interpolate the "
          "unwrapped values",
          "[embedding]")
{
  // longitude = 360 i / N + 40 sin(j) in degrees, i wrapped
  constexpr int N = 8;
  const auto grid = pcms::test::MakeUnitSquareGrid(
    N, pcms::csys::CreateParameterSystem("xi", 2, {{0, 0.0, 1.0}}));
  const auto longitude = [](Real i, Real j) {
    return 360.0 * i / N + 40 * std::sin(j);
  };
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(grid, 1, 2),
    [&](const Real* xi, Real* u) {
      u[0] = longitude(N * xi[0], N * xi[1]);
      u[1] = 90 * xi[1];
    },
    pcms::csys::CreateParameterSystem("lonlat", 2, {{0, 0.0, 360.0}}));
  const auto points = RandomPoints(20, {(N - 1.0) / N, 0.0});
  const auto mapped = MapPoints(position, grid, points, 2);
  for (std::size_t p = 0; p < points.size() / 2; ++p) {
    CAPTURE(p);
    const Real tx = N * points[2 * p] - (N - 1);
    const Real y = N * points[2 * p + 1];
    const Real j = std::floor(y);
    const Real ty = y - j;
    const Real west =
      (1 - ty) * longitude(N - 1, j) + ty * longitude(N - 1, j + 1);
    const Real east = (1 - ty) * longitude(N, j) + ty * longitude(N, j + 1);
    CHECK_THAT(mapped.x(p, 0), WithinAbs((1 - tx) * west + tx * east, 1e-9));
    CHECK_THAT(mapped.jacobian(p, 0), WithinAbs(N * (east - west), 1e-9));
  }
}

TEST_CASE("Coordinates at a node wrapped along two axes interpolate through "
          "the image nearest each neighbour",
          "[embedding]")
{
  // u = (2 x + 2 y, x - y) on a grid periodic in x and y, both of period 1
  constexpr int N = 6;
  const std::vector<pcms::Identification> periodic = {{0, 0.0, 1.0},
                                                    {1, 0.0, 1.0}};
  const auto grid = pcms::test::MakeUnitSquareGrid(
    N, pcms::csys::CreateParameterSystem("xi", 2, periodic));
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(grid, 1, 2),
    [](const Real* xi, Real* u) {
      u[0] = 2 * xi[0] + 2 * xi[1];
      u[1] = xi[0] - xi[1];
    },
    pcms::csys::CreateParameterSystem("u", 2, periodic));
  const auto points = RandomPoints(20, {(N - 1.0) / N, (N - 1.0) / N});
  const auto mapped = MapPoints(position, grid, points, 2);
  for (std::size_t p = 0; p < points.size() / 2; ++p) {
    CAPTURE(p);
    const Real x = points[2 * p];
    const Real y = points[2 * p + 1];
    CHECK_THAT(mapped.x(p, 0), WithinAbs(2 * x + 2 * y, 1e-9));
    CHECK_THAT(mapped.x(p, 1), WithinAbs(x - y, 1e-9));
    const std::array<Real, 4> expected = {2.0, 2.0, 1.0, -1.0};
    for (int k = 0; k < 4; ++k) {
      CHECK_THAT(mapped.jacobian(p, k), WithinAbs(expected[k], 1e-9));
    }
  }
}

namespace
{

constexpr Real Pi = std::numbers::pi;
constexpr Real R0 = 3.0;

KOKKOS_INLINE_FUNCTION Real SafetyFactor(Real r)
{
  return 1 + r;
}

// Xi = (r, theta, z) -> U = (r, theta, z + q(r) theta)
struct Shear
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& x,
                                  Kokkos::Array<Real, 3>& u,
                                  Kokkos::Array<Real, 9>& j) const
  {
    u = {x[0], x[1], x[2] + SafetyFactor(x[0]) * x[1]};
    j = {1, 0, 0, 0, 1, 0, x[1], SafetyFactor(x[0]), 1};
  }
};

// U = (r, theta, phi) -> Xi = (r, theta, phi - q(r) theta)
struct Unshear
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& u,
                                  Kokkos::Array<Real, 3>& x,
                                  Kokkos::Array<Real, 9>& j) const
  {
    x = {u[0], u[1], u[2] - SafetyFactor(u[0]) * u[1]};
    j = {1, 0, 0, 0, 1, 0, -u[1], -SafetyFactor(u[0]), 1};
  }
};

// T = (a, b, c) -> Xi = (0.6 + 0.2 a, 1 + 2 b, 2 pi c)
struct Box
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& t,
                                  Kokkos::Array<Real, 3>& x,
                                  Kokkos::Array<Real, 9>& j) const
  {
    x = {0.6 + 0.2 * t[0], 1 + 2 * t[1], 2 * Pi * t[2]};
    j = {0.2, 0, 0, 0, 2, 0, 0, 0, 2 * Pi};
  }
};

std::shared_ptr<const pcms::Geometry> MakeGrid(
  std::shared_ptr<const pcms::CoordinateSystem> system,
  Kokkos::Array<Real, 3> lo, Kokkos::Array<Real, 3> length,
  Kokkos::Array<pcms::LO, 3> divisions)
{
  pcms::UniformGrid<3> grid;
  grid.bot_left = lo;
  grid.edge_length = length;
  grid.divisions = divisions;
  return pcms::uniform_grid::CreateGeometry(grid, std::move(system));
}

// A torus of radii 0.5 < r < 1 around R0: U = (r, theta, phi) -> (R, phi, Z)
// = (R0 + r cos theta, phi, r sin theta), a P1 position periodic in theta and
// phi, and the field-aligned chart Xi, periodic in z, linked into U.
struct Torus
{
  std::shared_ptr<const pcms::Geometry> u_grid = MakeGrid(
    pcms::csys::CreateParameterSystem("U", 3, {{1, 0, 2 * Pi}, {2, 0, 2 * Pi}}),
    {0.5, 0, 0}, {0.5, 2 * Pi, 2 * Pi}, {4, 64, 4});
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakePositionFunction(
      pcms::CreateLagrangeSpace(u_grid, 1, 3),
      [](const Real* u, Real* x) {
        x[0] = R0 + u[0] * std::cos(u[1]);
        x[1] = u[2];
        x[2] = u[0] * std::sin(u[1]);
      },
      pcms::csys::CreateCylindricalRThetaZ());
  std::shared_ptr<const pcms::EmbeddedGeometry> u =
    pcms::CreateEmbeddedGeometry(u_grid, position);
  std::shared_ptr<const pcms::Geometry> xi_grid =
    MakeGrid(pcms::csys::CreateParameterSystem("Xi", 3, {{2, 0, 2 * Pi}}),
             {0.5, 0, 0}, {0.5, 2 * Pi, 2 * Pi}, {4, 64, 8});
  std::shared_ptr<const pcms::CoordinateMap> shear = pcms::CreateCoordinateMap(
    xi_grid->GetCoordinateSystem(), u->GetParameterSystem(), Shear{});
  std::shared_ptr<const pcms::EmbeddedGeometry> xi =
    pcms::CreateEmbeddedGeometry(xi_grid, shear, u);

  // the integral of 1 on `g`
  static Real Volume(const std::shared_ptr<const pcms::EmbeddedGeometry>& g)
  {
    const auto space = pcms::CreateEmbeddedSpace(
      g, pcms::CreateLagrangeSpace(g->GetParameterGeometry(), 0, 1));
    auto one = space->CreateFunction<Real>("one");
    pcms::test::SetFieldComponents(
      one, [](Real, Real, Real, Real* out) { out[0] = 1; });
    return pcms::L2Integral(space).Apply(one);
  }

  // the parameter coordinates of `g` at the images of the U points `u_points`
  // under the position: the identity on a P1 field of the parameters
  std::vector<Real> RoundTrip(
    const std::shared_ptr<const pcms::EmbeddedGeometry>& g,
    const std::vector<Real>& u_points) const
  {
    const auto u_view = pcms::test::CreateDeviceCoordinateView(
      u_points, u->GetParameterSystem(), 3);
    const auto x =
      pcms::CreateCoordinateMap(position)->Bind(u_view.coordinate_view);
    const auto bound = g->Bind(x->MappedPoints());
    const auto space = pcms::CreateEmbeddedSpace(
      g, pcms::CreateLagrangeSpace(g->GetParameterGeometry(), 1, 3));
    auto f = space->CreateFunction<Real>("parameters",
                                         {.kind = pcms::ValueKind::Array});
    pcms::test::SetFieldComponents(f, [](Real a, Real b, Real c, Real* out) {
      out[0] = a;
      out[1] = b;
      out[2] = c;
    });
    const auto host =
      pcms::test::EvaluateToHost(*space->CreatePointEvaluator<Real>(bound), f,
                                 f.GetData().GetValueBasis(), 3);
    return std::vector<Real>(host.data(), host.data() + host.size());
  }
};

// Random Xi points off the last cell in z, and their U images.
std::pair<std::vector<Real>, std::vector<Real>> XiPoints(int n)
{
  std::mt19937 rng(7);
  std::uniform_real_distribution<Real> uniform(0.0, 1.0);
  std::vector<Real> xi, u;
  for (int i = 0; i < n; ++i) {
    const Real r = 0.55 + 0.4 * uniform(rng);
    const Real theta = 0.1 + (2 * Pi - 0.2) * uniform(rng);
    const Real z = 0.05 + 1.6 * Pi * uniform(rng);
    xi.insert(xi.end(), {r, theta, z});
    u.insert(u.end(), {r, theta, z + SafetyFactor(r) * theta});
  }
  return {xi, u};
}

} // namespace

TEST_CASE("A chained embedded geometry composes the factors of its links",
          "[embedding][embedded_geometry][chain]")
{
  const Torus t;
  CHECK(t.xi->GetCoordinateSystem() == t.u->GetCoordinateSystem());
  CHECK(t.xi->GetParent() == t.u);
  CHECK(t.xi->GetPosition() == nullptr);
  // the volume 2 pi^2 R0 (r1^2 - r0^2) to the error of the P1 position
  CHECK_THAT(Torus::Volume(t.xi), WithinRel(1.5 * R0 * Pi * Pi, 3e-3));
  // at theta at the middle of U's cells, where the P1 Jacobian is closest
  const std::vector<Real> params = {0.6, 21.5 * Pi / 32, 1.0,
                                    0.9, 50.5 * Pi / 32, 4.0};
  const auto bound = t.xi->Bind(pcms::test::CreateDeviceCoordinateView(
                                  params, t.xi->GetParameterSystem(), 3)
                                  .coordinate_view);
  const auto j = bound->Jacobian();
  const auto g = bound->Metric();
  const auto v = bound->VolumeFactor();
  Kokkos::View<Real**, DeviceMemorySpace> device("factors", 2, 19);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 2),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int k = 0; k < 9; ++k) {
        device(p, k) = j(p, k / 3, k % 3);
        device(p, 9 + k) = g(p, k / 3, k % 3);
      }
      device(p, 18) = v(p);
    });
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), device);
  for (std::size_t p = 0; p < 2; ++p) {
    CAPTURE(p);
    const Real r = params[3 * p];
    const Real theta = params[3 * p + 1];
    const Real q = SafetyFactor(r);
    const Real big_r = R0 + r * std::cos(theta);
    // J pushes Xi components forward to cylindrical (R, phi, Z) ones
    const Real jx[3][3] = {{std::cos(theta), -r * std::sin(theta), 0},
                           {theta, q, 1},
                           {std::sin(theta), r * std::cos(theta), 0}};
    const Real metric[3] = {1, big_r * big_r, 1};
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        Real gab = 0;
        for (int k = 0; k < 3; ++k) {
          gab += jx[k][a] * metric[k] * jx[k][b];
        }
        CHECK_THAT(host(p, 3 * a + b), WithinAbs(jx[a][b], 2e-3));
        CHECK_THAT(host(p, 9 + 3 * a + b), WithinAbs(gab, 5e-3 * (1 + gab)));
      }
    }
    CHECK_THAT(host(p, 18), WithinRel(big_r * r, 3e-3));
  }
}

TEST_CASE("The located points of a chained embedded geometry are the images "
          "through its parent",
          "[embedding][embedded_geometry][chain]")
{
  const Torus t;
  const std::vector<Real> params = {0.6, 1.0, 1.0, 0.9, 2.0, 4.0};
  const auto bound = t.xi->Bind(pcms::test::CreateDeviceCoordinateView(
                                  params, t.xi->GetParameterSystem(), 3)
                                  .coordinate_view);
  const auto located = bound->LocatedPoints();
  CHECK(located.GetCoordinateSystem() == t.xi->GetCoordinateSystem());
  const auto host = pcms::test::CopyCoordinatesToHost(located.GetValues());
  REQUIRE(host.extent(1) == 3);
  for (std::size_t p = 0; p < 2; ++p) {
    CAPTURE(p);
    const Real r = params[3 * p];
    const Real theta = params[3 * p + 1];
    const Real phi = params[3 * p + 2] + SafetyFactor(r) * theta;
    CHECK_THAT(host(p, 0), WithinAbs(R0 + r * std::cos(theta), 2e-3));
    CHECK_THAT(std::remainder(host(p, 1) - phi, 2 * Pi), WithinAbs(0.0, 1e-12));
    CHECK_THAT(host(p, 2), WithinAbs(r * std::sin(theta), 2e-3));
  }
}

TEST_CASE("A chained embedded geometry binds physical points through its "
          "parent",
          "[embedding][embedded_geometry][chain]")
{
  const Torus t;
  const auto unshear = pcms::CreateCoordinateMap(
    t.u->GetParameterSystem(), t.xi->GetParameterSystem(), Unshear{});
  const auto inverted =
    pcms::CreateEmbeddedGeometry(t.xi_grid, t.shear, t.u, {.inverse = unshear});
  CHECK_THROWS_WITH(
    pcms::CreateEmbeddedGeometry(t.xi_grid, t.shear, t.u, {.inverse = t.shear}),
    ContainsSubstring("the inverse must map"));
  CHECK_THROWS_WITH(pcms::CreateEmbeddedGeometry(t.xi_grid, unshear, t.u),
                    ContainsSubstring("the link must map"));
  const auto [xi, u] = XiPoints(20);
  for (const auto& g : {t.xi, inverted}) {
    CAPTURE(g == inverted);
    const auto found = t.RoundTrip(g, u);
    for (std::size_t i = 0; i < xi.size(); ++i) {
      CAPTURE(i);
      CHECK_THAT(found[i], WithinAbs(xi[i], 1e-8));
    }
  }
}

TEST_CASE("A chained embedded geometry uses a given volume factor",
          "[embedding][embedded_geometry][chain]")
{
  const Torus t;
  auto volume = std::make_shared<pcms::Function<Real>>(
    pcms::CreateLagrangeSpace(t.xi_grid, 0, 1)->CreateFunction<Real>("volume"));
  pcms::test::SetFieldComponents(
    *volume, [](Real, Real, Real, Real* out) { out[0] = 2; });
  const auto xi =
    pcms::CreateEmbeddedGeometry(t.xi_grid, t.shear, t.u, {.volume = volume});
  // 2 times the parameter volume 0.5 (2 pi)^2
  CHECK_THAT(Torus::Volume(xi), WithinRel(4 * Pi * Pi, 1e-12));
}

TEST_CASE("Three chained embedded geometries compose",
          "[embedding][embedded_geometry][chain]")
{
  const Torus t;
  const auto box = pcms::csys::CreateParameterSystem("T", 3);
  const auto tube = pcms::CreateEmbeddedGeometry(
    MakeGrid(box, {0, 0, 0}, {1, 1, 1}, {2, 16, 2}),
    pcms::CreateCoordinateMap(box, t.xi->GetParameterSystem(), Box{}), t.xi);
  CHECK(tube->GetCoordinateSystem() == t.u->GetCoordinateSystem());
  // 2 pi (R0 (r1^2 - r0^2) / 2 (theta1 - theta0) + (r1^3 - r0^3) / 3
  // (sin theta1 - sin theta0)) over 0.6 < r < 0.8, 1 < theta < 3
  const Real volume =
    2 * Pi *
    (R0 * 0.28 + (0.512 - 0.216) / 3 * (std::sin(3.0) - std::sin(1.0)));
  CHECK_THAT(Torus::Volume(tube), WithinRel(volume, 3e-3));
  const std::vector<Real> params = {0.3, 0.6, 0.4, 0.8, 0.1, 0.7};
  std::vector<Real> u;
  for (std::size_t i = 0; i < params.size(); i += 3) {
    const Real r = 0.6 + 0.2 * params[i];
    const Real theta = 1 + 2 * params[i + 1];
    u.insert(u.end(),
             {r, theta, 2 * Pi * params[i + 2] + SafetyFactor(r) * theta});
  }
  const auto found = t.RoundTrip(tube, u);
  for (std::size_t i = 0; i < params.size(); ++i) {
    CHECK_THAT(found[i], WithinAbs(params[i], 1e-8));
  }
}

namespace
{

// The affine map Affine into a parameter system, which has no metric.
struct MetricLess
{
  std::shared_ptr<const pcms::Geometry> grid =
    pcms::test::MakeUnitSquareGrid(4);
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakePositionFunction(
      pcms::CreateLagrangeSpace(grid, 1, 2),
      [](const Real* xi, Real* x) {
        const auto y = Affine(xi[0], xi[1]);
        x[0] = y[0];
        x[1] = y[1];
      },
      pcms::csys::CreateParameterSystem("s", 2));

  // the integral of 1 on `g`
  static Real Volume(const std::shared_ptr<const pcms::Geometry>& g,
                     const std::shared_ptr<const pcms::Geometry>& parameter)
  {
    const auto plain = pcms::CreateLagrangeSpace(parameter, 0, 1);
    const auto embedded =
      std::dynamic_pointer_cast<const pcms::EmbeddedGeometry>(g);
    const auto space =
      embedded != nullptr ? pcms::CreateEmbeddedSpace(embedded, plain) : plain;
    auto one = space->CreateFunction<Real>("one");
    pcms::test::SetFieldComponents(one,
                                   [](Real, Real, Real* out) { out[0] = 1; });
    return pcms::L2Integral(space).Apply(one);
  }

  [[nodiscard]] std::shared_ptr<const pcms::BoundPoints> Bind(
    const pcms::Geometry& g) const
  {
    return g.Bind(pcms::test::CreateDeviceCoordinateView(
                    Params, grid->GetCoordinateSystem(), 2)
                    .coordinate_view);
  }
};

} // namespace

TEST_CASE("Without a metric, an embedded geometry measures by a given volume "
          "factor",
          "[embedding][embedded_geometry][measure]")
{
  const MetricLess m;
  auto volume = std::make_shared<pcms::Function<Real>>(
    pcms::CreateLagrangeSpace(m.grid, 0, 1)->CreateFunction<Real>("volume"));
  pcms::test::SetFieldComponents(*volume,
                                 [](Real, Real, Real* out) { out[0] = 3.5; });
  const auto g =
    pcms::CreateEmbeddedGeometry(m.grid, m.position, {.volume = volume});
  CHECK_THAT(MetricLess::Volume(g, m.grid), WithinRel(3.5, 1e-12));
  CHECK_THROWS_WITH(m.Bind(*g)->Metric(), ContainsSubstring("no metric"));
}

TEST_CASE("Without a metric in its coordinate system, an embedded geometry "
          "uses a given metric",
          "[embedding][embedded_geometry][measure]")
{
  const MetricLess m;
  const std::array<Real, 4> given = {2.0, 0.5, 0.5, 3.0};
  const auto g = pcms::CreateEmbeddedGeometry(
    m.grid, m.position, {.metric = ConstantMetric(m.grid, given, false)});
  const auto factors = FactorsOnHost(*m.Bind(*g));
  for (int p = 0; p < 4; ++p) {
    CAPTURE(p);
    for (int k = 0; k < 4; ++k) {
      CHECK_THAT(factors.metric(p, k), WithinAbs(given[k], 1e-12));
    }
    CHECK_THAT(factors.volume(p), WithinAbs(std::sqrt(5.75), 1e-12));
  }
}

TEST_CASE("Without a metric, an embedded geometry of a square Jacobian "
          "measures by |det J|",
          "[embedding][embedded_geometry][measure]")
{
  const MetricLess m;
  const auto g = pcms::CreateEmbeddedGeometry(m.grid, m.position);
  const Real det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  CHECK_THAT(MetricLess::Volume(g, m.grid), WithinRel(det, 1e-12));
  CHECK_THROWS_WITH(m.Bind(*g)->Metric(), ContainsSubstring("no metric"));
}

TEST_CASE("A surface in a coordinate system without a metric has no measure",
          "[embedding][embedded_geometry][measure]")
{
  const MetricLess m;
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(m.grid, 1, 3),
    [](const Real* xi, Real* x) {
      x[0] = xi[0];
      x[1] = xi[1];
      x[2] = xi[0] * xi[1];
    },
    pcms::csys::CreateParameterSystem("s3", 3));
  CHECK_THROWS_WITH(pcms::CreateEmbeddedGeometry(m.grid, position),
                    ContainsSubstring("no measure"));
}

TEST_CASE("A geometry in a parameter system measures by its coordinates",
          "[embedding][measure]")
{
  const MetricLess m;
  CHECK_THAT(MetricLess::Volume(m.grid, m.grid), WithinRel(1.0, 1e-12));
  const auto bound = m.Bind(*m.grid);
  const auto volume = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(),
    Kokkos::View<const Real*, DeviceMemorySpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      bound->VolumeFactor().data_handle(), 1));
  CHECK(volume(0) == 1.0);
  CHECK_THROWS_WITH(bound->Metric(), ContainsSubstring("no metric"));
}

namespace
{

// (a, b, c) -> (a, b)
struct DropLast
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& u,
                                  Kokkos::Array<Real, 2>& xi,
                                  Kokkos::Array<Real, 6>& j) const
  {
    xi = {u[0], u[1]};
    j = {1, 0, 0, 0, 1, 0};
  }
};

} // namespace

TEST_CASE("An embedding does not collapse parameter directions",
          "[embedding][embedded_geometry]")
{
  const auto box = pcms::csys::CreateParameterSystem("B", 3);
  const auto grid = MakeGrid(box, {0, 0, 0}, {1, 1, 1}, {2, 2, 2});
  const auto flat = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(grid, 1, 2),
    [](const Real* u, Real* x) {
      x[0] = 2 + u[0];
      x[1] = u[1] + u[2];
    },
    pcms::csys::CreateCylindricalRZ());
  CHECK_THROWS_WITH(pcms::CreateEmbeddedGeometry(grid, flat),
                    ContainsSubstring("dim(target) >= dim(parameters)"));
  // a surface in 3D
  const auto square = pcms::test::MakeUnitSquareGrid(4);
  const auto surface = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(square, 1, 3),
    [](const Real* xi, Real* x) {
      x[0] = xi[0];
      x[1] = xi[1];
      x[2] = 0.5 * xi[0] + 0.25 * xi[1];
    },
    pcms::csys::CreateCartesian(3));
  CHECK_NOTHROW(pcms::CreateEmbeddedGeometry(square, surface));
  // a link from 3 parameters into the 2 of an annulus
  const EmbeddedFixture f;
  const auto lower = pcms::CreateCoordinateMap<3, 2>(
    box, f.embedded->GetParameterSystem(), DropLast{});
  CHECK_THROWS_WITH(pcms::CreateEmbeddedGeometry(grid, lower, f.embedded),
                    ContainsSubstring("dim(target) >= dim(parameters)"));
}
