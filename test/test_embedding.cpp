#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "embedding_test_utils.hpp"
#include <numbers>
#include <random>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
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
  CHECK(bound->ValueTransform().OutputBasis(coordinates) == coordinates);
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

TEST_CASE("Parameter-system vectors and covectors come out in the coordinate "
          "system of an "
          "embedded geometry",
          "[embedding][embedded_geometry]")
{
  const EmbeddedFixture f;
  const auto grid_space = pcms::CreateLagrangeSpace(f.grid, 1, 2);
  const auto space = pcms::CreateEmbeddedSpace(f.embedded, grid_space);
  const auto bound =
    f.embedded->Bind(pcms::test::CreateDeviceCoordinateView(
                       Params, f.embedded->GetParameterSystem(), 2)
                       .coordinate_view);
  for (const auto variance : {pcms::values::Vector, pcms::values::Covector}) {
    // the constant field e_xi1 (or dxi1)
    // components on the parameter system (the space's default is the coordinate
    // system of the embedded geometry)
    auto field = space->CreateFunction<Real>("e1", variance,
                                             f.embedded->GetParameterSystem());
    pcms::test::SetFieldComponents(field, [](Real, Real, Real* out) {
      out[0] = 1.0;
      out[1] = 0.0;
    });
    const auto evaluator = space->CreatePointEvaluator<Real>(bound);
    const auto basis = evaluator->OutputBasis(field.GetData().GetValueBasis());
    CHECK(basis.system->Kind() == "cylindrical_rz");
    CHECK(basis.scaling == pcms::ComponentScaling::Physical);
    Kokkos::View<Real**, DeviceMemorySpace> out("out", 4, 2);
    evaluator->Evaluate(field, pcms::ValueView<Real, DeviceMemorySpace>(
                                 basis, pcms::MakeRank2View(out)));
    const auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
    const bool covariant = variance == pcms::values::Covector;
    for (int p = 0; p < 4; ++p) {
      const Real theta = Annulus::Theta(Params[2 * p + 1]);
      // the vector e_xi1 is d x / d xi1 = 0.5 (cos, sin); the covector
      // dxi1 has the physical components 2 (cos, sin)
      const Real scale = covariant ? 2.0 : 0.5;
      CAPTURE(p, covariant);
      CHECK_THAT(host(p, 0), WithinAbs(scale * std::cos(theta), 1e-4));
      CHECK_THAT(host(p, 1), WithinAbs(scale * std::sin(theta), 1e-4));
    }
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
  // a scalar field of two components, and vectors and covectors in
  // parameter-system components, which come out in the coordinate system's
  for (const auto variance :
       {pcms::values::Scalar, pcms::values::Vector, pcms::values::Covector}) {
    auto basis_in = pcms::detail::MakeValueBasis(
      space->GetCoordinateSystem(), variance, f.embedded->GetParameterSystem(),
      std::nullopt);
    if (variance.Rank() == 0) {
      basis_in.kind = pcms::ValueKind::Array;
    }
    auto x = space->CreateFunction<Real>("x", basis_in);
    auto functional = space->CreateFunction<Real>("functional", basis_in);
    const auto basis = evaluator->OutputBasis(x.GetData().GetValueBasis());
    CHECK((variance == pcms::values::Scalar ||
           basis.system->Kind() == "cylindrical_rz"));
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
  const auto quadrature = pcms::CreateGaussQuadrature(f.embedded, 4);
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
