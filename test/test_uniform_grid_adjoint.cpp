#include "periodic_test_utils.hpp"
#include "embedding_test_utils.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/embedding.hpp>
#include <pcms/geometry/product_geometry.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include <cmath>
#include <numbers>

namespace
{

using pcms::DeviceMemorySpace;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::test::GetDofs;
using pcms::test::Pseudo;
using pcms::test::SetDofs;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

std::shared_ptr<const pcms::Geometry> Grid2(
  std::shared_ptr<const pcms::CoordinateSystem> system =
    pcms::csys::CreateCartesian(2),
  int nx = 5, int ny = 4)
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {-1.0, 0.5};
  grid.edge_length = {3.0, 2.0};
  grid.divisions = {nx, ny};
  return pcms::uniform_grid::CreateGeometry(grid, std::move(system));
}

std::shared_ptr<const pcms::Geometry> Grid3(
  std::shared_ptr<const pcms::CoordinateSystem> system =
    pcms::csys::CreateCartesian(3))
{
  pcms::UniformGrid<3> grid;
  grid.bot_left = {-1.0, 0.5, 0.0};
  grid.edge_length = {3.0, 2.0, 1.5};
  grid.divisions = {4, 3, 5};
  return pcms::uniform_grid::CreateGeometry(grid, std::move(system));
}

/// @brief n points spread over the box [lo, lo + length).
std::vector<std::vector<Real>> Points(int n, const std::vector<Real>& lo,
                                      const std::vector<Real>& length)
{
  std::vector<std::vector<Real>> points;
  for (int i = 0; i < n; ++i) {
    std::vector<Real> x;
    for (std::size_t d = 0; d < lo.size(); ++d) {
      x.push_back(
        lo[d] +
        length[d] *
          std::fmod(0.1234 + 0.6180339 * (i + 1) * (d + 1) + 0.31 * d, 1.0));
    }
    points.push_back(x);
  }
  return points;
}

std::shared_ptr<const pcms::BoundPoints> BindPoints(
  const pcms::Geometry& geometry, const std::vector<std::vector<Real>>& points,
  pcms::OutOfBoundsPolicy policy = {})
{
  const auto device = pcms::test::DevicePoints(points);
  return geometry.Bind(
    pcms::CoordinateView<DeviceMemorySpace>(geometry.GetCoordinateSystem(),
                                            pcms::MakeConstRank2View(device)),
    policy);
}

Kokkos::View<Real**, DeviceMemorySpace> DeviceValues(
  const std::vector<Real>& values, std::size_t columns)
{
  Kokkos::View<Real**, HostMemorySpace> host("values", values.size() / columns,
                                             columns);
  for (std::size_t i = 0; i < values.size(); ++i) {
    host(i / columns, i % columns) = values[i];
  }
  return pcms::detail::ToDeviceView(host);
}

std::vector<Real> HostValues(const Kokkos::View<Real**, DeviceMemorySpace>& v)
{
  const auto host = Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, v);
  std::vector<Real> out;
  for (std::size_t i = 0; i < host.extent(0); ++i) {
    for (std::size_t j = 0; j < host.extent(1); ++j) {
      out.push_back(host(i, j));
    }
  }
  return out;
}

std::vector<Real> AllDofs(const pcms::Field<Real>& field)
{
  const auto values = field.GetDOFHolderDataHost().GetValues();
  std::vector<Real> out;
  for (std::size_t i = 0; i < values.extent(0); ++i) {
    for (std::size_t j = 0; j < values.extent(1); ++j) {
      out.push_back(values(i, j));
    }
  }
  return out;
}

Real Dot(const std::vector<Real>& a, const std::vector<Real>& b)
{
  REQUIRE(a.size() == b.size());
  Real sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

/// @return E^T y for the evaluator of `space` at `bound`
std::vector<Real> Transpose(const pcms::FunctionSpace& space,
                            std::shared_ptr<const pcms::BoundPoints> bound,
                            const std::vector<Real>& y)
{
  const auto evaluator = space.CreatePointEvaluator<Real>(std::move(bound));
  auto functional = space.CreateFunction<Real>(
    "functional", pcms::test::ComponentsBasis(space));
  const auto device_y = DeviceValues(
    y, static_cast<std::size_t>(space.GetLayout()->GetNumComponents()));
  evaluator->EvaluateTranspose(
    pcms::ValueView<const Real, DeviceMemorySpace>(
      functional.GetData().GetValueBasis(), pcms::MakeConstRank2View(device_y)),
    functional);
  return AllDofs(functional);
}

/// @return E x for the evaluator of `space` at `bound`
std::vector<Real> Apply(const pcms::FunctionSpace& space,
                        std::shared_ptr<const pcms::BoundPoints> bound,
                        const std::vector<Real>& x)
{
  const auto evaluator = space.CreatePointEvaluator<Real>(bound);
  const int components = space.GetLayout()->GetNumComponents();
  auto u = space.CreateFunction<Real>("u", pcms::test::ComponentsBasis(space));
  const auto device_x = DeviceValues(x, components);
  u.SetDOFHolderData(pcms::ValueView<const Real, DeviceMemorySpace>(
    u.GetData().GetValueBasis(), pcms::MakeConstRank2View(device_x)));
  Kokkos::View<Real**, DeviceMemorySpace> out("Ex", bound->NumPoints(),
                                              components);
  evaluator->Evaluate(u, pcms::test::TagLike(u, out));
  return HostValues(out);
}

/// @return u = (x + 2 y) (1 + 0.5 cos z - 0.3 sin 2 z) on `space`, the
///         product of `plane` and the Fourier modes {0, 1, 2}
pcms::Function<Real> ProductFunction(const pcms::FunctionSpace& plane,
                                     const pcms::FunctionSpace& space)
{
  std::vector<Real> c;
  const auto xy = pcms::test::CopyCoordinatesToHost(
    plane.GetLayout()->GetDOFHolderCoordinates().GetValues());
  for (std::size_t i = 0; i < xy.extent(0); ++i) {
    for (const Real f : {1.0, 0.5, 0.0, 0.0, -0.3}) {
      c.push_back((xy(i, 0) + 2 * xy(i, 1)) * f);
    }
  }
  auto u = space.CreateFunction<Real>("u", pcms::test::ComponentsBasis(space));
  SetDofs(u, c);
  return u;
}

} // namespace

TEST_CASE("Uniform-grid Lagrange spaces: <E x, y> = <x, E^T y> and the "
          "matrix is E",
          "[uniform_grid][adjoint]")
{
  for (const int order : {0, 1}) {
    CAPTURE(order);
    pcms::test::CheckAdjointAndMatrix(
      pcms::CreateLagrangeSpace(Grid2(), order, 1),
      Points(23, {-1.0, 0.5}, {3.0, 2.0}));
    pcms::test::CheckAdjointAndMatrix(
      pcms::CreateLagrangeSpace(Grid3(), order, 1),
      Points(41, {-1.0, 0.5, 0.0}, {3.0, 2.0, 1.5}));
  }
}

TEST_CASE("Uniform-grid Lagrange spaces: the transpose of several "
          "components, and of points outside the grid",
          "[uniform_grid][adjoint]")
{
  // the last two points are outside the grid
  auto points = Points(17, {-1.0, 0.5}, {3.0, 2.0});
  points.push_back({-1.5, 0.7});
  points.push_back({1.0, 3.1});
  for (const int order : {0, 1}) {
    const auto space = pcms::CreateLagrangeSpace(Grid2(), order, 2);
    const auto num_values =
      static_cast<std::size_t>(space->GetLayout()->GetNumDofHolder()) * 2;
    const auto x = Pseudo(num_values, 3);
    const auto y = Pseudo(points.size() * 2, 4);
    for (const auto mode : {pcms::OutOfBoundsMode::FILL,
                            pcms::OutOfBoundsMode::NEAREST_BOUNDARY}) {
      CAPTURE(order, static_cast<int>(mode));
      const auto bound =
        BindPoints(*space->GetGeometry(), points, {.mode = mode});
      auto ex = Apply(*space, bound, x);
      if (mode == pcms::OutOfBoundsMode::FILL) {
        for (std::size_t k = ex.size() - 4; k < ex.size(); ++k) {
          ex[k] = 0;
        }
      }
      const auto ety = Transpose(*space, bound, y);
      CHECK(Dot(ex, y) ==
            Catch::Approx(Dot(x, ety)).epsilon(1e-12).margin(1e-12));
    }
  }
}

TEST_CASE("Uniform-grid Lagrange spaces: the matrix of the gradient is that "
          "of the partials",
          "[uniform_grid][adjoint]")
{
  const auto space = pcms::CreateLagrangeSpace(Grid3(), 1, 1);
  const auto points = Points(13, {-1.0, 0.5, 0.0}, {3.0, 2.0, 1.5});
  const auto bound = BindPoints(*space->GetGeometry(), points);
  const auto evaluator = space->CreatePointEvaluator<Real>(
    bound, pcms::DifferentialOperator::Gradient);
  REQUIRE(evaluator->SupportsMatrix());
  const auto x =
    Pseudo(static_cast<std::size_t>(space->GetLayout()->GetNumDofHolder()), 5);
  auto u =
    space->CreateFunction<Real>("u", pcms::test::ComponentsBasis(*space));
  SetDofs(u, x);
  const auto gradient = pcms::test::EvaluateAt(
    *space, u, points, pcms::DifferentialOperator::Gradient);
  const auto matrix = evaluator->GetMatrix();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.row_offsets);
  const auto columns =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.columns);
  const auto values =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.values);
  REQUIRE(offsets.extent(0) == points.size() * 3 + 1);
  for (std::size_t p = 0; p < points.size(); ++p) {
    for (int j = 0; j < 3; ++j) {
      Real sum = 0;
      for (pcms::LO k = offsets(p * 3 + j); k < offsets(p * 3 + j + 1); ++k) {
        sum += values(k) * x[columns(k)];
      }
      CHECK(sum == Catch::Approx(gradient[p][j]).epsilon(1e-12).margin(1e-12));
    }
  }
}

TEST_CASE("Product of a uniform grid and a periodic axis: transpose and "
          "gradient, with the loads of wrapped points on the identified DOFs",
          "[uniform_grid][product][adjoint]")
{
  const auto plane = pcms::CreateLagrangeSpace(Grid2(), 1, 1);
  const auto zs =
    pcms::CreateFourierSpace(pcms::test::PeriodicGrid(5), {0, 1, 2});
  const auto space = pcms::CreateProductSpace(plane, zs);
  auto points = Points(31, {-1.0, 0.5, -10.0}, {3.0, 2.0, 20.0});
  pcms::test::CheckAdjointAndMatrix(space, points);

  const auto y = Pseudo(points.size(), 6);
  const auto at =
    Transpose(*space, BindPoints(*space->GetGeometry(), points), y);
  for (auto& p : points) {
    p[2] += 3 * TwoPi;
  }
  const auto wrapped =
    Transpose(*space, BindPoints(*space->GetGeometry(), points), y);
  for (std::size_t i = 0; i < at.size(); ++i) {
    CHECK(wrapped[i] == Catch::Approx(at[i]).margin(1e-12));
  }

  const auto u = ProductFunction(*plane, *space);
  const auto pts = Points(7, {-1.0, 0.5, -3.0}, {3.0, 2.0, 6.0});
  const auto gradient = pcms::test::EvaluateAt(
    *space, u, pts, pcms::DifferentialOperator::Gradient);
  for (std::size_t p = 0; p < pts.size(); ++p) {
    const Real x = pts[p][0], yy = pts[p][1], z = pts[p][2];
    const Real f = 1 + 0.5 * std::cos(z) - 0.3 * std::sin(2 * z);
    const Real df = -0.5 * std::sin(z) - 0.6 * std::cos(2 * z);
    CHECK(gradient[p][0] == Catch::Approx(f).margin(1e-12));
    CHECK(gradient[p][1] == Catch::Approx(2 * f).margin(1e-12));
    CHECK(gradient[p][2] == Catch::Approx((x + 2 * yy) * df).margin(1e-12));
  }
}

namespace
{

/// @brief A multilinear BOUT++-like position (R, z + shift(x, y), Z) of the
///        index coordinates xi.
const Real Wedge = TwoPi / 5;

void Shifted(const Real* xi, Real* x)
{
  x[0] = 3.5 + 0.4 * xi[0] + 0.1 * xi[1] + 0.05 * xi[0] * xi[1];
  x[1] = xi[2] + 0.5 * xi[0] + 0.3 * xi[1] + 0.4 * xi[0] * xi[1];
  x[2] = 0.1 * xi[0] + 0.3 * xi[1] + 0.05 * xi[0] * xi[1];
}

// constant, so that every grid represents it exactly
constexpr Real GInv[9] = {0.25, 0.1, 0.0, 0.1, 4.0, 0.0, 0.0, 0.0, 1.0};

// the parameter system the grids of the run share
const std::shared_ptr<const pcms::CoordinateSystem>& BoutIndex()
{
  static const auto system = pcms::csys::CreateParameterSystem("bout_index", 3);
  return system;
}

struct Bout
{
  std::shared_ptr<const pcms::CoordinateSystem> logical = BoutIndex();
  std::shared_ptr<const pcms::Geometry> grid;
  std::shared_ptr<const pcms::EmbeddedGeometry> geometry;
  std::shared_ptr<const pcms::FunctionSpace> fields;

  explicit Bout(int n)
  {
    pcms::UniformGrid<3> g;
    g.bot_left = {0.0, 0.0, 0.0};
    g.edge_length = {1.0, 1.0, Wedge};
    g.divisions = {n, n, 2 * n};
    grid = pcms::uniform_grid::CreateGeometry(g, logical);
    const auto space = pcms::CreateLagrangeSpace(grid, 1, 1);
    const auto position = pcms::test::MakePositionFunction(
      pcms::CreateLagrangeSpace(grid, 1, 3), Shifted,
      pcms::csys::CreateCylindricalRThetaZ(Wedge));
    auto g_inv = std::make_shared<pcms::Function<Real>>(
      pcms::CreateLagrangeSpace(grid, 1, 9)
        ->CreateFunction<Real>(
          "g^ij", pcms::ValueBasis{logical, pcms::ComponentScaling::Natural,
                                   pcms::values::Tensor}));
    pcms::test::SetFieldComponents(*g_inv, [](Real, Real, Real, Real* out) {
      for (int k = 0; k < 9; ++k) {
        out[k] = GInv[k];
      }
    });
    geometry = pcms::CreateEmbeddedGeometry(grid, position, {.metric = g_inv});
    fields = pcms::CreateEmbeddedSpace(geometry, space);
  }
};

} // namespace

TEST_CASE("BOUT++-like embedded grid: the transpose at points translated by "
          "the symmetry loads the same DOFs, and filled points load none",
          "[uniform_grid][embedding][adjoint]")
{
  const Bout bout(3);
  const std::vector<Real> parameters = {0.3, 0.6, 0.5, 0.7, 0.2, 0.2,
                                        0.9, 0.9, 1.1, 0.5, 0.5, 0.05};
  const auto num_dofs =
    static_cast<std::size_t>(bout.fields->GetLayout()->GetNumDofHolder());
  const auto x = Pseudo(num_dofs, 9);
  const auto y = Pseudo(5, 10);
  std::vector<Real> reference;
  for (const int periods : {0, 1, -3, 7}) {
    CAPTURE(periods);
    std::vector<Real> points;
    for (std::size_t i = 0; i < parameters.size(); i += 3) {
      Real p[3];
      Shifted(&parameters[i], p);
      points.insert(points.end(), {p[0], p[1] + periods * Wedge, p[2]});
    }
    points.insert(points.end(), {10.0, 0.0, 10.0});
    const auto bound = pcms::test::Bind(
      *bout.geometry, points,
      pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
    REQUIRE(bound->MissingPoints().extent(0) == 1);
    auto ex = Apply(*bout.fields, bound, x);
    CHECK(ex.back() == -7.0);
    ex.back() = 0;
    const auto ety = Transpose(*bout.fields, bound, y);
    CHECK(Dot(ex, y) ==
          Catch::Approx(Dot(x, ety)).epsilon(1e-12).margin(1e-12));
    if (reference.empty()) {
      reference = ety;
    }
    for (std::size_t i = 0; i < num_dofs; ++i) {
      CHECK(ety[i] == Catch::Approx(reference[i]).margin(1e-9));
    }
  }
}

namespace
{

/// @return the integral of x^a over [lo, hi]
Real MonomialIntegral(int a, Real lo, Real hi)
{
  return (std::pow(hi, a + 1) - std::pow(lo, a + 1)) / (a + 1);
}

/// @return the sum of the weights times f at the points of `quadrature`
template <typename F>
Real Integrate(const pcms::Quadrature& quadrature, F f)
{
  const auto x = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  const auto w =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, quadrature.weights);
  Real sum = 0;
  for (std::size_t p = 0; p < w.extent(0); ++p) {
    sum += w(p) * f(&x(p, 0));
  }
  return sum;
}

/// @return sqrt(det g_ij) of GInv
Real VolumeFactor()
{
  const Real det = GInv[0] * (GInv[4] * GInv[8] - GInv[5] * GInv[7]) -
                   GInv[1] * (GInv[3] * GInv[8] - GInv[5] * GInv[6]) +
                   GInv[2] * (GInv[3] * GInv[7] - GInv[4] * GInv[6]);
  return 1 / std::sqrt(det);
}

} // namespace

TEST_CASE("Uniform-grid quadrature is exact for polynomials of its degree in "
          "each coordinate",
          "[uniform_grid][quadrature]")
{
  const std::vector<Real> lo = {-1.0, 0.5, 0.0};
  const std::vector<Real> hi = {2.0, 2.5, 1.5};
  for (const auto& system : {pcms::csys::CreateCartesian(3),
                             pcms::csys::CreateParameterSystem("xi", 3)}) {
    const auto grid = Grid3(system);
    for (int order = 1; order <= 5; ++order) {
      CAPTURE(order, system->Kind());
      const auto quadrature = pcms::CreateGaussQuadrature(grid, order);
      const int per_axis = order / 2 + 1;
      CHECK(quadrature.points->NumPoints() ==
            4 * 3 * 5 * per_axis * per_axis * per_axis);
      CHECK(quadrature.points->GetGeometry() == grid);
      for (int a = 0; a <= order; ++a) {
        for (int c = 0; c <= order; c += 2) {
          const int b = order - a;
          const Real expected = MonomialIntegral(a, lo[0], hi[0]) *
                                MonomialIntegral(b, lo[1], hi[1]) *
                                MonomialIntegral(c, lo[2], hi[2]);
          const Real sum = Integrate(quadrature, [&](const Real* x) {
            return std::pow(x[0], a) * std::pow(x[1], b) * std::pow(x[2], c);
          });
          CHECK(sum == Catch::Approx(expected).epsilon(1e-12).margin(1e-12));
        }
      }
    }
  }
  const auto rz = Grid2(pcms::csys::CreateCylindricalRZ());
  const auto quadrature = pcms::CreateGaussQuadrature(rz, 3);
  const Real sum = Integrate(
    quadrature, [](const Real* x) { return x[0] * x[0] * x[1] * x[1] * x[1]; });
  CHECK(sum == Catch::Approx(MonomialIntegral(3, -1.0, 2.0) *
                             MonomialIntegral(3, 0.5, 2.5))
                 .epsilon(1e-12));
  CHECK_THROWS(pcms::CreateMonteCarloQuadrature(rz, 2));
  CHECK_THROWS(pcms::CreateGaussQuadrature(rz, -1));
  CHECK(pcms::CreateLagrangeSpace(Grid2(), 0, 1)
          ->CreateQuadrature()
          .points->NumPoints() == 20);
  CHECK(pcms::CreateLagrangeSpace(Grid2(), 1, 1)
          ->CreateQuadrature()
          .points->NumPoints() == 80);
}

TEST_CASE("Product of a uniform grid and a periodic axis: the tensor "
          "quadrature",
          "[uniform_grid][product][quadrature]")
{
  const auto plane = pcms::CreateLagrangeSpace(Grid2(), 1, 1);
  const auto space = pcms::CreateProductSpace(
    plane, pcms::CreateFourierSpace(pcms::test::PeriodicGrid(5), {0, 1, 2}));
  const auto u = ProductFunction(*plane, *space);
  const Real expected = (MonomialIntegral(1, -1.0, 2.0) * 2.0 +
                         2 * 3.0 * MonomialIntegral(1, 0.5, 2.5)) *
                        TwoPi;
  CHECK(pcms::L2Integral(space).Apply(u) ==
        Catch::Approx(expected).epsilon(1e-12));
}

TEST_CASE("Conservative projection into a uniform grid keeps the integral",
          "[uniform_grid][projection]")
{
  // nested cells: the target rules integrate source x target basis exactly
  const auto system = pcms::csys::CreateCylindricalRZ();
  pcms::UniformGrid<2> grid;
  grid.bot_left = {1.0, -0.5};
  grid.edge_length = {2.0, 1.0};
  grid.divisions = {6, 4};
  const auto source_space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, system), 1, 1);
  grid.divisions = {12, 8};
  const auto target_geometry = pcms::uniform_grid::CreateGeometry(grid, system);
  auto source = source_space->CreateFunction<Real>(
    "n", pcms::test::ComponentsBasis(*source_space));
  SetDofs(source, Pseudo(static_cast<std::size_t>(
                           source_space->GetLayout()->GetNumDofHolder()),
                         7));
  const Real expected = pcms::L2Integral(source_space).Apply(source);
  for (const int order : {0, 1}) {
    const auto target_space =
      pcms::CreateLagrangeSpace(target_geometry, order, 1);
    const pcms::L2Integral integral(target_space);
    for (const auto mass :
         {pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped}) {
      CAPTURE(order, static_cast<int>(mass));
      pcms::ProjectionOptions options;
      options.mass_matrix_type = mass;
      options.solve.relative_tolerance = 1e-14;
      auto target = target_space->CreateFunction<Real>(
        "n", pcms::test::ComponentsBasis(*target_space));
      pcms::CreateProjection(source_space, target_space, options)
        ->Apply(source, target);
      CHECK(integral.Apply(target) == Catch::Approx(expected).epsilon(1e-12));
    }
  }
}

TEST_CASE("BOUT++-like embedded grid with a given metric: quadrature and "
          "conservative projection",
          "[uniform_grid][embedding][projection]")
{
  const Bout coarse(2);
  const Bout fine(4);
  const auto quadrature = fine.fields->CreateQuadrature();
  CHECK(quadrature.points->GetGeometry() == fine.geometry);
  CHECK(Integrate(quadrature, [](const Real*) { return 1.0; }) ==
        Catch::Approx(Wedge * VolumeFactor()).epsilon(1e-12));
  const auto rule = pcms::MapQuadrature(
    fine.geometry,
    pcms::CreateGaussQuadrature(fine.geometry->GetParameterGeometry(), 2));
  CHECK(
    Integrate(rule,
              [](const Real* x) { return x[0] * x[0] * x[1] * x[2] * x[2]; }) ==
    Catch::Approx(VolumeFactor() / 3 / 2 * Wedge * Wedge * Wedge / 3)
      .epsilon(1e-12));

  auto source = coarse.fields->CreateFunction<Real>("P");
  SetDofs(source, Pseudo(static_cast<std::size_t>(
                           coarse.fields->GetLayout()->GetNumDofHolder()),
                         8));
  const Real expected = pcms::L2Integral(coarse.fields).Apply(source);
  for (const auto mass :
       {pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped}) {
    CAPTURE(static_cast<int>(mass));
    pcms::ProjectionOptions options;
    options.mass_matrix_type = mass;
    options.solve.relative_tolerance = 1e-14;
    options.quadrature = quadrature;
    auto target = fine.fields->CreateFunction<Real>("P");
    pcms::CreateProjection(coarse.fields, fine.fields, options)
      ->Apply(source, target);
    CHECK(pcms::L2Integral(fine.fields).Apply(target) ==
          Catch::Approx(expected).epsilon(1e-12));
  }
}
