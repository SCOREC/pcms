// A BOUT++-like synthetic grid: a uniform logical grid embedded in (R, Z) (a
// shaped torus cross-section) and in (R, phi, Z) with a toroidal shift, by a
// position function, as a reader of a BOUT++ grid file would build it.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "embedding_test_utils.hpp"
#include <cmath>
#include <numbers>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using pcms::DeviceMemorySpace;
using pcms::Real;

namespace
{

// The cross-section: minor radius rho in [a0, a1] (the first logical
// coordinate), poloidal angle theta in [0, 2 pi] (the second), elongated and
// triangular: R = R0 + rho cos(theta + delta sin theta), Z = kappa rho sin
// theta.
struct ShapedTorus
{
  static constexpr Real R0 = 3.0;
  static constexpr Real A0 = 0.4;
  static constexpr Real A1 = 1.0;
  static constexpr Real Kappa = 1.4;
  static constexpr Real Delta = 0.3;

  static std::array<Real, 2> Position(Real xi1, Real xi2)
  {
    const Real rho = A0 + (A1 - A0) * xi1;
    const Real theta = 2 * std::numbers::pi * xi2;
    return {R0 + rho * std::cos(theta + Delta * std::sin(theta)),
            Kappa * rho * std::sin(theta)};
  }
  // J(i, j) = d x_j / d xi_i
  static std::array<std::array<Real, 2>, 2> Jacobian(Real xi1, Real xi2)
  {
    const Real rho = A0 + (A1 - A0) * xi1;
    const Real theta = 2 * std::numbers::pi * xi2;
    const Real u = theta + Delta * std::sin(theta);
    const Real du = 1 + Delta * std::cos(theta);
    const Real drho = A1 - A0;
    const Real dtheta = 2 * std::numbers::pi;
    return {{{drho * std::cos(u), drho * Kappa * std::sin(theta)},
             {-dtheta * rho * std::sin(u) * du,
              dtheta * Kappa * rho * std::cos(theta)}}};
  }
  // the integral of R |det J| over the unit square: the area of the section
  // times R, the volume of revolution per radian. The integrand is smooth and
  // periodic in theta: the midpoint rule converges to roundoff in theta and
  // quadratically in rho.
  static Real MomentOfR()
  {
    const int n = 1500;
    Real sum = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const Real xi1 = (i + 0.5) / n;
        const Real xi2 = (j + 0.5) / n;
        const auto jac = Jacobian(xi1, xi2);
        sum += Position(xi1, xi2)[0] *
               std::abs(jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]);
      }
    }
    return sum / (static_cast<Real>(n) * n);
  }
};

constexpr int Cells = 64;
// the position is a spline of 64 x 64 cells
constexpr Real SplineTol = 1e-5;

struct Section
{
  std::shared_ptr<const pcms::Geometry> grid =
    pcms::test::MakeUnitSquareGrid(Cells);
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakeSplinePosition(grid, ShapedTorus::Position);
  // the seam theta = 0 is the boundary of the logical grid: a point on it is
  // found if its parameters are within the accuracy of the position of the
  // grid, 1e-4 of a cell
  std::shared_ptr<const pcms::EmbeddedGeometry> embedded =
    pcms::CreateEmbeddedGeometry(grid, position,
                                 {.search = {.tolerance = 1e-4}});
};

// Parameters well inside the grid and near its seam (theta = 0).
const std::vector<Real> Parameters = {0.5,  0.5,   0.2,  0.3,   0.8, 0.7,
                                      0.35, 0.003, 0.65, 0.997, 0.5, 0.25};

std::vector<Real> PhysicalPoints(const std::vector<Real>& parameters)
{
  std::vector<Real> points;
  for (size_t i = 0; i < parameters.size(); i += 2) {
    const auto x = ShapedTorus::Position(parameters[i], parameters[i + 1]);
    points.push_back(x[0]);
    points.push_back(x[1]);
  }
  return points;
}

std::vector<Real> Evaluate(
  const pcms::FunctionSpace& space, const pcms::Function<Real>& function,
  const std::shared_ptr<const pcms::BoundPoints>& bound)
{
  const auto evaluator = space.CreatePointEvaluator<Real>(bound);
  Kokkos::View<Real**, DeviceMemorySpace> out("out", bound->NumPoints(), 1);
  evaluator->Evaluate(function, pcms::test::TagLike(function, out));
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
  return std::vector<Real>(host.data(), host.data() + host.extent(0));
}

} // namespace

TEST_CASE("BOUT++-like 2D grid: physical Bind recovers the parameters, "
          "parameter Bind needs no search",
          "[embedding][acceptance]")
{
  const Section section;
  // the cubic spline space of the logical grid, on the embedded geometry
  const auto space = pcms::CreateEmbeddedSpace(
    section.embedded, pcms::CreateSplineSpace(section.grid));
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, KOKKOS_LAMBDA(Real xi1, Real xi2) { return xi1 + 3 * xi2; });

  const auto physical =
    pcms::test::Bind(*section.embedded, PhysicalPoints(Parameters));
  const auto from_physical = Evaluate(*space, f, physical);
  const auto parametric = section.embedded->Bind(
    pcms::test::CreateDeviceCoordinateView(
      Parameters, section.embedded->GetParameterSystem(), 2)
      .coordinate_view);
  const auto from_parameters = Evaluate(*space, f, parametric);
  for (size_t i = 0; i < Parameters.size() / 2; ++i) {
    CAPTURE(i);
    const Real expected = Parameters[2 * i] + 3 * Parameters[2 * i + 1];
    // located by the search: to the error of the position's spline
    CHECK_THAT(from_physical[i], WithinAbs(expected, 10 * SplineTol));
    // located directly: exact
    CHECK_THAT(from_parameters[i], WithinAbs(expected, 1e-10));
  }
}

TEST_CASE("BOUT++-like 2D grid: an analytic field through the embedded "
          "geometry",
          "[embedding][acceptance]")
{
  const Section section;
  const auto analytic = [](Real r, Real z) {
    return std::sin(1.3 * r) * std::cos(2.1 * z);
  };
  // the field's DOFs are its values at the images of the grid vertices
  auto source = pcms::CreateSplineSpace(section.grid);
  auto field = source->CreateFunction<Real>("u");
  pcms::test::SetFieldComponents(field, [&](Real xi1, Real xi2, Real* out) {
    const auto x = ShapedTorus::Position(xi1, xi2);
    out[0] = analytic(x[0], x[1]);
  });
  const auto space = pcms::CreateEmbeddedSpace(section.embedded, source);
  auto u = space->CreateFunction<Real>("u");
  pcms::test::SetFieldComponents(u, [&](Real xi1, Real xi2, Real* out) {
    const auto x = ShapedTorus::Position(xi1, xi2);
    out[0] = analytic(x[0], x[1]);
  });

  SECTION("at points")
  {
    const auto points = PhysicalPoints(Parameters);
    const auto values =
      Evaluate(*space, u, pcms::test::Bind(*section.embedded, points));
    for (size_t i = 0; i < points.size() / 2; ++i) {
      CAPTURE(i);
      CHECK_THAT(values[i],
                 WithinAbs(analytic(points[2 * i], points[2 * i + 1]), 1e-3));
    }
  }
  SECTION("onto a mesh of the (R, Z) plane inside the section")
  {
    Omega_h::Library lib;
    auto mesh =
      Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0);
    // the box R in [R0 + 0.6, R0 + 0.9], Z in [-0.1, 0.1], inside the section
    const auto coords = mesh.coords();
    Omega_h::Write<Real> placed(coords.size());
    Omega_h::parallel_for(
      mesh.nverts(), OMEGA_H_LAMBDA(Omega_h::LO v) {
        placed[2 * v] = ShapedTorus::R0 + 0.6 + 0.3 * coords[2 * v];
        placed[2 * v + 1] = -0.1 + 0.2 * coords[2 * v + 1];
      });
    mesh.set_coords(Omega_h::Reals(placed));
    const auto target = pcms::test::MakeP1Space(
      pcms::omega_h::CreateGeometry(mesh, pcms::csys::CreateCylindricalRZ()));
    const auto interpolator = pcms::CreateInterpolator<Real>(space, target);
    auto v = target->CreateFunction<Real>("v");
    interpolator->Apply(u, v);
    const auto coordinates = pcms::test::CopyCoordinatesToHost(
      target->GetLayout()->GetDOFHolderCoordinates().GetValues());
    const auto result = v.GetDOFHolderDataHost().GetValues();
    for (size_t i = 0; i < coordinates.extent(0); ++i) {
      CAPTURE(i);
      CHECK_THAT(
        result(i, 0),
        WithinAbs(analytic(coordinates(i, 0), coordinates(i, 1)), 1e-3));
    }
  }
}

TEST_CASE("BOUT++-like 2D grid: quadrature, induced metric and curvilinear "
          "components",
          "[embedding][acceptance]")
{
  const Section section;
  SECTION("quadrature integrates the volume of revolution per radian")
  {
    const auto quadrature = pcms::MapQuadrature(
      section.embedded,
      pcms::CreateGaussQuadrature(section.embedded->GetParameterGeometry(), 4));
    const auto weights = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace(), quadrature.weights);
    Real sum = 0;
    for (size_t i = 0; i < weights.extent(0); ++i) {
      sum += weights(i);
    }
    CHECK_THAT(sum, WithinRel(ShapedTorus::MomentOfR(), 1e-5));
  }

  const auto parameters = pcms::test::CreateDeviceCoordinateView(
    Parameters, section.embedded->GetParameterSystem(), 2);
  const auto n = static_cast<int>(Parameters.size() / 2);
  const auto bound = section.embedded->Bind(parameters.coordinate_view);
  SECTION("the induced metric is J^T J")
  {
    const auto metric = bound->Metric();
    Kokkos::View<Real***, DeviceMemorySpace> g("g", n, 2, 2);
    Kokkos::parallel_for(
      "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
      KOKKOS_LAMBDA(const pcms::LO p) {
        for (int a = 0; a < 2; ++a) {
          for (int b = 0; b < 2; ++b) {
            g(p, a, b) = metric(p, a, b);
          }
        }
      });
    const auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), g);
    for (int p = 0; p < n; ++p) {
      const auto j = ShapedTorus::Jacobian(Parameters[std::size_t{2} * p],
                                           Parameters[std::size_t{2} * p + 1]);
      for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
          CAPTURE(p, a, b);
          const Real expected = j[a][0] * j[b][0] + j[a][1] * j[b][1];
          CHECK_THAT(host(p, a, b), WithinAbs(expected, 5e-3));
        }
      }
    }
  }
  SECTION("curvilinear vectors are not yet supported")
  {
    const auto space = pcms::CreateEmbeddedSpace(
      section.embedded, pcms::CreateLagrangeSpace(section.grid, 1, 2));
    const auto e2 = space->CreateFunction<Real>(
      "e2", {.system = section.embedded->GetParameterSystem(),
             .scaling = pcms::ComponentScaling::Natural,
             .variance = pcms::values::Vector});
    CHECK_THROWS_WITH(space->CreatePointEvaluator<Real>(bound)->OutputBasis(
                        e2.GetData().GetValueBasis()),
                      ContainsSubstring("not yet supported"));
  }
}

namespace
{

// A logical grid (x, y, z) with z over a period, embedded in (R, phi, Z) with
// phi = z + shift(x, y): the toroidal shift of a field-aligned grid. All
// components are multilinear, so the trilinear position reproduces them.
struct Shifted
{
  static constexpr Real R0 = 3.0;
  static void Position(const Real* xi, Real* x)
  {
    x[0] = R0 + 0.5 + 0.4 * xi[0] + 0.1 * xi[1] + 0.05 * xi[0] * xi[1];
    x[1] = xi[2] + 0.5 * xi[0] + 0.3 * xi[1] + 0.4 * xi[0] * xi[1];
    x[2] = 0.1 * xi[0] + 0.3 * xi[1] + 0.05 * xi[0] * xi[1];
  }
  // the inverse g^ab of the induced metric g = J^T G J, G = diag(1, R^2, 1)
  static void InverseMetric(const Real* xi, Real* g_inv)
  {
    Real x[3];
    Position(xi, x);
    // d x_j / d xi_a
    const Real j[3][3] = {{0.4 + 0.05 * xi[1], 0.1 + 0.05 * xi[0], 0.0},
                          {0.5 + 0.4 * xi[1], 0.3 + 0.4 * xi[0], 1.0},
                          {0.1 + 0.05 * xi[1], 0.3 + 0.05 * xi[0], 0.0}};
    const Real scale[3] = {1.0, x[0] * x[0], 1.0};
    Real g[3][3] = {};
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        for (int k = 0; k < 3; ++k) {
          g[a][b] += scale[k] * j[k][a] * j[k][b];
        }
      }
    }
    const Real det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                     g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                     g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        // the cofactor of g[b][a]
        const int r0 = (b + 1) % 3, r1 = (b + 2) % 3;
        const int c0 = (a + 1) % 3, c1 = (a + 2) % 3;
        g_inv[3 * a + b] =
          (g[r0][c0] * g[r1][c1] - g[r0][c1] * g[r1][c0]) / det;
      }
    }
  }
};

Real TotalWeight(const std::shared_ptr<const pcms::EmbeddedGeometry>& geometry)
{
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(),
    pcms::MapQuadrature(geometry, pcms::CreateGaussQuadrature(
                                    geometry->GetParameterGeometry(), 2))
      .weights);
  Real sum = 0;
  for (std::size_t i = 0; i < weights.extent(0); ++i) {
    sum += weights(i);
  }
  return sum;
}

} // namespace

TEST_CASE("BOUT++-like 3D grid: toroidal shift, periodic phi and a wedge "
          "symmetry",
          "[embedding][acceptance]")
{
  // a 1/5 toroidal wedge, as BOUT++ with zperiod = 5: z covers one period
  const Real wedge = 2 * std::numbers::pi / 5;
  pcms::UniformGrid<3> grid;
  grid.bot_left = {0.0, 0.0, 0.0};
  grid.edge_length = {1.0, 1.0, wedge};
  grid.divisions = {8, 8, 8};
  const auto logical = pcms::uniform_grid::CreateGeometry(
    grid, pcms::csys::CreateParameterSystem("bout_index", 3));
  const auto space = pcms::CreateLagrangeSpace(logical, 1, 1);
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(logical, 1, 3), Shifted::Position,
    pcms::csys::CreateCylindricalRThetaZ(wedge));
  const auto embedded = pcms::CreateEmbeddedGeometry(logical, position);
  CHECK(embedded->TopologicalDimension() == 3);
  const auto fields = pcms::CreateEmbeddedSpace(embedded, space);
  auto g = fields->CreateFunction<Real>("g");
  pcms::test::SetField(
    g, KOKKOS_LAMBDA(Real x, Real y, Real z) { return x + 2 * y + 3 * z; });

  // (x, y, z) of points in the wedge, including phi in the overlap of the two
  // windows (the shift makes the image wider than the wedge)
  const std::vector<Real> parameters = {0.3, 0.6, 0.5, 0.7, 0.2, 0.2,
                                        0.9, 0.9, 1.1, 0.5, 0.5, 0.05};
  for (const int periods : {0, 1, -3, 7}) {
    CAPTURE(periods);
    std::vector<Real> points;
    for (size_t i = 0; i < parameters.size(); i += 3) {
      Real x[3];
      Shifted::Position(&parameters[i], x);
      // outside the wedge: translated by the symmetry (and, for the
      // cylindrical system, 2 pi is a multiple of it)
      points.insert(points.end(), {x[0], x[1] + periods * wedge, x[2]});
    }
    const auto bound = pcms::test::Bind(*embedded, points);
    // the located points are reduced: phi is within two wedges of the image
    const auto located =
      pcms::test::CopyCoordinatesToHost(bound->LocatedPoints().GetValues());
    for (size_t i = 0; i < located.extent(0); ++i) {
      CHECK(located(i, 1) >= 0);
      CHECK(located(i, 1) <= 2 * wedge + 1.2);
    }
    const auto evaluator = fields->CreatePointEvaluator<Real>(bound);
    Kokkos::View<Real**, DeviceMemorySpace> out("out", 4, 1);
    evaluator->Evaluate(g, pcms::test::TagLike(g, out));
    const auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
    for (int i = 0; i < 4; ++i) {
      CAPTURE(i);
      // the point is the image of (x, y, z) and of its translates by the
      // symmetry (z + k wedge): the field differs by 3 k wedge
      const Real expected = parameters[std::size_t{3} * i] +
                            2 * parameters[std::size_t{3} * i + 1] +
                            3 * parameters[std::size_t{3} * i + 2];
      const Real k = (host(i, 0) - expected) / (3 * wedge);
      CHECK_THAT(k, WithinAbs(std::round(k), 1e-6));
    }
  }
}

TEST_CASE("BOUT++-like 3D grid: a given metric g^ij", "[embedding][acceptance]")
{
  // as a reader of a grid file builds it: the logical grid, the position
  // (R, z + shift, Z) and the contravariant metric in index components
  const Real wedge = 2 * std::numbers::pi / 5;
  pcms::UniformGrid<3> grid;
  grid.bot_left = {0.0, 0.0, 0.0};
  grid.edge_length = {1.0, 1.0, wedge};
  grid.divisions = {8, 8, 8};
  const auto logical = pcms::csys::CreateParameterSystem("bout_index", 3);
  const auto geometry = pcms::uniform_grid::CreateGeometry(grid, logical);
  const auto position = pcms::test::MakePositionFunction(
    pcms::CreateLagrangeSpace(geometry, 1, 3), Shifted::Position,
    pcms::csys::CreateCylindricalRThetaZ(wedge));
  const auto metric_space = pcms::CreateLagrangeSpace(geometry, 1, 9);
  auto g_inv =
    std::make_shared<pcms::Function<Real>>(metric_space->CreateFunction<Real>(
      "g^ij", {.scaling = pcms::ComponentScaling::Natural,
               .variance = pcms::values::Tensor}));
  pcms::test::SetFieldComponents(*g_inv, [](Real x, Real y, Real z, Real* out) {
    const Real xi[3] = {x, y, z};
    Shifted::InverseMetric(xi, out);
  });
  const auto given =
    pcms::CreateEmbeddedGeometry(geometry, position, {.metric = g_inv});
  const auto induced = pcms::CreateEmbeddedGeometry(geometry, position);
  CHECK(given->GetMetric() == g_inv);
  CHECK(given->GetPosition() == position);
  CHECK(*given->GetParameterSystem() == *logical);
  // at the vertices the given metric is the induced one: so are the volume
  // factors (between them it is interpolated, which is not the induced one)
  const std::vector<Real> vertices = {0.25, 0.5, 0.0, 0.75, 0.125, 0.5};
  const auto at = pcms::test::CreateDeviceCoordinateView(vertices, logical, 3);
  const auto given_bound = given->Bind(at.coordinate_view);
  const auto induced_bound = induced->Bind(at.coordinate_view);
  Kokkos::View<Real**, DeviceMemorySpace> volumes("volumes", 2, 2);
  const auto given_volume = given_bound->VolumeFactor();
  const auto induced_volume = induced_bound->VolumeFactor();
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 2),
    KOKKOS_LAMBDA(const pcms::LO p) {
      volumes(p, 0) = given_volume(p);
      volumes(p, 1) = induced_volume(p);
    });
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), volumes);
  for (int p = 0; p < 2; ++p) {
    CAPTURE(p);
    CHECK_THAT(host(p, 0), WithinRel(host(p, 1), 1e-10));
  }
  // and the quadrature uses the given metric: the volume of the wedge
  // differs from the induced one by the interpolation of g^ij
  CHECK(std::abs(TotalWeight(given) - TotalWeight(induced)) >
        1e-3 * TotalWeight(induced));
}
