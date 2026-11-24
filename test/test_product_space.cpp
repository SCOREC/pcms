#include "periodic_test_utils.hpp"
#include "field_test_utils.hpp"
#include <Omega_h_library.hpp>
#include <pcms/field/periodic_axis.hpp>
#include <pcms/geometry/product_geometry.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include <pcms/configuration.hpp>
#include <pcms/transfer/forms.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <cmath>
#include <numbers>

namespace
{

using pcms::Real;
using pcms::test::EvaluateAt;
using pcms::test::SetDofs;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

std::shared_ptr<const pcms::Geometry> Periodic(int cells, Real period = TwoPi)
{
  return pcms::test::PeriodicGrid(cells, period);
}

std::vector<std::vector<Real>> Points(int n, int dim, Real lo, Real hi)
{
  std::vector<std::vector<Real>> points;
  for (int i = 0; i < n; ++i) {
    std::vector<Real> x;
    for (int d = 0; d < dim; ++d) {
      x.push_back(
        lo +
        (hi - lo) *
          std::fmod(0.1234 + 0.6180339 * (i + 1) * (d + 1) + 0.31 * d, 1.0));
    }
    points.push_back(x);
  }
  return points;
}

// the DOFs (i, k) = a_i b_k of the product of two spaces
std::vector<Real> Outer(const std::vector<Real>& a, const std::vector<Real>& b)
{
  std::vector<Real> c;
  for (const Real x : a) {
    for (const Real y : b) {
      c.push_back(x * y);
    }
  }
  return c;
}

// f(z) = 1 + 0.5 cos z - 0.3 sin 2 z, in the DOFs of modes {0, 1, 2}
const std::vector<Real> ZCoefficients = {1.0, 0.5, 0.0, 0.0, -0.3};
Real Z(Real z)
{
  return 1 + 0.5 * std::cos(z) - 0.3 * std::sin(2 * z);
}
Real DZ(Real z)
{
  return -0.5 * std::sin(z) - 0.6 * std::cos(2 * z);
}

// The nodal values of f on the DOF holders of a layout with coordinates
template <typename F>
std::vector<Real> NodalValues(const pcms::FieldLayout& layout, F f)
{
  const auto device = layout.GetDOFHolderCoordinates().GetValues();
  Kokkos::View<Real**, pcms::DeviceMemorySpace> copy("x", device.extent(0),
                                                     device.extent(1));
  pcms::detail::CopyDeviceRank2ViewToDeviceView(copy, device);
  const auto x =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, copy);
  std::vector<Real> values;
  for (std::size_t i = 0; i < x.extent(0); ++i) {
    values.push_back(f(x(i, 0), x.extent(1) > 1 ? x(i, 1) : 0.0));
  }
  return values;
}

// An evaluation-only source: an analytic function of the 3D points
struct Source
{
  KOKKOS_INLINE_FUNCTION Real operator()(Real x, Real y, Real z) const
  {
    return (x + 2 * y) * (1 + 0.5 * Kokkos::cos(z) - 0.3 * Kokkos::sin(2 * z));
  }
};

class SourceEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  explicit SourceEvaluator(std::shared_ptr<const pcms::BoundPoints> bound)
    : PointEvaluator(std::move(bound))
  {
  }
  void EvaluateImpl(
    const pcms::Field<Real>&,
    pcms::ValueView<Real, pcms::DeviceMemorySpace> out) const override
  {
    const auto values = out.GetValues();
    const auto points = Bound().LocatedPoints().GetValues();
    pcms::detail::ForEach(
      "source", static_cast<pcms::LO>(points.extent(0)),
      KOKKOS_LAMBDA(const pcms::LO p) {
        values(p, 0) = Source{}(points(p, 0), points(p, 1), points(p, 2));
      });
  }
};

class SourceFactory final : public pcms::PointEvaluatorFactory
{
public:
  explicit SourceFactory(std::shared_ptr<const pcms::Geometry> geometry)
    : geometry_(std::move(geometry))
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
    pcms::DifferentialOperator) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<SourceEvaluator>(std::move(bound)));
  }

private:
  std::shared_ptr<const pcms::Geometry> geometry_;
};

} // namespace

TEST_CASE("Product space: the layout's global ids and sizes",
          "[product][space]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 2);
  const auto plane = pcms::CreateLagrangeSpace(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
    1, {.backend = pcms::LagrangeBackend::OmegaH});
  const auto fourier = pcms::CreateFourierSpace(Periodic(4), {0, 1});
  const auto space = pcms::CreateProductSpace(plane, fourier);
  const auto& layout = *space->GetLayout();
  const auto na = plane->GetLayout()->GetNumDofHolder();
  REQUIRE(layout.GetNumDofHolder() == na * 3);
  CHECK(layout.GetNumGlobalDofHolder() == static_cast<pcms::GO>(na) * 3);
  CHECK(layout.GetNumComponents() == 1);
  CHECK(layout.GetDimension() == 3);
  const auto gids = layout.GetGidsHost();
  const auto gid_a = plane->GetLayout()->GetGidsHost();
  for (pcms::LO i = 0; i < na; ++i) {
    for (int k = 0; k < 3; ++k) {
      CHECK(gids(i * 3 + k) == gid_a(i) * 3 + k);
    }
  }
  CHECK_THROWS(pcms::CreateProductSpace(plane, nullptr));
  CHECK_THROWS(pcms::CreateProductSpace(
    plane, fourier,
    pcms::CreateProductGeometry(plane->GetGeometry(), Periodic(4))));
}

TEST_CASE("Product space: value and gradient of representable functions",
          "[product][space]")
{
  // Fourier x Fourier: u = cos(x) sin(2 z)... as DOFs a_i b_k
  {
    const Real period = 2.0;
    const auto xs = pcms::CreateFourierSpace(Periodic(4, period), {0, 1}, 1);
    const auto zs = pcms::CreateFourierSpace(Periodic(5), {0, 1, 2});
    const auto space = pcms::CreateProductSpace(xs, zs);
    const std::vector<Real> a = {0.5, 1.5, -0.5};
    const auto dofs = Outer(a, ZCoefficients);
    auto u = space->CreateFunction<Real>();
    SetDofs(u, dofs);
    const Real w = TwoPi / period;
    const auto fx = [&](Real x) {
      return a[0] + a[1] * std::cos(w * x) + a[2] * std::sin(w * x);
    };
    const auto dfx = [&](Real x) {
      return w * (-a[1] * std::sin(w * x) + a[2] * std::cos(w * x));
    };
    const auto points = Points(17, 2, -3.0, 7.0);
    const auto value = EvaluateAt(*space, u, points);
    const auto gradient =
      EvaluateAt(*space, u, points, pcms::DifferentialOperator::Gradient);
    for (std::size_t p = 0; p < points.size(); ++p) {
      const Real x = points[p][0];
      const Real z = points[p][1];
      CHECK(value[p][0] == Catch::Approx(fx(x) * Z(z)).margin(1e-11));
      REQUIRE(gradient[p].size() == 2);
      CHECK(gradient[p][0] == Catch::Approx(dfx(x) * Z(z)).margin(1e-10));
      CHECK(gradient[p][1] == Catch::Approx(fx(x) * DZ(z)).margin(1e-10));
    }
  }
  // piecewise linear x Fourier, and nested: (Lagrange x Fourier) x Fourier
  {
    const auto lagrange = pcms::CreateLagrangeSpace(Periodic(6), 1, 1);
    const auto zs = pcms::CreateFourierSpace(Periodic(5), {0, 1, 2});
    const auto space = pcms::CreateProductSpace(lagrange, zs);
    const std::vector<Real> a = {1, 4, -2, 0.5, 3, 7};
    auto u = space->CreateFunction<Real>();
    SetDofs(u, Outer(a, ZCoefficients));
    const Real h = TwoPi / 6;
    // at a node the Lagrange factor is its value; its slope is constant on a
    // cell
    const auto nodes = EvaluateAt(*space, u, {{2 * h, 1.1}, {5 * h, 0.2}});
    CHECK(nodes[0][0] == Catch::Approx(a[2] * Z(1.1)).margin(1e-12));
    CHECK(nodes[1][0] == Catch::Approx(a[5] * Z(0.2)).margin(1e-12));
    const auto slope = EvaluateAt(*space, u, {{2.25 * h, 0.7}},
                                  pcms::DifferentialOperator::Gradient);
    CHECK(slope[0][0] ==
          Catch::Approx((a[3] - a[2]) / h * Z(0.7)).margin(1e-11));
    CHECK(slope[0][1] ==
          Catch::Approx((0.75 * a[2] + 0.25 * a[3]) * DZ(0.7)).margin(1e-11));
    // nested: the third factor multiplies by 1 + 2 cos(y)
    const auto ys = pcms::CreateFourierSpace(Periodic(3), {0, 1});
    const auto nested = pcms::CreateProductSpace(space, ys);
    auto v = nested->CreateFunction<Real>();
    SetDofs(v, Outer(Outer(a, ZCoefficients), {1.0, 2.0, 0.0}));
    const auto n = EvaluateAt(*nested, v, {{2 * h, 1.1, 0.9}});
    CHECK(n[0][0] ==
          Catch::Approx(a[2] * Z(1.1) * (1 + 2 * std::cos(0.9))).margin(1e-11));
  }
  // Omega_h P1 x Fourier: u = (x + 2 y) f(z)
  {
    Omega_h::Library lib;
    auto mesh = pcms::test::BuildSquareMesh(lib, 3);
    const auto plane = pcms::CreateLagrangeSpace(
      pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
      1, {.backend = pcms::LagrangeBackend::OmegaH});
    const auto zs = pcms::CreateFourierSpace(Periodic(5), {0, 1, 2});
    const auto space = pcms::CreateProductSpace(plane, zs);
    auto u = space->CreateFunction<Real>();
    SetDofs(u, Outer(NodalValues(*plane->GetLayout(),
                                 [](Real x, Real y) { return x + 2 * y; }),
                     ZCoefficients));
    const auto points = Points(13, 3, 0.05, 0.95);
    auto pts = points;
    for (auto& p : pts) {
      p[2] = (p[2] - 0.5) * 20;
    }
    const auto value = EvaluateAt(*space, u, pts);
    for (std::size_t p = 0; p < pts.size(); ++p) {
      CHECK(value[p][0] ==
            Catch::Approx((pts[p][0] + 2 * pts[p][1]) * Z(pts[p][2]))
              .margin(1e-11));
    }
    // the plane space has no gradient, so neither has the product
    CHECK_FALSE(space->Supports(pcms::DifferentialOperator::Gradient));
  }
}

TEST_CASE("Product space: <E x, y> = <x, E^T y> and the matrix is E",
          "[product][adjoint]")
{
  const auto periodic = Periodic(5, 3.0);
  const auto fourier = pcms::CreateFourierSpace(Periodic(4), {0, 1, 3});
  pcms::test::CheckAdjointAndMatrix(
    pcms::CreateProductSpace(pcms::CreateLagrangeSpace(periodic, 1, 1),
                             fourier),
    Points(23, 2, -4.0, 8.0));
  pcms::test::CheckAdjointAndMatrix(
    pcms::CreateProductSpace(pcms::CreateSplineSpace(periodic),
                             pcms::CreatePeriodicHermiteSpace(Periodic(4))),
    Points(23, 2, -4.0, 8.0));
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto plane = pcms::CreateLagrangeSpace(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
    1, {.backend = pcms::LagrangeBackend::OmegaH});
  auto pts = Points(29, 3, 0.02, 0.98);
  for (auto& p : pts) {
    p[2] = (p[2] - 0.5) * 20;
  }
  pcms::test::CheckAdjointAndMatrix(pcms::CreateProductSpace(plane, fourier),
                                    pts);
}

TEST_CASE("Product space: the natural quadrature is the tensor product and "
          "the mass matrix is the Kronecker product of the factors'",
          "[product][quadrature]")
{
  const auto xs = pcms::CreateFourierSpace(Periodic(4, 2.0), {0, 1});
  const auto zs = pcms::CreateFourierSpace(Periodic(5), {0, 1, 2});
  const auto space = pcms::CreateProductSpace(xs, zs);
  const pcms::transfer::detail::BilinearForm mass(space,
                                                  space->CreateQuadrature());
  const std::vector<Real> mass_x = {2.0, 1.0, 1.0};
  const std::vector<Real> mass_z = {TwoPi, TwoPi / 2, TwoPi / 2, TwoPi / 2,
                                    TwoPi / 2};
  const std::size_t n = mass_x.size() * mass_z.size();
  for (std::size_t j = 0; j < n; ++j) {
    std::vector<Real> unit(n, 0.0);
    unit[j] = 1.0;
    auto e = space->CreateFunction<Real>();
    SetDofs(e, unit);
    auto column = space->CreateFunction<Real>();
    mass.Apply(e, column);
    const auto values = pcms::test::GetDofs(column);
    for (std::size_t i = 0; i < n; ++i) {
      const Real expected =
        i != j ? 0.0 : mass_x[i / mass_z.size()] * mass_z[i % mass_z.size()];
      CHECK(values[i] == Catch::Approx(expected).margin(1e-12));
    }
  }
}

#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
TEST_CASE("Product space: projection of an analytic source (matrix-free)",
          "[product][projection]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto plane = pcms::CreateLagrangeSpace(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
    1, {.backend = pcms::LagrangeBackend::OmegaH});
  const auto periodic = Periodic(5);
  const auto geometry =
    pcms::CreateProductGeometry(plane->GetGeometry(), periodic);
  // a rule of the product of the plane's and 7 equispaced points in z, exact
  // for the source times any target below
  const auto rich = pcms::CreateProductSpace(
    plane, pcms::CreateFourierSpace(periodic, {0, 1, 2, 3}), geometry);
  for (const int modes : {1, 3}) {
    CAPTURE(modes);
    std::vector<int> list;
    for (int m = 0; m < modes; ++m) {
      list.push_back(m);
    }
    const auto space = pcms::CreateProductSpace(
      plane, pcms::CreateFourierSpace(periodic, list), geometry);
    pcms::ProjectionOptions options;
    options.quadrature = rich->CreateQuadrature();
    auto target = space->CreateFunction<Real>();
    const auto source = std::make_shared<SourceFactory>(geometry);
    pcms::CreateProjection(source, space, options)
      ->Apply(space->CreateFunction<Real>(), target);
    auto pts = Points(15, 3, 0.05, 0.95);
    for (auto& p : pts) {
      p[2] = (p[2] - 0.5) * 20;
    }
    const auto value = EvaluateAt(*space, target, pts);
    for (std::size_t p = 0; p < pts.size(); ++p) {
      // modes = 3 reproduces the source; modes = 1 its average over z
      const Real f = modes == 3 ? Z(pts[p][2]) : 1.0;
      CHECK(value[p][0] ==
            Catch::Approx((pts[p][0] + 2 * pts[p][1]) * f).margin(1e-9));
    }
  }
}
#endif

TEST_CASE("Product space: twist-shift periodic axis wraps the plane's y into "
          "its interval and shifts z by the number of periods times a function "
          "of x",
          "[product][twist]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto plane = pcms::CreateLagrangeSpace(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
    1, {.backend = pcms::LagrangeBackend::OmegaH});
  const auto zs = pcms::CreateFourierSpace(Periodic(5), {0, 1, 2});
  // the shift of z across y = 1: 0.3 + 0.8 cos(2 pi x), a function of x
  const auto xs = pcms::CreateFourierSpace(Periodic(4, 1.0), {0, 1});
  auto shift = xs->CreateFunction<Real>();
  SetDofs(shift, {0.3, 0.8, 0.0});
  const auto delta = [](Real x) { return 0.3 + 0.8 * std::cos(TwoPi * x); };
  const auto geometry = pcms::CreateProductGeometry(
    plane->GetGeometry(), zs->GetGeometry(),
    pcms::csys::CreateParameterSystem(
      "twisted", 3,
      {{1, 0.0, 1.0,
        pcms::CreateTwistShiftTransform(
          2, std::make_shared<const pcms::Function<Real>>(std::move(shift)),
          {0})},
       {2, 0.0, TwoPi}}));
  const auto space = pcms::CreateProductSpace(plane, zs, geometry);
  auto u = space->CreateFunction<Real>();
  SetDofs(u, Outer(NodalValues(*plane->GetLayout(),
                               [](Real x, Real y) { return x + 2 * y; }),
                   ZCoefficients));
  // (x, y + n, z) is (x, y, z + n delta(x)): u there is (x + 2 y) f(z + n
  // delta)
  for (const int n : {-2, -1, 0, 1, 3}) {
    CAPTURE(n);
    std::vector<std::vector<Real>> points;
    for (const Real x : {0.1, 0.45, 0.8}) {
      points.push_back({x, 0.37 + n, 1.9});
    }
    const auto value = EvaluateAt(*space, u, points);
    for (std::size_t p = 0; p < points.size(); ++p) {
      const Real x = points[p][0];
      CHECK(
        value[p][0] ==
        Catch::Approx((x + 2 * 0.37) * Z(1.9 + n * delta(x))).margin(1e-10));
    }
  }
  // without the periodic axis a point outside the plane is an error
  const auto plain = pcms::CreateProductSpace(plane, zs);
  CHECK_THROWS(
    EvaluateAt(*plain, plain->CreateFunction<Real>(), {{0.5, 1.37, 0.0}}));
  // the adjoint identity holds for points across the wrap
  auto pts = Points(31, 3, 0.02, 0.98);
  for (auto& p : pts) {
    p[1] = (p[1] - 0.5) * 6;
    p[2] = (p[2] - 0.5) * 20;
  }
  pcms::test::CheckAdjointAndMatrix(space, pts);
  // invalid periodic axes
  CHECK_THROWS(pcms::CreateTwistShiftTransform(2, nullptr, {0}));
  // the twist-shift of a second factor would shift the wrong product axes
  CHECK_THROWS(pcms::CreateProductGeometry(zs->GetGeometry(), geometry));
}
