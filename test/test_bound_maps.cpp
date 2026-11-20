#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/backends/omega_h/omega_h_geometry.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/transfer/interpolator.hpp>
#include "field_test_utils.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace pcms
{
using omega_h::CreateGeometry;
} // namespace pcms

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::ComponentScaling;
using pcms::CoordinateMap;
using pcms::DeviceMemorySpace;
using pcms::FunctionSpace;
using pcms::GO;
using pcms::HostMemorySpace;
using pcms::LagrangeBackend;
using pcms::LagrangeOptions;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;
using pcms::ValueBasis;
using pcms::ValueView;
using pcms::detail::SameValueBasis;
namespace values = pcms::values;

namespace
{

std::shared_ptr<const FunctionSpace> BuildCylindricalSource(
  Omega_h::Library& lib, Omega_h::Mesh& mesh, int num_components)
{
  mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1, 8, 8,
                            6, false);
  return CreateLagrangeSpace(
    pcms::CreateGeometry(mesh, pcms::csys::CreateCylindricalRThetaZ()), 1,
    num_components, LagrangeOptions{.backend = LagrangeBackend::OmegaH});
}

std::shared_ptr<const FunctionSpace> BuildCartesianTarget(Omega_h::Library& lib,
                                                          Omega_h::Mesh& mesh,
                                                          int num_components,
                                                          double x_extent = 1.0)
{
  mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, x_extent, 1.0, 1.0, 4,
                            4, 4, false);
  return CreateLagrangeSpace(
    pcms::CreateGeometry(mesh, pcms::csys::CreateCartesian(3)), 1,
    num_components, LagrangeOptions{.backend = LagrangeBackend::OmegaH});
}

ValueBasis CylindricalVector()
{
  return ValueBasis{pcms::csys::CreateCylindricalRThetaZ(),
                    ComponentScaling::Physical, values::Vector};
}

ValueBasis CartesianVector()
{
  return ValueBasis{pcms::csys::CreateCartesian(3), ComponentScaling::Physical,
                    values::Vector};
}

std::shared_ptr<const CoordinateMap> CartToCyl()
{
  return pcms::CreateCartesianToCylindrical();
}

/// Binds the DOF holders of `target` on the geometry of `source` through
/// `map`, as an Interpolator across coordinate systems needs.
std::shared_ptr<const pcms::BoundPoints> BindTargetSites(
  const pcms::FunctionSpace& source, const pcms::FunctionSpace& target,
  std::shared_ptr<const CoordinateMap> map, OutOfBoundsPolicy policy = {})
{
  const auto layout = target.GetLayout();
  return source.GetGeometry()->Bind(layout->GetDOFHolderCoordinates(), policy,
                                    std::move(map));
}

} // namespace

TEST_CASE("Bind with a map: scalar evaluation maps the query points")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 1);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 1);
  auto f = src->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real r, Real theta, Real z) {
      return 2.0 * r + 3.0 * theta - z;
    });

  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  const auto bound =
    src->GetGeometry()->Bind(pts, OutOfBoundsPolicy{}, CartToCyl());
  // the BoundPoints keeps the mapped points, which the evaluators read
  REQUIRE(bound->GetBoundMap() != nullptr);
  REQUIRE(pcms::detail::SameCoordinateSystem(bound->GetCoordinateSystem(),
                                             tgt->GetCoordinateSystem()));
  REQUIRE(static_cast<int>(bound->LocatedPoints().GetValues().extent(0)) == n);
  auto evaluator = src->CreatePointEvaluator<Real>(bound);
  REQUIRE(pcms::detail::SameValueBasis(evaluator->OutputBasis(ValueBasis{}),
                                       ValueBasis{}));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto located =
    pcms::test::CopyCoordinatesToHost(bound->LocatedPoints().GetValues());
  const auto result =
    pcms::test::EvaluateToHost(*evaluator, f, ValueBasis{}, 1);
  for (int i = 0; i < n; ++i) {
    const Real x = coords(i, 0), y = coords(i, 1), z = coords(i, 2);
    CAPTURE(i);
    REQUIRE_THAT(located(i, 0),
                 WithinAbs(std::sqrt(x * x + y * y), pcms::test::EvalTol));
    REQUIRE_THAT(located(i, 1),
                 WithinAbs(std::atan2(y, x), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 0), WithinAbs(2.0 * std::sqrt(x * x + y * y) +
                                           3.0 * std::atan2(y, x) - z,
                                         pcms::test::EvalTol));
  }
}

TEST_CASE("Bind with a map: vector evaluation rotates into the target basis")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src->CreatePointEvaluator<Real>(
    src->GetGeometry()->Bind(pts, OutOfBoundsPolicy{}, CartToCyl()));
  REQUIRE(pcms::detail::SameValueBasis(
    evaluator->OutputBasis(CylindricalVector()), CartesianVector()));
  REQUIRE(pcms::detail::SameValueBasis(
    evaluator->OutputBasis(CartesianVector()), CartesianVector()));
  const ValueBasis natural_cylindrical{pcms::csys::CreateCylindricalRThetaZ(),
                                       ComponentScaling::Natural,
                                       values::Vector};
  // natural components are re-expressed too (the map's Jacobian)
  REQUIRE(pcms::detail::SameValueBasis(
    evaluator->OutputBasis(natural_cylindrical), CartesianVector()));
  REQUIRE_THROWS_WITH(evaluator->OutputBasis(
                        ValueBasis{pcms::csys::CreateCylindricalRZ(),
                                   ComponentScaling::Physical, values::Vector}),
                      ContainsSubstring("takes components in"));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result =
    pcms::test::EvaluateToHost(*evaluator, b, CartesianVector(), 3);
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(coords(i, 0), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(coords(i, 1), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), pcms::test::EvalTol));
  }

  SECTION("evaluate gates the field's stored basis and the output tag")
  {
    auto natural =
      src->CreateFunction<Real>("c", values::Vector, ComponentScaling::Natural);
    Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 3);
    // natural components pass through the map's Jacobian
    REQUIRE_NOTHROW(evaluator->Evaluate(
      natural, ValueView<Real, DeviceMemorySpace>(CartesianVector(),
                                                  pcms::MakeRank2View(out))));
    REQUIRE_THROWS_WITH(
      evaluator->Evaluate(b, ValueView<Real, DeviceMemorySpace>(
                               CylindricalVector(), pcms::MakeRank2View(out))),
      ContainsSubstring("does not match"));
  }

  SECTION("evaluate rejects an output with the wrong number of points")
  {
    REQUIRE(evaluator->NumPoints() == n);
    Kokkos::View<Real**, DeviceMemorySpace> short_out("short_out", n - 1, 3);
    REQUIRE_THROWS_AS(evaluator->Evaluate(b, ValueView<Real, DeviceMemorySpace>(
                                               CartesianVector(),
                                               pcms::MakeRank2View(short_out))),
                      std::invalid_argument);
  }
}

TEST_CASE("Bind with a map: FILL rows survive the rotation untouched")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3, 3.0);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  const Real fill = -7.0;
  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src->CreatePointEvaluator<Real>(src->GetGeometry()->Bind(
    pts, OutOfBoundsPolicy{OutOfBoundsMode::FILL, fill}, CartToCyl()));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result =
    pcms::test::EvaluateToHost(*evaluator, b, CartesianVector(), 3);
  int num_filled = 0;
  for (int i = 0; i < n; ++i) {
    const Real x = coords(i, 0), y = coords(i, 1);
    CAPTURE(i);
    if (std::sqrt(x * x + y * y) > 2.0) {
      ++num_filled;
      REQUIRE(result(i, 0) == fill);
      REQUIRE(result(i, 1) == fill);
      REQUIRE(result(i, 2) == fill);
    } else {
      REQUIRE_THAT(result(i, 0), WithinAbs(x, pcms::test::EvalTol));
      REQUIRE_THAT(result(i, 1), WithinAbs(y, pcms::test::EvalTol));
      REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), pcms::test::EvalTol));
    }
  }
  REQUIRE(num_filled > 0);
  REQUIRE(static_cast<int>(evaluator->FilledPoints().extent(0)) == num_filled);
}

TEST_CASE("Bind with a map: in-bounds values equal to the fill value are "
          "still rotated")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3, 3.0);
  const Real fill = 1.0;
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [fill](Real, Real, Real, Real* out) {
    out[0] = fill;
    out[1] = fill;
    out[2] = fill;
  });

  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src->CreatePointEvaluator<Real>(src->GetGeometry()->Bind(
    pts, OutOfBoundsPolicy{OutOfBoundsMode::FILL, fill}, CartToCyl()));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result =
    pcms::test::EvaluateToHost(*evaluator, b, CartesianVector(), 3);
  int num_rotated = 0;
  for (int i = 0; i < n; ++i) {
    const Real x = coords(i, 0), y = coords(i, 1);
    const Real r = std::sqrt(x * x + y * y);
    CAPTURE(i);
    if (r > 2.0) {
      REQUIRE(result(i, 0) == fill);
      REQUIRE(result(i, 1) == fill);
    } else if (r > 0.0) {
      const Real c = x / r, s = y / r;
      REQUIRE_THAT(result(i, 0),
                   WithinAbs(fill * (c - s), pcms::test::EvalTol));
      REQUIRE_THAT(result(i, 1),
                   WithinAbs(fill * (s + c), pcms::test::EvalTol));
      ++num_rotated;
    }
    REQUIRE_THAT(result(i, 2), WithinAbs(fill, pcms::test::EvalTol));
  }
  REQUIRE(num_rotated > 0);
}

TEST_CASE("Bind with a map: the map must take the query system to the "
          "geometry's")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  const std::vector<Real> pts = {0.5, 0.5, 0.5, 0.2, 0.9, 0.1};
  auto cart_query = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::CreateCartesian(3), 3);
  auto cyl_query = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::CreateCylindricalRThetaZ(), 3);
  const auto& geometry = *src->GetGeometry();

  // a map ending in another system than the geometry's
  REQUIRE_THROWS_WITH(
    geometry.Bind(cart_query.coordinate_view, OutOfBoundsPolicy{},
                  pcms::CreateCylindricalToCartesian()),
    ContainsSubstring("it must take the query coordinate system to the "
                      "geometry's 'cylindrical_rthetaz'"));
  // query points in another system than the map's source
  REQUIRE_THROWS_WITH(
    geometry.Bind(cyl_query.coordinate_view, OutOfBoundsPolicy{}, CartToCyl()),
    ContainsSubstring("the query points are in 'cylindrical_rthetaz'"));
  // no map: the query points must be in the geometry's system
  REQUIRE_THROWS_WITH(geometry.Bind(cart_query.coordinate_view),
                      ContainsSubstring("as Bind's map argument"));
  // a null map is no map
  REQUIRE(
    geometry.Bind(cyl_query.coordinate_view, OutOfBoundsPolicy{}, nullptr)
      ->GetBoundMap() == nullptr);
}

TEST_CASE("Bind with a map: a rebind binds anew with the map")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 1);
  const auto& geometry = *src->GetGeometry();
  const auto map = CartToCyl();
  auto first = pcms::test::CreateDeviceCoordinateView(
    {0.5, 0.5, 0.5, 0.2, 0.9, 0.1}, pcms::csys::CreateCartesian(3), 3);
  auto second = pcms::test::CreateDeviceCoordinateView(
    {0.9, 0.3, 0.4, 0.3, 0.4, 0.8, 0.6, 0.6, 0.2},
    pcms::csys::CreateCartesian(3), 3);
  auto bound = geometry.Bind(first.coordinate_view, OutOfBoundsPolicy{}, map);
  const std::weak_ptr<const pcms::BoundPoints> previous = bound;
  src.reset();
  REQUIRE(bound->GetGeometry().use_count() == 1);
  bound =
    bound->GetGeometry()->Rebind(std::move(bound), second.coordinate_view);
  REQUIRE(previous.expired());
  REQUIRE(bound->GetGeometry().get() == &geometry);
  REQUIRE(bound->GetCoordinateMap() == map);
  const auto fresh =
    geometry.Bind(second.coordinate_view, OutOfBoundsPolicy{}, map);
  const auto located =
    pcms::test::CopyCoordinatesToHost(bound->LocatedPoints().GetValues());
  const auto expected =
    pcms::test::CopyCoordinatesToHost(fresh->LocatedPoints().GetValues());
  REQUIRE(located.extent(0) == 3);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t d = 0; d < 3; ++d) {
      REQUIRE(located(i, d) == expected(i, d));
    }
  }
}

TEST_CASE("Bind with a map: evaluator creation errors")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  const std::vector<Real> cart_pts = {0.5, 0.5, 0.5, 0.2, 0.9, 0.1};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cart_pts, pcms::csys::CreateCartesian(3), 3);

  // No in-tree backend evaluates integer fields at points, so a mapped bind
  // surfaces the backend's refusal; the T != Real guard in
  // TransformedPointEvaluator::Evaluate is defensive.
  REQUIRE_THROWS_WITH(
    src->CreatePointEvaluator<GO>(src->GetGeometry()->Bind(
      query.coordinate_view, OutOfBoundsPolicy{}, CartToCyl())),
    ContainsSubstring("only supports fields of its own value type"));
}

TEST_CASE("Bind with a map: NEAREST_BOUNDARY follows the geometry")
{
  auto cart_to_cyl = CartToCyl();
  // (2, 0, 0.5) has r = 2, beyond the grid's r range
  const std::vector<Real> cart_pts = {0.5, 0.5, 0.5, 2.0, 0.0, 0.5};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cart_pts, pcms::csys::CreateCartesian(3), 3);
  const OutOfBoundsPolicy nearest{OutOfBoundsMode::NEAREST_BOUNDARY};

  SECTION("a uniform grid extrapolates")
  {
    pcms::UniformGrid<3> grid;
    grid.bot_left = {0.1, -3.2, 0.0};
    grid.edge_length = {1.5, 6.4, 1.0};
    grid.divisions = {6, 16, 4};
    auto src =
      CreateLagrangeSpace(pcms::uniform_grid::CreateGeometry(
                            grid, pcms::csys::CreateCylindricalRThetaZ()),
                          1, 1);
    auto bound =
      src->GetGeometry()->Bind(query.coordinate_view, nearest, cart_to_cyl);
    REQUIRE(bound->MissingPoints().extent(0) == 0);
  }

  SECTION("an Omega_h mesh rejects the mode")
  {
    auto lib = Omega_h::Library{};
    Omega_h::Mesh src_mesh(&lib);
    auto src = BuildCylindricalSource(lib, src_mesh, 1);
    REQUIRE_THROWS_WITH(
      src->GetGeometry()->Bind(query.coordinate_view, nearest, cart_to_cyl),
      ContainsSubstring("NEAREST_BOUNDARY"));
  }
}

TEST_CASE("Bind with a map: one map serves many fields, spaces and targets")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_a_mesh(&lib), tgt_b_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt_a = BuildCartesianTarget(lib, tgt_a_mesh, 3);
  auto tgt_b = BuildCartesianTarget(lib, tgt_b_mesh, 3, 0.8);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  auto e =
    src->CreateFunction<Real>("e", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });
  pcms::test::SetFieldComponents(e, [](Real r, Real, Real z, Real* out) {
    out[0] = 2.0 * r;
    out[1] = 0.0;
    out[2] = -z;
  });

  // One map object for every target; one bind per target serves every field
  // on the source space.
  const auto map = CartToCyl();
  const auto to_a = pcms::CreateInterpolator<Real>(
    src, tgt_a, BindTargetSites(*src, *tgt_a, map));
  const auto to_b = pcms::CreateInterpolator<Real>(
    src, tgt_b, BindTargetSites(*src, *tgt_b, map));

  auto b_a = tgt_a->CreateFunction<Real>("b", values::Vector);
  auto e_a = tgt_a->CreateFunction<Real>("e", values::Vector);
  auto b_b = tgt_b->CreateFunction<Real>("b", values::Vector);
  to_a->Apply(b, b_a);
  to_a->Apply(e, e_a);
  to_b->Apply(b, b_b);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    tgt_a->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto rb = b_a.GetDOFHolderDataHost().GetValues();
  const auto re = e_a.GetDOFHolderDataHost().GetValues();
  for (int i = 0; i < static_cast<int>(coords.extent(0)); ++i) {
    CAPTURE(i);
    REQUIRE_THAT(rb(i, 0), WithinAbs(coords(i, 0), pcms::test::EvalTol));
    REQUIRE_THAT(re(i, 0), WithinAbs(2.0 * coords(i, 0), pcms::test::EvalTol));
    REQUIRE_THAT(re(i, 2), WithinAbs(-coords(i, 2), pcms::test::EvalTol));
  }

  SECTION("a mapped bind serves every space on the source geometry")
  {
    auto scalars =
      CreateLagrangeSpace(src->GetGeometry(), 1, 1,
                          LagrangeOptions{.backend = LagrangeBackend::OmegaH});
    auto f = scalars->CreateFunction<Real>("f");
    pcms::test::SetField(
      f, OMEGA_H_LAMBDA(Real r, Real, Real z) { return r + z; });
    const auto bound = BindTargetSites(*src, *tgt_a, map);
    const auto vectors = pcms::test::EvaluateToHost(
      *src->CreatePointEvaluator<Real>(bound), b, CartesianVector(), 3);
    const auto sums = pcms::test::EvaluateToHost(
      *scalars->CreatePointEvaluator<Real>(bound), f, ValueBasis{}, 1);
    for (int i = 0; i < static_cast<int>(coords.extent(0)); ++i) {
      const Real x = coords(i, 0), y = coords(i, 1), z = coords(i, 2);
      CAPTURE(i);
      REQUIRE_THAT(vectors(i, 1), WithinAbs(y, pcms::test::EvalTol));
      REQUIRE_THAT(sums(i, 0), WithinAbs(std::sqrt(x * x + y * y) + z,
                                         pcms::test::EvalTol));
    }
  }
}

TEST_CASE("Bind with a map: maps compose as separate binds")
{
  // Formerly a nested view (cylindrical -> Cartesian -> cylindrical): the
  // cylindrical query points are mapped to Cartesian ones by one bind, those
  // are bound on the source geometry through the second map, and the first
  // bound map takes the Cartesian values back to cylindrical components.
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  const std::vector<Real> cyl_pts = {0.7, 0.3, 0.4, 1.0, 0.9,
                                     0.9, 1.5, 1.2, 0.2};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cyl_pts, pcms::csys::CreateCylindricalRThetaZ(), 3);
  const int n = 3;
  const auto to_cartesian =
    pcms::CreateCylindricalToCartesian()->Bind(query.coordinate_view);
  auto evaluator = src->CreatePointEvaluator<Real>(src->GetGeometry()->Bind(
    to_cartesian->MappedPoints(), OutOfBoundsPolicy{}, CartToCyl()));
  Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 3);
  evaluator->Evaluate(b, ValueView<Real, DeviceMemorySpace>(
                           CartesianVector(), pcms::MakeRank2View(out)));
  const auto& transform = to_cartesian->ValueTransform();
  const auto cylindrical = transform.OutputBasis(CartesianVector());
  REQUIRE(pcms::detail::SameValueBasis(cylindrical, CylindricalVector()));
  transform.Apply(
    ValueView<Real, DeviceMemorySpace>(CartesianVector(),
                                       pcms::MakeRank2View(out)),
    ValueView<Real, DeviceMemorySpace>(cylindrical, pcms::MakeRank2View(out)));
  auto result = Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, out);
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    const auto k = 3 * static_cast<size_t>(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(cyl_pts[k], pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(0.0, pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(cyl_pts[k + 2], pcms::test::EvalTol));
  }
}

TEST_CASE("Interpolator over a mapped bind transfers across systems")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3);
  auto b_src =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  auto b_tgt = tgt->CreateFunction<Real>("b", values::Vector);
  pcms::test::SetFieldComponents(b_src, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  const auto op = pcms::CreateInterpolator<Real>(
    src, tgt, BindTargetSites(*src, *tgt, CartToCyl()));
  op->Apply(b_src, b_tgt);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    tgt->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto result = b_tgt.GetDOFHolderDataHost().GetValues();
  const int n = static_cast<int>(coords.extent(0));
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(coords(i, 0), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(coords(i, 1), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), pcms::test::EvalTol));
  }

  SECTION("the target must declare the basis the evaluator writes")
  {
    auto wrong = tgt->CreateFunction<Real>(
      "w", values::Vector, pcms::csys::CreateCylindricalRThetaZ(),
      ComponentScaling::Physical);
    REQUIRE_THROWS_WITH(op->Apply(b_src, wrong),
                        ContainsSubstring("target's declared basis"));
  }

  SECTION("the BoundPoints must be of the target's DOF holders")
  {
    // as many DOF holders as the target, at other coordinates
    Omega_h::Mesh other_mesh(&lib);
    auto other = BuildCartesianTarget(lib, other_mesh, 3, 0.9);
    REQUIRE(other->GetLayout()->GetNumDofHolder() ==
            tgt->GetLayout()->GetNumDofHolder());
    REQUIRE_THROWS_WITH(pcms::CreateInterpolator<Real>(
                          src, other, BindTargetSites(*src, *tgt, CartToCyl())),
                        ContainsSubstring("target's DOF-holder sites"));
  }

  SECTION("a plain Interpolator refuses to bridge coordinate systems")
  {
    REQUIRE_THROWS_WITH(pcms::CreateInterpolator<Real>(src, tgt),
                        ContainsSubstring("coordinate system"));
  }
}

TEST_CASE("Interpolator over a mapped bind passes through values already in "
          "the target basis")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3);
  // Cartesian components on a cylindrical space; linear in (r, theta, z) so
  // the interpolation is exact and any stray rotation would show.
  auto b_src = src->CreateFunction<Real>("b", values::Vector,
                                         pcms::csys::CreateCartesian(3));
  auto b_tgt = tgt->CreateFunction<Real>("b", values::Vector);
  pcms::test::SetFieldComponents(b_src,
                                 [](Real r, Real theta, Real z, Real* out) {
                                   out[0] = r;
                                   out[1] = theta;
                                   out[2] = z;
                                 });

  const auto op = pcms::CreateInterpolator<Real>(
    src, tgt, BindTargetSites(*src, *tgt, CartToCyl()));
  op->Apply(b_src, b_tgt);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    tgt->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto result = b_tgt.GetDOFHolderDataHost().GetValues();
  for (int i = 0; i < static_cast<int>(coords.extent(0)); ++i) {
    const Real x = coords(i, 0), y = coords(i, 1);
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0),
                 WithinAbs(std::sqrt(x * x + y * y), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 1),
                 WithinAbs(std::atan2(y, x), pcms::test::EvalTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), pcms::test::EvalTol));
  }
}
