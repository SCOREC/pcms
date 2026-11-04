#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/geometry/bound_source.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/backends/omega_h/geometry.hpp>
#include <pcms/backends/point_cloud/geometry.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/omega_h/entity_layout.hpp>
#include <pcms/field/query_points.hpp>
#include "field_test_utils.hpp"
#include <memory>
#include <random>
#include <type_traits>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::BoundSource;
using pcms::DeviceMemorySpace;
using pcms::FunctionSpace;
using pcms::HostMemorySpace;
using pcms::LO;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;

namespace
{

template <typename T>
auto ToHost(const Kokkos::View<T, DeviceMemorySpace>& view)
{
  return Kokkos::create_mirror_view_and_copy(HostMemorySpace(), view);
}

/// Checks that rebinding `previous`, which nothing else holds, to `pts` reuses
/// it and gives what a new bind of them gives: the same found and missing
/// points and the same values of `f`.
std::shared_ptr<const BoundSource> CheckRebind(
  const FunctionSpace& space, const pcms::Function<Real>& f,
  std::shared_ptr<const BoundSource> previous, const std::vector<Real>& pts)
{
  const auto& geometry = *space.GetGeometry();
  auto coords =
    pcms::test::CreateDeviceCoordinateView(pts, geometry.GetCoordinateSystem());
  const auto fresh = geometry.Bind(coords.coordinate_view, previous->Policy());
  const BoundSource* address = previous.get();
  auto rebound = geometry.Rebind(std::move(previous), coords.coordinate_view);
  REQUIRE(rebound.get() == address);
  REQUIRE(rebound->NumPoints() == fresh->NumPoints());
  for (const auto& [a, b] :
       {std::pair{rebound->FoundPoints(), fresh->FoundPoints()},
        std::pair{rebound->MissingPoints(), fresh->MissingPoints()}}) {
    const auto rebound_points = ToHost(a);
    const auto fresh_points = ToHost(b);
    REQUIRE(rebound_points.extent(0) == fresh_points.extent(0));
    for (size_t i = 0; i < rebound_points.extent(0); ++i) {
      REQUIRE(rebound_points(i) == fresh_points(i));
    }
  }
  const auto basis = f.GetData().GetValueBasis();
  const auto rebound_values = pcms::test::EvaluateToHost(
    *space.CreatePointEvaluator<Real>(rebound), f, basis);
  const auto fresh_values = pcms::test::EvaluateToHost(
    *space.CreatePointEvaluator<Real>(fresh), f, basis);
  for (LO i = 0; i < fresh->NumPoints(); ++i) {
    CAPTURE(i);
    REQUIRE_THAT(rebound_values(i, 0),
                 WithinAbs(fresh_values(i, 0), pcms::test::ExactTol));
  }
  return rebound;
}

} // namespace

TEST_CASE("Evaluators created from one bound source share its points")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto field = space->CreateFunction<Real>("f");
  const auto linear = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return pcms::test::linear_f(x, y);
  };
  pcms::test::SetField(field, linear);

  const auto pts = pcms::test::StandardEvalCoords2D();
  const size_t n = pts.size() / 2;
  auto bound = [&] {
    // The bound source keeps what it needs, so the caller's points may go.
    return pcms::test::Bind(*space->GetGeometry(), pts);
  }();
  REQUIRE(bound->NumPoints() == static_cast<pcms::LO>(n));
  REQUIRE(bound->Policy().mode == OutOfBoundsMode::ERROR);

  auto first = space->CreatePointEvaluator<Real>(bound);
  auto second = space->CreatePointEvaluator<Real>(bound);
  pcms::test::CheckEvaluation(*first, field, pts, linear);
  pcms::test::CheckEvaluation(*second, field, pts, linear);
}

TEST_CASE("A bound source is accepted only by sources on its geometry")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto other = pcms::test::MakeP1Space(mesh);
  auto bound =
    pcms::test::Bind(*space->GetGeometry(), pcms::test::StandardEvalCoords2D());

  REQUIRE_THROWS_WITH(other->CreatePointEvaluator<Real>(bound),
                      ContainsSubstring("another geometry"));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(nullptr),
                      ContainsSubstring("must not be null"));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(
                        bound, pcms::DifferentialOperator::Divergence),
                      ContainsSubstring("Divergence is not supported"));
}

TEST_CASE("Binding checks the query coordinate system")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pcms::test::StandardEvalCoords2D(), pcms::csys::CylindricalRZ::Create());
  REQUIRE_THROWS_WITH(space->GetGeometry()->Bind(coords.coordinate_view),
                      ContainsSubstring("is not the geometry's coordinate"));
}

TEST_CASE("Sources without derivatives reject the gradient")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto bound =
    pcms::test::Bind(*space->GetGeometry(), pcms::test::StandardEvalCoords2D());
  REQUIRE_FALSE(space->Supports(pcms::DifferentialOperator::Gradient));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(
                        bound, pcms::DifferentialOperator::Gradient),
                      ContainsSubstring("Gradient is not supported"));
}

TEST_CASE("Binding splits found and missing points and applies the policy")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto field = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    field,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  // the second point is outside the unit square
  const std::vector<Real> pts = {0.5, 0.5, 2.0, 2.0, 0.25, 0.75};

  REQUIRE_THROWS_WITH(pcms::test::Bind(*space->GetGeometry(), pts),
                      ContainsSubstring("outside the source's domain"));
  REQUIRE_THROWS_WITH(
    pcms::test::Bind(*space->GetGeometry(), pts,
                     OutOfBoundsPolicy{OutOfBoundsMode::NEAREST_BOUNDARY}),
    ContainsSubstring("NEAREST_BOUNDARY"));

  auto bound = pcms::test::Bind(*space->GetGeometry(), pts,
                                OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
  const auto found = Kokkos::create_mirror_view_and_copy(HostMemorySpace(),
                                                         bound->FoundPoints());
  const auto missing = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), bound->MissingPoints());
  REQUIRE(found.extent(0) == 2);
  REQUIRE(found(0) == 0);
  REQUIRE(found(1) == 2);
  REQUIRE(missing.extent(0) == 1);
  REQUIRE(missing(0) == 1);

  pcms::test::CheckEvaluationWithFill(
    *space->CreatePointEvaluator<Real>(bound), field, pts, {true, false, true},
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
    -1.0);
}

TEST_CASE("A rebind gives what a new bind of the points gives")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  const auto cartesian = pcms::csys::Cartesian::Create(2);
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {4, 4};
  auto vertices = std::make_shared<pcms::OmegaHEntityLayout>(
    pcms::test::MakeOmegaHGeometry(mesh), pcms::Vertex, 1);
  const std::vector<std::shared_ptr<const FunctionSpace>> spaces{
    pcms::test::MakeP1Space(pcms::CreateOmegaHGeometry(mesh, cartesian)),
    pcms::CreateLagrangeSpace(pcms::CreateUniformGridGeometry(grid, cartesian),
                              1, 1),
    pcms::CreatePolynomialReconstructionSpace(
      pcms::CreatePointCloudGeometry(vertices), pcms::MLSOptions{.degree = 1})};
  // points inside and outside the unit square, the second set of a different
  // count; then fewer and more points than the first bind's
  const std::vector<Real> first{0.25, 0.5, 1.5, 0.5, 0.75, 0.1};
  const std::vector<Real> second{0.6, 0.6, 0.3, 0.2, -0.5, 0.5, 0.9, 0.4};
  const std::vector<Real> fewer{0.4, 0.7};
  const std::vector<Real> more = [] {
    auto pts = pcms::test::StandardEvalCoords2D();
    for (Real x : pcms::test::StandardOutsideCoords2D()) {
      pts.push_back(x);
    }
    return pts;
  }();
  for (const auto& space : spaces) {
    auto f = space->CreateFunction<Real>("f");
    pcms::test::SetField(
      f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
    auto bound =
      pcms::test::Bind(*space->GetGeometry(), first,
                       OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
    for (const auto* pts : {&second, &fewer, &more, &first}) {
      bound = CheckRebind(*space, f, std::move(bound), *pts);
    }
  }
}

TEST_CASE("A rebind reuses the bound source only when nothing else holds it")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  const auto geometry = space->GetGeometry();
  auto f = space->CreateFunction<Real>("f");
  const auto linear = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return pcms::test::linear_f(x, y);
  };
  pcms::test::SetField(f, linear);
  const std::vector<Real> a{0.25, 0.5, 0.75, 0.1};
  const std::vector<Real> b{0.6, 0.6, 0.3, 0.2};
  const std::vector<Real> fewer{0.4, 0.7};
  const std::vector<Real> more = pcms::test::StandardEvalCoords2D();
  const auto coords = [&](const std::vector<Real>& pts) {
    return pcms::test::CreateDeviceCoordinateView(
      pts, geometry->GetCoordinateSystem());
  };
  const auto rebind = [&](std::shared_ptr<const BoundSource> bound,
                          const std::vector<Real>& pts) {
    auto query = coords(pts);
    return geometry->Rebind(std::move(bound), query.coordinate_view);
  };

  auto bound = pcms::test::Bind(*geometry, a);
  const BoundSource* address = bound.get();
  const Real* storage = bound->LocatedPointValues().data();
  // released into the rebind: the same source, written in place
  bound = rebind(std::move(bound), b);
  REQUIRE(bound.get() == address);
  REQUIRE(bound->LocatedPointValues().data() == storage);
  pcms::test::CheckEvaluation(*space->CreatePointEvaluator<Real>(bound), f, b,
                              linear);

  // the caller keeps a view of its points: a new source, the view unchanged
  {
    const auto kept = bound->LocatedPointValues();
    std::weak_ptr<const BoundSource> previous = bound;
    bound = rebind(std::move(bound), a);
    REQUIRE(previous.expired());
    const auto kept_host = ToHost(kept);
    for (size_t i = 0; i < b.size(); ++i) {
      REQUIRE(kept_host(i / 2, i % 2) == b[i]);
    }
  }

  // an evaluator holds it: a new source, the evaluator unchanged
  {
    const auto evaluator = space->CreatePointEvaluator<Real>(bound);
    address = bound.get();
    bound = rebind(std::move(bound), b);
    REQUIRE(bound.get() != address);
    pcms::test::CheckEvaluation(*evaluator, f, a, linear);
  }

  // fewer points reuse the storage (on a layout whose first rows are
  // contiguous); more points reallocate it, in the same source
  address = bound.get();
  storage = bound->LocatedPointValues().data();
  bound = rebind(std::move(bound), fewer);
  REQUIRE(bound.get() == address);
  if constexpr (std::is_same_v<
                  Kokkos::View<Real**, DeviceMemorySpace>::array_layout,
                  Kokkos::LayoutRight>) {
    REQUIRE(bound->LocatedPointValues().data() == storage);
  }
  pcms::test::CheckEvaluation(*space->CreatePointEvaluator<Real>(bound), f,
                              fewer, linear);
  storage = bound->LocatedPointValues().data();
  bound = rebind(std::move(bound), more);
  REQUIRE(bound.get() == address);
  REQUIRE(bound->LocatedPointValues().data() != storage);
  pcms::test::CheckEvaluation(*space->CreatePointEvaluator<Real>(bound), f,
                              more, linear);

  // a bound source of another geometry, or none
  const auto other = pcms::test::MakeP1Space(mesh)->GetGeometry();
  auto query = coords(a);
  REQUIRE_THROWS_WITH(other->Rebind(bound, query.coordinate_view),
                      ContainsSubstring("belongs to another geometry"));
  REQUIRE_THROWS_WITH(geometry->Rebind(nullptr, query.coordinate_view),
                      ContainsSubstring("must not be null"));
}

TEST_CASE("ForEach writes every entry of small and large batches")
{
  for (const pcms::LO n : {1, 3, 256, 257, 10000}) {
    Kokkos::View<pcms::LO*, pcms::DeviceMemorySpace> out(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "out"), n);
    pcms::detail::ForEach(
      "test", n, KOKKOS_LAMBDA(const pcms::LO i) { out(i) = 2 * i + 1; });
    const auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
    for (pcms::LO i = 0; i < n; ++i) {
      REQUIRE(host(i) == 2 * i + 1);
    }
  }
}

namespace
{

/// Rebinds `previous` to `pts` and checks that the binding has the found and
/// missing points and search records of a new bind.
std::shared_ptr<const BoundSource> CheckOmegaHRebind(
  const pcms::Geometry& geometry, std::shared_ptr<const BoundSource> previous,
  const std::vector<Real>& pts)
{
  const int dim = geometry.GetCoordinateSystem()->Dimension();
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry.GetCoordinateSystem(), dim);
  const auto fresh = std::dynamic_pointer_cast<const pcms::OmegaHMeshBinding>(
    geometry.Bind(coords.coordinate_view, previous->Policy()));
  const BoundSource* address = previous.get();
  const auto rebound = std::dynamic_pointer_cast<const pcms::OmegaHMeshBinding>(
    geometry.Rebind(std::move(previous), coords.coordinate_view));
  REQUIRE(rebound.get() == address);
  REQUIRE(rebound->NumPoints() == fresh->NumPoints());
  const auto found = ToHost(rebound->FoundPoints());
  const auto fresh_found = ToHost(fresh->FoundPoints());
  const auto missing = ToHost(rebound->MissingPoints());
  REQUIRE(found.extent(0) == fresh_found.extent(0));
  REQUIRE(missing.extent(0) == ToHost(fresh->MissingPoints()).extent(0));
  // the same entities, with the same reduced coordinates up to rounding (the
  // guessed element computes them from its own barycentric coordinates)
  pcms::test::CheckSameRecords(pcms::test::FoundRecords(*rebound),
                               pcms::test::FoundRecords(*fresh), 1e-12);
  for (size_t k = 0; k < found.extent(0); ++k) {
    CAPTURE(k, found(k));
    REQUIRE(found(k) == fresh_found(k));
  }
  return rebound;
}

/// `n` points in [lo, hi]^dim, interleaved.
std::vector<Real> RandomPoints(std::mt19937& rng, int n, int dim, Real lo,
                               Real hi)
{
  std::uniform_real_distribution<Real> uniform(lo, hi);
  std::vector<Real> pts(static_cast<size_t>(n) * dim);
  for (auto& x : pts) {
    x = uniform(rng);
  }
  return pts;
}

} // namespace

TEST_CASE("An Omega_h rebind from the last elements gives a new bind's")
{
  auto lib = Omega_h::Library{};
  const int dim = GENERATE(2, 3);
  auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 16)
                       : pcms::test::BuildUnitCube(lib, 4);
  const auto geometry =
    pcms::CreateOmegaHGeometry(mesh, pcms::csys::Cartesian::Create(dim));
  std::mt19937 rng(7);
  const int n = 200;
  auto pts = RandomPoints(rng, n, dim, 0.0, 1.0);
  std::shared_ptr<const BoundSource> bound = pcms::test::Bind(
    *geometry, pts, OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});

  SECTION("small moves, as of tracers")
  {
    std::uniform_real_distribution<Real> step(-0.02, 0.02);
    for (int s = 0; s < 10; ++s) {
      for (auto& x : pts) {
        x += step(rng);
      }
      bound = CheckOmegaHRebind(*geometry, std::move(bound), pts);
    }
  }

  SECTION("jumps across the mesh, and other point counts")
  {
    for (int count : {n, n / 2, n, 2 * n}) {
      pts = RandomPoints(rng, count, dim, 0.0, 1.0);
      bound = CheckOmegaHRebind(*geometry, std::move(bound), pts);
    }
  }

  SECTION("points onto vertices and edges, leaving and re-entering")
  {
    // vertices and the midpoints of edges of the mesh (h = 1/16 or 1/4),
    // where the search picks among the elements sharing them
    const Real h = dim == 2 ? 1.0 / 16 : 1.0 / 4;
    std::vector<Real> shared;
    for (int i = 0; i < n; ++i) {
      for (int d = 0; d < dim; ++d) {
        const int steps = (i * (d + 3)) % (dim == 2 ? 9 : 5);
        shared.push_back(h * steps / (i % 2 == 0 ? 1.0 : 2.0));
      }
    }
    bound = CheckOmegaHRebind(*geometry, std::move(bound), shared);
    // half the points leave the domain, then come back
    auto outside = pts;
    for (size_t i = 0; i < outside.size(); i += static_cast<size_t>(2 * dim)) {
      outside[i] += 2.0;
    }
    bound = CheckOmegaHRebind(*geometry, std::move(bound), outside);
    REQUIRE(bound->MissingPoints().extent(0) == static_cast<size_t>(n / 2));
    bound = CheckOmegaHRebind(*geometry, std::move(bound), pts);
    REQUIRE(bound->MissingPoints().extent(0) == 0);
  }
}

TEST_CASE("An Omega_h rebind applies the out-of-bounds policy")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto f = space->CreateFunction<Real>("f");
  const auto linear = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return pcms::test::linear_f(x, y);
  };
  pcms::test::SetField(f, linear);
  const std::vector<Real> inside{0.25, 0.5, 0.75, 0.1};
  const std::vector<Real> leaving{0.26, 0.5, 1.75, 0.1};
  const std::vector<Real> back{0.27, 0.51, 0.74, 0.11};
  const auto coords = [&](const std::vector<Real>& pts) {
    return pcms::test::CreateDeviceCoordinateView(
      pts, space->GetGeometry()->GetCoordinateSystem());
  };

  pcms::QueryPoints query(space->GetGeometry(),
                          OutOfBoundsPolicy{OutOfBoundsMode::ERROR});
  query.Locate(coords(inside).coordinate_view);
  REQUIRE_THROWS_WITH(query.Locate(coords(leaving).coordinate_view),
                      ContainsSubstring("outside the source's domain"));
  query.Locate(coords(back).coordinate_view);
  pcms::test::CheckEvaluation(*space->CreatePointEvaluator<Real>(query.Bound()),
                              f, back, linear);

  // under FILL the point that leaves is missing, and found when it is back
  auto bound = pcms::test::Bind(*space->GetGeometry(), inside,
                                OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
  bound = CheckOmegaHRebind(*space->GetGeometry(), std::move(bound), leaving);
  REQUIRE(ToHost(bound->MissingPoints())(0) == 1);
  bound = CheckOmegaHRebind(*space->GetGeometry(), std::move(bound), back);
  REQUIRE(bound->MissingPoints().extent(0) == 0);
}

namespace
{

int num_allocations = 0;

void CountAllocation(const Kokkos::Profiling::SpaceHandle, const char*,
                     const void*, const uint64_t)
{
  ++num_allocations;
}

} // namespace

TEST_CASE("An Omega_h rebind of as many points allocates nothing")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 16);
  const auto geometry =
    pcms::CreateOmegaHGeometry(mesh, pcms::csys::Cartesian::Create(2));
  // 64 points inside triangles, farther than the moves from their edges, so
  // the guesses hold and no point is searched for on the grid
  std::vector<Real> pts;
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      pts.push_back((2 * i + 0.3) / 16);
      pts.push_back((2 * j + 0.6) / 16);
    }
  }
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem());
  auto moved = pts;
  for (auto& x : moved) {
    x += 1e-3;
  }
  auto moved_coords = pcms::test::CreateDeviceCoordinateView(
    moved, geometry->GetCoordinateSystem());
  std::shared_ptr<const BoundSource> bound =
    geometry->Bind(coords.coordinate_view);
  bound = geometry->Rebind(std::move(bound), moved_coords.coordinate_view);
  const auto* binding =
    dynamic_cast<const pcms::OmegaHMeshBinding*>(bound.get());
  REQUIRE(binding != nullptr);
  const auto* records = binding->Results<2>().data();

  num_allocations = 0;
  Kokkos::Tools::Experimental::set_allocate_data_callback(CountAllocation);
  bound = geometry->Rebind(std::move(bound), coords.coordinate_view);
  bound = geometry->Rebind(std::move(bound), moved_coords.coordinate_view);
  const int rebind_allocations = num_allocations;
  // the count sees allocations: a new bind allocates
  const auto fresh = geometry->Bind(coords.coordinate_view);
  Kokkos::Tools::Experimental::set_allocate_data_callback(nullptr);
  REQUIRE(rebind_allocations == 0);
  REQUIRE(num_allocations > 0);
  REQUIRE(bound.get() == binding);
  REQUIRE(binding->Results<2>().data() == records);
}
