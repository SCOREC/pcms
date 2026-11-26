// Point clouds and MLS reconstruction in 4..6 dimensions (phase space).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/field/query_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/parameter.hpp"
#include "field_test_utils.hpp"

#include <algorithm>
#include <atomic>
#include <memory>
#include <cmath>
#include <random>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::HostMemorySpace;
using pcms::LO;
using pcms::Real;

namespace
{

constexpr double Tol = 1e-8;

// phase-space coordinates: labels, no metric
std::shared_ptr<const pcms::CoordinateSystem> Phase(int dim)
{
  return pcms::csys::CreateParameterSystem("phase", dim);
}

std::vector<Real> Random(int n, int dim, unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<Real> dist(0.0, 1.0);
  std::vector<Real> pts(static_cast<size_t>(n) * dim);
  for (auto& p : pts) {
    p = dist(gen);
  }
  return pts;
}

// a linear function of the coordinates of point i
Real Linear(const std::vector<Real>& pts, int dim, int i)
{
  Real v = 0.5;
  for (int d = 0; d < dim; ++d) {
    v += (d + 1.0) * pts[static_cast<size_t>(i) * dim + d];
  }
  return v;
}

void SetLinear(pcms::Function<Real>& f, const std::vector<Real>& sites, int dim)
{
  const int n = static_cast<int>(sites.size()) / dim;
  Kokkos::View<Real*, HostMemorySpace> data("data", n);
  for (int i = 0; i < n; ++i) {
    data(i) = Linear(sites, dim, i);
  }
  f.SetDOFHolderDataHost(pcms::ValueView<const Real, HostMemorySpace>(
    pcms::ValueBasis{},
    pcms::Rank2View<const Real, HostMemorySpace>(data.data(), n, 1)));
}

pcms::CoordinateView<HostMemorySpace> HostCoords(std::vector<Real>& pts,
                                                 int dim)
{
  return pcms::CoordinateView<HostMemorySpace>(
    Phase(dim), pcms::Rank2View<Real, HostMemorySpace>(
                  pts.data(), static_cast<LO>(pts.size()) / dim, dim));
}

// Locates `targets` on `geometry` and checks the linear function f at them.
void CheckLinear(const std::shared_ptr<const pcms::FunctionSpace>& space,
                 const pcms::Function<Real>& f, pcms::QueryPoints& query,
                 const std::vector<Real>& targets, int dim, double tol)
{
  const int n = static_cast<int>(targets.size()) / dim;
  auto coords = pcms::test::CreateDeviceCoordinateView(
    targets, space->GetCoordinateSystem(), dim);
  query.Locate(coords.coordinate_view);
  REQUIRE(query.Bound()->NumPoints() == n);
  REQUIRE(query.Bound()->MissingPoints().extent(0) == 0);
  Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 1);
  query.Evaluate(f, pcms::test::TagLike(f, out));
  auto host = Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out);
  for (int i = 0; i < n; ++i) {
    REQUIRE_THAT(host(i, 0), WithinAbs(Linear(targets, dim, i), tol));
  }
}

void RunPointCloud(int dim, int num_sites, int min_supports, double radius)
{
  CAPTURE(dim);
  auto sites = Random(num_sites, dim, 7 + dim);
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(
      HostCoords(sites, dim),
      pcms::NeighbourhoodOptions{.radius = radius,
                                 .min_supports = min_supports}),
    pcms::MLSOptions{.degree = 1, .tol = 1e-12});
  REQUIRE(space->GetCoordinateSystem()->Dimension() == dim);
  auto f = space->CreateFunction<Real>();
  SetLinear(f, sites, dim);
  pcms::QueryPoints query(space->GetGeometry());
  // two batches of different points: the second one takes the rebind path
  CheckLinear(space, f, query, Random(50, dim, 100 + dim), dim, Tol);
  CheckLinear(space, f, query, Random(30, dim, 200 + dim), dim, Tol);
}

// A FieldLayout of an external kind: its DOF holders are the sites of
// `geometry`
class ExternalSitesLayout final : public pcms::FieldLayout
{
public:
  ExternalSitesLayout(std::shared_ptr<const pcms::Geometry> geometry,
                      const std::vector<Real>& sites, int dim)
    : n_(static_cast<LO>(sites.size()) / dim), dim_(dim)
  {
    SetGeometry(std::move(geometry));
    SetCoordinateSystem(Phase(dim));
    sites_ = pcms::test::CreateDeviceRank2View(sites, dim);
    owned_ = std::make_unique<bool[]>(n_);
    std::fill_n(owned_.get(), n_, true);
    gids_.resize(n_);
    classification_.assign(n_, 0);
    ids_.resize(n_);
    for (LO i = 0; i < n_; ++i) {
      gids_[i] = i;
      ids_[i] = i;
    }
  }
  int GetNumComponents() const override { return 1; }
  LO GetNumDofHolder() const override { return n_; }
  pcms::GO GetNumGlobalDofHolder() const override { return n_; }
  bool IsDistributed() const override { return false; }
  pcms::Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override
  {
    return pcms::Rank1View<const bool, HostMemorySpace>(owned_.get(), n_);
  }
  pcms::GlobalIDView<HostMemorySpace> GetGidsHost() const override
  {
    return pcms::GlobalIDView<HostMemorySpace>(gids_.data(), gids_.size());
  }
  pcms::EntOffsetsArray GetEntOffsets() const override
  {
    pcms::EntOffsetsArray offsets{};
    for (auto& o : offsets) {
      o = n_;
    }
    offsets[0] = 0;
    return offsets;
  }
  pcms::CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates()
    const override
  {
    return pcms::test::MakeCoords(sites_, GetCoordinateSystem());
  }
  int GetDimension() const override { return dim_; }
  // one DOF holder per vertex of the geometry: the supports are vertex ids
  int GetDOFHolderEntityDim() const override { return 0; }
  pcms::Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override
  {
    return pcms::Rank1View<const LO, HostMemorySpace>(classification_.data(),
                                                      classification_.size());
  }
  pcms::Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override
  {
    return pcms::Rank1View<const LO, HostMemorySpace>(ids_.data(), ids_.size());
  }

private:
  LO n_;
  int dim_;
  Kokkos::View<Real**, DeviceMemorySpace> sites_;
  std::unique_ptr<bool[]> owned_;
  std::vector<pcms::GO> gids_;
  std::vector<LO> ids_;
  std::vector<LO> classification_;
};

// An external geometry kind: the supports of every query point are the
// `count` first sites (a stand-in for a kind's own search); it cannot be bound.
class FixedSupportsGeometry final : public pcms::Geometry
{
public:
  FixedSupportsGeometry(int dim, LO count) : system_(Phase(dim)), count_(count)
  {
  }
  const std::shared_ptr<const pcms::CoordinateSystem>& GetCoordinateSystem()
    const noexcept override
  {
    return system_;
  }
  pcms::SupportResults Neighbourhood(
    pcms::CoordinateView<DeviceMemorySpace> coords, int /*site_dim*/,
    const pcms::NeighbourhoodOptions& /*options*/,
    const pcms::OwnerHints& /*owners*/ = {}) const override
  {
    ++calls;
    const LO n = static_cast<LO>(coords.GetValues().extent(0));
    Kokkos::View<LO*, HostMemorySpace> ptr("ptr", n + 1);
    Kokkos::View<LO*, HostMemorySpace> idx("idx", std::size_t(n) * count_);
    Kokkos::View<Real*, HostMemorySpace> radii2("radii2", n);
    for (LO i = 0; i < n; ++i) {
      ptr(i + 1) = (i + 1) * count_;
      radii2(i) = 25.0;
      for (LO k = 0; k < count_; ++k) {
        idx(i * count_ + k) = k;
      }
    }
    Kokkos::View<LO*, DeviceMemorySpace> d_ptr("ptr", n + 1);
    Kokkos::View<LO*, DeviceMemorySpace> d_idx("idx", std::size_t(n) * count_);
    Kokkos::View<Real*, DeviceMemorySpace> d_radii2("radii2", n);
    Kokkos::deep_copy(d_ptr, ptr);
    Kokkos::deep_copy(d_idx, idx);
    Kokkos::deep_copy(d_radii2, radii2);
    return {d_ptr, d_idx, d_radii2};
  }
  mutable std::atomic<int> calls{0};

protected:
  std::shared_ptr<pcms::BoundPoints> BindImpl(
    const pcms::BindRequest&) const override
  {
    throw pcms::pcms_error("FixedSupportsGeometry cannot be bound");
  }

private:
  std::shared_ptr<const pcms::CoordinateSystem> system_;
  LO count_;
};

} // namespace

TEST_CASE("csys::CreateParameterSystem: coordinates of dimension 1..6")
{
  for (int dim = 1; dim <= 6; ++dim) {
    auto system = Phase(dim);
    REQUIRE(system->Dimension() == dim);
    REQUIRE(*system == *system);
    REQUIRE_FALSE(*system == *Phase(dim));
  }
  REQUIRE_FALSE(*Phase(3) == *Phase(4));
  REQUIRE_FALSE(*Phase(3) == *pcms::csys::CreateCartesian(3));
  REQUIRE_THROWS_WITH(Phase(0), ContainsSubstring("1..6"));
  REQUIRE_THROWS_WITH(Phase(7), ContainsSubstring("1..6"));
  REQUIRE_THROWS_AS(Phase(-1), pcms::pcms_error);
  // Cartesian stays 1..3
  REQUIRE_THROWS_AS(pcms::csys::CreateCartesian(4), pcms::pcms_error);
}

TEST_CASE("point-cloud MLS of degree 1 reproduces a linear function in 5D")
{
  RunPointCloud(5, 400, 30, 1.0);
}

TEST_CASE("point-cloud MLS in 4D and 6D")
{
  RunPointCloud(4, 300, 25, 1.0);
  RunPointCloud(6, 500, 40, 1.2);
}

TEST_CASE("a point cloud of more than 6 coordinates is refused")
{
  std::vector<Real> pts(std::size_t(7) * 10, 0.0);
  pcms::Rank2View<Real, HostMemorySpace> view(pts.data(), 10, 7);
  REQUIRE_THROWS_AS(
    pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::detail::DeferredCartesian(), view)),
    pcms::pcms_error);
}

TEST_CASE("a phase-space point cloud is refused by the physical-space maps")
{
  auto sites = Random(60, 3, 3);
  auto geometry = pcms::point_cloud::CreateGeometry(HostCoords(sites, 3));
  auto targets = Random(5, 3, 4);
  auto coords = pcms::test::CreateDeviceCoordinateView(
    targets, geometry->GetCoordinateSystem(), 3);
  REQUIRE_THROWS_AS(geometry->Bind(coords.coordinate_view, {},
                                   pcms::CreateCartesianToCylindrical()),
                    pcms::pcms_error);
}

TEST_CASE("an external geometry kind provides the supports of a 5D cloud")
{
  constexpr int Dim = 5;
  auto sites = Random(200, Dim, 11);
  auto kind = std::make_shared<FixedSupportsGeometry>(Dim, 200);
  auto layout = std::make_shared<ExternalSitesLayout>(kind, sites, Dim);
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(layout),
    pcms::MLSOptions{.degree = 1, .tol = 1e-12});
  auto f = space->CreateFunction<Real>();
  SetLinear(f, sites, Dim);
  pcms::QueryPoints query(space->GetGeometry());
  CheckLinear(space, f, query, Random(20, Dim, 12), Dim, Tol);
  REQUIRE(kind->calls == 1);
  CheckLinear(space, f, query, Random(10, Dim, 13), Dim, Tol);
  REQUIRE(kind->calls == 2);
}

namespace
{

constexpr Real TwoPi = 6.283185307179586;

// the sawtooth of theta: linear (slope 3) on the branch (-pi, pi], so it is
// linear in the coordinate unwrapped around a point near 0
Real Sawtooth(Real theta, Real x)
{
  return 3.0 * (theta > TwoPi / 2 ? theta - TwoPi : theta) + 2.0 * x + 1.0;
}

std::shared_ptr<const pcms::FunctionSpace> MakeSpace(
  std::vector<Real>& sites,
  std::shared_ptr<const pcms::CoordinateSystem> system,
  pcms::NeighbourhoodOptions options)
{
  return pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(
      pcms::CoordinateView<HostMemorySpace>(
        std::move(system), pcms::Rank2View<Real, HostMemorySpace>(
                             sites.data(), sites.size() / 2, 2)),
      options),
    pcms::MLSOptions{.degree = 1, .tol = 1e-12});
}

} // namespace

TEST_CASE("parameter system axes: periods, validation and equality")
{
  auto system =
    pcms::csys::CreateParameterSystem("theta_energy", 2, {TwoPi, std::nullopt});
  REQUIRE(system->Dimension() == 2);
  REQUIRE(system->Kind() == "theta_energy");
  REQUIRE(system->Period(0) == std::optional<Real>(TwoPi));
  REQUIRE_FALSE(system->Period(1).has_value());
  REQUIRE_FALSE(*system == *pcms::csys::CreateParameterSystem(
                             "theta_energy", 2, {TwoPi, std::nullopt}));
  REQUIRE_FALSE(Phase(2)->HasMetric());
  REQUIRE_THROWS_AS(pcms::csys::CreateParameterSystem("a", 1, {-1.0}),
                    pcms::pcms_error);
}

TEST_CASE("point-cloud MLS uses the supports across a period")
{
  // sites only on the far side of the period: theta in [2pi - 0.6, 2pi)
  std::mt19937 gen(5);
  std::uniform_real_distribution<Real> dist(0.0, 1.0);
  std::vector<Real> sites;
  for (int i = 0; i < 200; ++i) {
    sites.push_back(TwoPi - 0.6 * dist(gen));
    sites.push_back(dist(gen));
  }
  // queries near theta = 0 (and one just below the period)
  std::vector<Real> targets{0.05, 0.5, 0.1, 0.3, 0.2, 0.8, TwoPi - 0.02, 0.4};
  const pcms::NeighbourhoodOptions options{
    .radius = 0.5, .min_supports = 8, .adapt_radius = false};
  auto values_at = [&](std::shared_ptr<const pcms::CoordinateSystem> system,
                       int& missing) {
    auto space = MakeSpace(sites, system, options);
    auto f = space->CreateFunction<Real>();
    Kokkos::View<Real*, HostMemorySpace> data("data", 200);
    for (std::size_t i = 0; i < 200; ++i) {
      data(i) = Sawtooth(sites[2 * i], sites[2 * i + 1]);
    }
    f.SetDOFHolderDataHost(pcms::ValueView<const Real, HostMemorySpace>(
      pcms::ValueBasis{},
      pcms::Rank2View<const Real, HostMemorySpace>(data.data(), 200, 1)));
    pcms::QueryPoints query(
      space->GetGeometry(),
      pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -1e30});
    auto coords = pcms::test::CreateDeviceCoordinateView(targets, system, 2);
    query.Locate(coords.coordinate_view);
    missing = static_cast<int>(query.Bound()->MissingPoints().extent(0));
    Kokkos::View<Real**, DeviceMemorySpace> out("out", 4, 1);
    query.Evaluate(f, pcms::test::TagLike(f, out));
    return Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out);
  };

  int missing = -1;
  auto periodic = values_at(
    pcms::csys::CreateParameterSystem("theta_x", 2, {TwoPi, std::nullopt}),
    missing);
  REQUIRE(missing == 0);
  for (std::size_t i = 0; i < 4; ++i) {
    REQUIRE_THAT(periodic(i, 0),
                 WithinAbs(Sawtooth(targets[2 * i], targets[2 * i + 1]), 1e-8));
  }
  // without the period no site is near the queries near 0
  values_at(pcms::csys::CreateParameterSystem("theta_x", 2), missing);
  REQUIRE(missing == 3);
}
