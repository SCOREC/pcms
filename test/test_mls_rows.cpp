#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "pcms/backends/point_cloud/layout.hpp"
#include "pcms/backends/uniform_grid/geometry.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/backends/point_cloud/point_cloud_geometry.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "field_test_utils.hpp"
#include <Kokkos_Timer.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstdio>
#include <random>
#include <vector>
#include "pcms/field/point_rows.hpp"

// Timing and checksum of repeated MLS evaluation, to compare the per-call
// solve with the cached coefficient rows.
TEST_CASE("MLS repeated evaluation timing", "[mls][timing]")
{
  const int n = 40;
  std::vector<pcms::Real> sources;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      sources.push_back(static_cast<pcms::Real>(i) / (n - 1));
      sources.push_back(static_cast<pcms::Real>(j) / (n - 1));
    }
  }
  pcms::Rank2View<pcms::Real, pcms::HostMemorySpace> view(sources.data(), n * n,
                                                          2);
  auto fs = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(
      pcms::CoordinateView<pcms::HostMemorySpace>(
        pcms::detail::DeferredCartesian(), view),
      pcms::NeighbourhoodOptions{.radius = 0.12, .min_supports = 12}),
    pcms::MLSOptions{.degree = 2});
  auto field = fs->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    field, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y) {
      return Kokkos::sin(3 * x) * Kokkos::cos(2 * y) + x * y;
    });
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> uniform(0.05, 0.95);
  std::vector<pcms::Real> points(std::size_t{4000} * 2);
  for (auto& p : points) {
    p = uniform(rng);
  }
  const auto bound = pcms::test::Bind(*fs->GetGeometry(), points);
  Kokkos::Timer setup_timer;
  const auto evaluator = fs->CreatePointEvaluator<pcms::Real>(bound);
  Kokkos::fence();
  const double setup = setup_timer.seconds();
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("out", 4000, 1);
  const auto values = pcms::test::TagLike(field, out);
  evaluator->Evaluate(field, values); // warm up
  Kokkos::fence();
  Kokkos::Timer timer;
  const int repeats = 10;
  for (int r = 0; r < repeats; ++r) {
    evaluator->Evaluate(field, values);
  }
  Kokkos::fence();
  const double seconds = timer.seconds() / repeats;
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, out);
  double sum = 0, weighted = 0;
  for (int i = 0; i < 4000; ++i) {
    sum += host(i, 0);
    weighted += (i + 1) * host(i, 0);
  }
  std::printf("[mls-timing] %.4f s per Evaluate (4000 targets, 1600 sources, "
              "degree 2); evaluator creation %.4f s; checksum %.17g %.17g\n",
              seconds, setup, sum, weighted);
}

namespace
{

using Host2D = Kokkos::View<pcms::Real**, pcms::HostMemorySpace>;

std::shared_ptr<const pcms::FunctionSpace> MakeCloudSpace(
  int n, int num_components, int degree,
  pcms::NeighbourhoodOptions supports = {.radius = 0.45, .min_supports = 8})
{
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> coords(
    "cloud", static_cast<std::size_t>(n) * n, 2);
  auto host = Kokkos::create_mirror_view(coords);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      host(j * n + i, 0) = static_cast<pcms::Real>(i) / (n - 1);
      host(j * n + i, 1) = static_cast<pcms::Real>(j) / (n - 1);
    }
  }
  Kokkos::deep_copy(coords, host);
  auto layout = std::make_shared<pcms::detail::PointCloudLayout>(
    2, coords, pcms::detail::DeferredCartesian(), num_components);
  return pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(layout, supports),
    pcms::MLSOptions{.degree = static_cast<unsigned>(degree)});
}

void SetDofs(pcms::Field<pcms::Real>& field, const Host2D& x)
{
  field.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
      field.GetData().GetValueBasis(), pcms::MakeConstRank2View(x)));
}

Host2D RandomDofs(const pcms::FieldLayout& layout, unsigned seed)
{
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  Host2D x("x", layout.GetNumDofHolder(), layout.GetNumComponents());
  for (std::size_t i = 0; i < x.extent(0); ++i) {
    for (std::size_t c = 0; c < x.extent(1); ++c) {
      x(i, c) = uniform(rng);
    }
  }
  return x;
}

} // namespace

TEST_CASE("MLS evaluator: several components, the same rows for each",
          "[mls][rows]")
{
  const auto space = MakeCloudSpace(9, 3, 1);
  const auto scalar_space = MakeCloudSpace(9, 1, 1);
  auto field = space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    field, [](pcms::Real x, pcms::Real y, pcms::Real* out) {
      out[0] = 1 + x - 2 * y;
      out[1] = 3 * x + y;
      out[2] = -2.0;
    });
  std::mt19937 rng(8);
  std::uniform_real_distribution<double> uniform(0.1, 0.9);
  std::vector<pcms::Real> points(200);
  for (auto& p : points) {
    p = uniform(rng);
  }
  const auto evaluator = space->CreatePointEvaluator<pcms::Real>(
    pcms::test::Bind(*space->GetGeometry(), points));
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("out", 100, 3);
  evaluator->Evaluate(field, pcms::test::TagLike(field, out));
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, out);
  for (std::size_t i = 0; i < 100; ++i) {
    const double x = points[2 * i], y = points[2 * i + 1];
    CHECK(host(i, 0) == Catch::Approx(1 + x - 2 * y).margin(1e-10));
    CHECK(host(i, 1) == Catch::Approx(3 * x + y).margin(1e-10));
    CHECK(host(i, 2) == Catch::Approx(-2.0).margin(1e-10));
  }
  // the scalar space's evaluator gives the first component identically
  auto scalar = scalar_space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    scalar, [](pcms::Real x, pcms::Real y, pcms::Real* out) {
      out[0] = 1 + x - 2 * y;
    });
  const auto scalar_evaluator = scalar_space->CreatePointEvaluator<pcms::Real>(
    pcms::test::Bind(*scalar_space->GetGeometry(), points));
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> scalar_out("out", 100, 1);
  scalar_evaluator->Evaluate(scalar, pcms::test::TagLike(scalar, scalar_out));
  const auto scalar_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, scalar_out);
  for (int i = 0; i < 100; ++i) {
    CHECK(host(i, 0) == Catch::Approx(scalar_host(i, 0)).margin(1e-14));
  }
}

TEST_CASE("MLS evaluator: the transpose is the adjoint and the matrix is the "
          "evaluation",
          "[mls][rows][adjoint]")
{
  const int components = GENERATE(1, 2);
  const auto space = MakeCloudSpace(9, components, 2);
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> uniform(0.1, 0.9);
  std::vector<pcms::Real> points(120);
  for (auto& p : points) {
    p = uniform(rng);
  }
  const int num_points = 60;
  const auto evaluator = space->CreatePointEvaluator<pcms::Real>(
    pcms::test::Bind(*space->GetGeometry(), points));
  REQUIRE(evaluator->SupportsTranspose());
  const auto* source = pcms::detail::AsRowsSource(*evaluator);
  REQUIRE(source != nullptr);
  const auto x = RandomDofs(*space->GetLayout(), 1);
  Host2D y("y", num_points, components);
  {
    const auto r = RandomDofs(*space->GetLayout(), 2);
    for (std::size_t i = 0; i < y.extent(0); ++i) {
      for (std::size_t c = 0; c < y.extent(1); ++c) {
        y(i, c) = r(i % r.extent(0), c);
      }
    }
  }
  auto field = space->CreateFunction<pcms::Real>();
  SetDofs(field, x);
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> e_x("Ex", num_points,
                                                          components);
  evaluator->Evaluate(field, pcms::test::TagLike(field, e_x));
  const auto e_x_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, e_x);
  // the device layout is the one ValueView expects (LayoutLeft on a GPU)
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> y_device("y", num_points,
                                                               components);
  {
    auto y_mirror = Kokkos::create_mirror_view(y_device);
    for (int i = 0; i < num_points; ++i) {
      for (int c = 0; c < components; ++c) {
        y_mirror(i, c) = y(i, c);
      }
    }
    Kokkos::deep_copy(y_device, y_mirror);
  }
  auto functional = space->CreateFunction<pcms::Real>();
  evaluator->EvaluateTranspose(
    pcms::ValueView<const pcms::Real, pcms::DeviceMemorySpace>(
      field.GetData().GetValueBasis(), pcms::MakeConstRank2View(y_device)),
    functional);
  double lhs = 0, scale = 0;
  for (int i = 0; i < num_points; ++i) {
    for (int c = 0; c < components; ++c) {
      lhs += e_x_host(i, c) * y(i, c);
      scale += std::abs(e_x_host(i, c) * y(i, c));
    }
  }
  CHECK(std::abs(pcms::detail::Apply(functional, field) - lhs) <=
        1e-12 * std::max(1.0, scale));

  const auto& rows = source->Rows();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, rows.offsets);
  const auto columns =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, rows.holders);
  const auto values =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, rows.values);
  for (int p = 0; p < num_points; ++p) {
    double row_sum = 0;
    for (int c = 0; c < components; ++c) {
      double sum = 0;
      for (int k = offsets(p); k < offsets(p + 1); ++k) {
        sum += values(k) * x(columns(k), c);
      }
      CHECK(sum == Catch::Approx(e_x_host(p, c)).margin(1e-12));
    }
    for (int k = offsets(p); k < offsets(p + 1); ++k) {
      row_sum += values(k);
    }
    // the polynomial basis has the constants: the weights sum to one
    CHECK(row_sum == Catch::Approx(1.0).margin(1e-10));
  }
}

TEST_CASE("Interpolator::ApplyTranspose with an MLS source",
          "[mls][rows][adjoint]")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {4, 4};
  const auto target_space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, pcms::csys::CreateCartesian(2)), 1,
    1);
  const auto source_space = MakeCloudSpace(9, 1, 1);
  const auto interpolator =
    pcms::CreateInterpolator<pcms::Real>(source_space, target_space);
  REQUIRE(interpolator->SupportsTranspose());
  CHECK_NOTHROW(
    pcms::CreateInterpolator<pcms::Real>(target_space, source_space));
  auto a = source_space->CreateFunction<pcms::Real>();
  auto b = target_space->CreateFunction<pcms::Real>();
  SetDofs(a, RandomDofs(*source_space->GetLayout(), 5));
  SetDofs(b, RandomDofs(*target_space->GetLayout(), 6));
  auto pa = target_space->CreateFunction<pcms::Real>();
  interpolator->Apply(a, pa);
  auto ptb = source_space->CreateFunction<pcms::Real>();
  interpolator->ApplyTranspose(b, ptb);
  CHECK(pcms::detail::DotProd(pa, b) ==
        Catch::Approx(pcms::detail::DotProd(a, ptb)).epsilon(1e-12));

  // loads: the interpolation reproduces constants, so the loads keep their sum
  auto loads = target_space->CreateFunction<pcms::Real>();
  Host2D ones("ones", target_space->GetLayout()->GetNumDofHolder(), 1);
  Kokkos::deep_copy(ones, 1.0);
  SetDofs(loads, ones);
  auto distributed = source_space->CreateFunction<pcms::Real>();
  interpolator->ApplyTranspose(loads, distributed);
  const auto values = distributed.GetDOFHolderDataHost().GetValues();
  double total = 0;
  for (std::size_t i = 0; i < values.extent(0); ++i) {
    total += values(i, 0);
  }
  CHECK(total ==
        Catch::Approx(static_cast<double>(ones.size())).epsilon(1e-10));
}
