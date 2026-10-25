#include "mdspan/mdspan.hpp"
#include "pcms.h"
#include <Kokkos_Core.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include "pcms/field/evaluator/spline_interpolator.hpp"
#include "field_test_utils.h"

using namespace pcms;
using pcms::test::CreateDeviceRank1View;

using TestMemorySpace = pcms::DeviceMemorySpace;
using TestExecutionSpace = TestMemorySpace::execution_space;

bool are_equal(double a, double b, double tolerance = 1e-7)
{
  return std::abs(a - b) < tolerance;
}

void createFlatGrid(Rank1View<double, TestMemorySpace> xvec,
                    Rank1View<double, TestMemorySpace> yvec,
                    Kokkos::View<double*, TestMemorySpace> X_flat,
                    Kokkos::View<double*, TestMemorySpace> Y_flat)
{
  size_t nx = xvec.size();
  size_t ny = yvec.size();

  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ny, nx}),
    KOKKOS_LAMBDA(const int i, const int j) {
      size_t idx = i * nx + j;
      X_flat[idx] = xvec[j];
      Y_flat[idx] = yvec[i];
    });
}

void tset(int nth, Kokkos::View<double*, TestMemorySpace> th,
          Kokkos::View<double*, TestMemorySpace> sth,
          Kokkos::View<double*, TestMemorySpace> cth, double thmin,
          double thmax)
{
  Kokkos::parallel_for(
    Kokkos::RangePolicy<TestExecutionSpace>(0, nth),
    KOKKOS_LAMBDA(const int ith) {
      th[ith] = thmin + static_cast<double>(ith) * (thmax - thmin) /
                          static_cast<double>(nth - 1);
      sth[ith] = 2.0 + std::sin(th[ith]);
      cth[ith] = std::cos(th[ith]);
    });
}

void xset(int nx, Kokkos::View<double*, TestMemorySpace> x,
          Kokkos::View<double*, TestMemorySpace> ex, double xmin, double xmax)
{
  if (nx < 2)
    return; // avoid division by zero

  Kokkos::parallel_for(
    Kokkos::RangePolicy<TestExecutionSpace>(0, nx),
    KOKKOS_LAMBDA(const int ix) {
      x[ix] = xmin + static_cast<double>(ix) * (xmax - xmin) /
                       static_cast<double>(nx - 1);
      ex[ix] = std::exp(2.0 * x[ix] - 1.0);
    });
}

void ffset(int num, Kokkos::View<double*, TestMemorySpace> xf,
           Kokkos::View<double*, TestMemorySpace> tf,
           Kokkos::View<double**, TestMemorySpace> f)
{
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num, num}),
    KOKKOS_LAMBDA(const int i, const int j) { f(i, j) = xf[i] * tf[j]; });
}

void bset(Rank1View<double, TestMemorySpace> fx, int nx,
          Rank1View<double, TestMemorySpace> fth, int nth,
          Rank1View<double, TestMemorySpace> bcx1,
          Rank1View<double, TestMemorySpace> bcx2,
          Rank1View<double, TestMemorySpace> bcth1,
          Rank1View<double, TestMemorySpace> bcth2)
{
  Kokkos::parallel_for(
    Kokkos::RangePolicy<TestExecutionSpace>(0, nth),
    KOKKOS_LAMBDA(const int ith) {
      bcx1[ith] = 2.0 * fx(0) * fth(ith);      // df/dx at x(1)
      bcx2[ith] = 2.0 * fx(nx - 1) * fth(ith); // df/dx at x(nx)
    });

  Kokkos::parallel_for(
    Kokkos::RangePolicy<TestExecutionSpace>(0, nx),
    KOKKOS_LAMBDA(const int ix) {
      bcth1(ix) = fx(ix); // df/dth at th = 0 (th[0])
      bcth2(ix) = fx(ix); // df/dth at th = 2π (th[nth - 1])
    });
}

void generate_gt_2d(Rank1View<double, TestMemorySpace> gt_splinv,
                    Rank1View<double, TestMemorySpace> fxtest,
                    Rank1View<double, TestMemorySpace> fthtest, int ntest)
{
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ntest, ntest}),
    KOKKOS_LAMBDA(const int j, const int i) {
      gt_splinv(j * ntest + i) = fxtest(i) * fthtest(j);
    });
}

void compare(Rank2View<double, TestMemorySpace> splinv,
             Rank1View<double, TestMemorySpace> gt_splinv,
             Kokkos::View<double*, HostMemorySpace> res, std::string label)
{
  // Compare the output of splinv with gt_splinv
  double result_fdif = 0.0;
  double result_fdifr = 0.0;
  Kokkos::parallel_reduce(
    "compare_splinv", splinv.extent(0),
    KOKKOS_LAMBDA(const int i, double& max_fdif, double& max_fdifr) {
      if (splinv(i, 0) != 0.0) {
        double dif = Kokkos::abs(splinv(i, 0) - gt_splinv(i));
        max_fdif = Kokkos::max(max_fdif, dif);
        max_fdifr = Kokkos::max(max_fdifr, dif / gt_splinv(i));
      }
    },
    Kokkos::Max<double>(result_fdif), Kokkos::Max<double>(result_fdifr));
  res(0) = result_fdif;
  res(1) = result_fdifr;
  std::cout << "2d " << label << " max absolute difference: " << res(0)
            << ", relative difference: " << res(1) << std::endl;
}

void compare_2d(Rank2View<double, TestMemorySpace> splinv,
                Rank1View<double, TestMemorySpace> gt_splinv,
                Kokkos::View<double*, HostMemorySpace> res, std::string label)
{
  // Compare the output of splinv with gt_splinv
  double result_fdif = 0.0;
  double result_fdifr = 0.0;
  Kokkos::parallel_reduce(
    "compare_splinv", splinv.extent(0),
    KOKKOS_LAMBDA(const int i, double& max_fdif, double& max_fdifr) {
      if (splinv(i, 0) != 0.0) {
        double dif = Kokkos::abs(splinv(i, 0) - gt_splinv(i));
        max_fdif = Kokkos::max(max_fdif, dif);
        max_fdifr =
          Kokkos::max(max_fdifr, dif / (0.5 * (splinv(i, 0) + gt_splinv(i))));
      }
    },
    Kokkos::Max<double>(result_fdif), Kokkos::Max<double>(result_fdifr));
  res(0) = result_fdif;
  res(1) = result_fdifr;
  std::cout << "2d " << label << " max absolute difference: " << res(0)
            << ", relative difference: " << res(1) << std::endl;
}

void reset_2dspan(Rank2View<double, TestMemorySpace> span)
{
  Kokkos::parallel_for(
    "zero_init",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                           {span.extent(0), span.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) { span(i, j) = 0.0; });
}

void dotest1(int ns, Rank1View<double, TestMemorySpace> x,
             Rank1View<double, TestMemorySpace> f, int nt,
             Rank1View<double, TestMemorySpace> xt,
             Rank1View<double, TestMemorySpace> ft,
             Kokkos::View<double*, HostMemorySpace> res_1d)
{
  using std::data;
  using std::size;

  Kokkos::View<double**, TestMemorySpace> splinv_view("splinv_view", nt, 3);
  auto splinv = MakeRank2View(splinv_view);
  reset_2dspan(splinv);

  // interpolator.cspline(x, ns, fspl, 1, 1, 1, 1, wk);
  ExplicitCubicSplineInterpolator<double, TestMemorySpace>
    explicit_interpolator(x, f, BoundaryCondition::FIRST_DERIVATIVE_MATCH, 1,
                          BoundaryCondition::FIRST_DERIVATIVE_MATCH, 1);
  explicit_interpolator.evaluate(xt, splinv);

  compare(splinv, ft, Kokkos::subview(res_1d, Kokkos::make_pair(0, 2)),
          "cspline");
  reset_2dspan(splinv);

  ExplicitCubicSplineInterpolator<double, TestMemorySpace>
    explicit_interpolator_periodic(x, f, BoundaryCondition::PERIODIC, 0,
                                   BoundaryCondition::PERIODIC, 0);
  explicit_interpolator_periodic.evaluate(xt, splinv);

  compare(splinv, ft, Kokkos::subview(res_1d, Kokkos::make_pair(2, 4)),
          "cspline_periodic");
  reset_2dspan(splinv);

  CompactCubicSplineInterpolator<double, TestMemorySpace> compact_interpolator(
    x, f, BoundaryCondition::FIRST_DERIVATIVE_MATCH, 1,
    BoundaryCondition::FIRST_DERIVATIVE_MATCH, 1);
  compact_interpolator.evaluate(xt, splinv);

  compare(splinv, ft, Kokkos::subview(res_1d, Kokkos::make_pair(4, 6)),
          "compact_cubic");
}

void pspltest1(Kokkos::View<double*, HostMemorySpace> res_1d)
{
  const double pi2 = 6.28318530718;
  const double zero = 0.0;
  int inum = 10;

  // Local arrays
  Kokkos::View<double*, TestMemorySpace> zdum_view("zdum_view", 1000);
  Kokkos::View<double*, TestMemorySpace> xtest_view("xtest_view", 1000);
  Kokkos::View<double*, TestMemorySpace> ftest_view("ftest_view", 1000);
  Kokkos::View<double*, TestMemorySpace> x_view("x_view", 10);
  Kokkos::View<double*, TestMemorySpace> zcos_view("zcos_view", 10);
  Kokkos::View<double*, TestMemorySpace> z2sin_view("z2sin_view", 10);

  // Prepare test data
  tset(1000, xtest_view, ftest_view, zdum_view, zero - 0.1, pi2 + 0.1);
  tset(inum, x_view, z2sin_view, zcos_view, zero, pi2);

  auto x = MakeRank1View(x_view);

  auto z2sin = MakeRank1View(z2sin_view);

  auto xtest = MakeRank1View(xtest_view);
  auto ftest = MakeRank1View(ftest_view);

  // Call test function
  dotest1(inum, x, z2sin, 1000, xtest, ftest, res_1d);
}

void pspltest1_long_periodic(Kokkos::View<double*, HostMemorySpace> res_1d)
{
  const double pi2 = 6.28318530718;
  const double zero = 0.0;
  int inum = 16;

  Kokkos::View<double*, TestMemorySpace> zdum_view("zdum_view_long", 1000);
  Kokkos::View<double*, TestMemorySpace> xtest_view("xtest_view_long", 1000);
  Kokkos::View<double*, TestMemorySpace> ftest_view("ftest_view_long", 1000);
  Kokkos::View<double*, TestMemorySpace> x_view("x_view_long", inum);
  Kokkos::View<double*, TestMemorySpace> zcos_view("zcos_view_long", inum);
  Kokkos::View<double*, TestMemorySpace> z2sin_view("z2sin_view_long", inum);

  tset(1000, xtest_view, ftest_view, zdum_view, zero - 0.1, pi2 + 0.1);
  tset(inum, x_view, z2sin_view, zcos_view, zero, pi2);

  auto x = MakeRank1View(x_view);
  auto z2sin = MakeRank1View(z2sin_view);
  auto xtest = MakeRank1View(xtest_view);
  auto ftest = MakeRank1View(ftest_view);

  dotest1(inum, x, z2sin, 1000, xtest, ftest, res_1d);
}

void dotest2(Rank1View<double, TestMemorySpace> x,
             Rank1View<double, TestMemorySpace> fx, int nx,
             Rank1View<double, TestMemorySpace> th,
             Rank1View<double, TestMemorySpace> fth,
             Rank1View<double, TestMemorySpace> dfth, int nth,
             Kokkos::View<double**, TestMemorySpace> f,
             Rank1View<double, TestMemorySpace> bcx1,
             Rank1View<double, TestMemorySpace> bcx2,
             Rank1View<double, TestMemorySpace> bcth1,
             Rank1View<double, TestMemorySpace> bcth2,
             Rank1View<double, TestMemorySpace> xtest,
             Rank1View<double, TestMemorySpace> fxtest,
             Rank1View<double, TestMemorySpace> thtest,
             Rank1View<double, TestMemorySpace> fthtest, int ntest,
             Kokkos::View<double*, HostMemorySpace> res_2d)
{

  Kokkos::View<double*, TestMemorySpace> values_view(
    "values_view", static_cast<std::size_t>(nth) * nx);
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nth, nx}),
    KOKKOS_LAMBDA(const int ith, const int ix) {
      values_view(ix * nth + ith) = f(ix, ith);
    });
  auto values = MakeRank1View(values_view);

  bset(fx, nx, fth, nth, bcx1, bcx2, bcth1, bcth2);
  Kokkos::View<double**, TestMemorySpace> splinv_view("splinv_view", 40000, 10);
  auto splinv = MakeRank2View(splinv_view);
  reset_2dspan(splinv);

  Kokkos::View<double*, TestMemorySpace> xtest_grid_view("xtest_grid_view",
                                                         40000);
  Kokkos::View<double*, TestMemorySpace> thtest_grid_view("thtest_grid_view",
                                                          40000);
  createFlatGrid(xtest, thtest, xtest_grid_view, thtest_grid_view);
  auto xtest_grid = MakeRank1View(xtest_grid_view);
  auto thtest_grid = MakeRank1View(thtest_grid_view);
  Kokkos::View<double*, TestMemorySpace> splinv_gt_view("splinv_view", 40000);
  auto splinv_gt = MakeRank1View(splinv_gt_view);
  generate_gt_2d(splinv_gt, fxtest, fthtest, ntest);

  ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>
    explicit_interpolator(x, th, values,
                          BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcx1,
                          BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcx2,
                          BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcth1,
                          BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcth2);
  explicit_interpolator.evaluate(xtest_grid, thtest_grid, splinv);

  compare_2d(splinv, splinv_gt,
             Kokkos::subview(res_2d, Kokkos::make_pair(0, 2)), "bcspline");
  reset_2dspan(splinv);

  CompactBiCubicSplineInterpolator<double, TestMemorySpace>
    compact_interpolator(x, th, values,
                         BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcx1,
                         BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcx2,
                         BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcth1,
                         BoundaryCondition::FIRST_DERIVATIVE_MATCH, bcth2);
  compact_interpolator.evaluate(xtest_grid, thtest_grid, splinv);

  compare_2d(splinv, splinv_gt,
             Kokkos::subview(res_2d, Kokkos::make_pair(2, 4)), "mkbicub");

  // CompactBiCubicSplineInterpolator<double, TestMemorySpace>
  //     test_interpolator(x, th, values);
}

void pspltest2(Kokkos::View<double*, HostMemorySpace> res_2d)
{
  const double pi2 = 6.28318530718;
  const double zero = 0.0;
  const double one = 1.0;

  Kokkos::View<double*, TestMemorySpace> x1_view("x1_view", 10);
  Kokkos::View<double*, TestMemorySpace> ex1_view("ex1_view", 10);
  Kokkos::View<double*, TestMemorySpace> t1_view("t1_view", 10);
  Kokkos::View<double*, TestMemorySpace> st1_view("st1_view", 10);
  Kokkos::View<double*, TestMemorySpace> ct1_view("ct1_view", 10);
  Kokkos::View<double*, TestMemorySpace> xtest_view("xtest_view", 200);
  Kokkos::View<double*, TestMemorySpace> extest_view("extest_view", 200);
  Kokkos::View<double*, TestMemorySpace> ttest_view("ttest_view", 200);
  Kokkos::View<double*, TestMemorySpace> stest_view("stest_view", 200);
  Kokkos::View<double*, TestMemorySpace> ctest_view("ctest_view", 200);
  Kokkos::View<double**, TestMemorySpace> f1("f1", 10, 10);
  Kokkos::View<double*, TestMemorySpace> bcx1_view("bcx1_view", 40);
  Kokkos::View<double*, TestMemorySpace> bcx2_view("bcx2_view", 40);
  Kokkos::View<double*, TestMemorySpace> bcth1_view("bcth1_view", 40);
  Kokkos::View<double*, TestMemorySpace> bcth2_view("bcth2_view", 40);

  xset(10, x1_view, ex1_view, zero, one);
  xset(200, xtest_view, extest_view, zero, one);

  tset(10, t1_view, st1_view, ct1_view, zero, pi2);
  tset(200, ttest_view, stest_view, ctest_view, zero, pi2);

  ffset(10, ex1_view, st1_view, f1);

  auto x1 = MakeRank1View(x1_view);
  auto ex1 = MakeRank1View(ex1_view);
  auto t1 = MakeRank1View(t1_view);
  auto st1 = MakeRank1View(st1_view);
  auto ct1 = MakeRank1View(ct1_view);
  auto bcx1 = MakeRank1View(bcx1_view);
  auto bcx2 = MakeRank1View(bcx2_view);
  auto bcth1 = MakeRank1View(bcth1_view);
  auto bcth2 = MakeRank1View(bcth2_view);

  auto xtest = MakeRank1View(xtest_view);
  auto extest = MakeRank1View(extest_view);
  auto ttest = MakeRank1View(ttest_view);
  auto stest = MakeRank1View(stest_view);

  dotest2(x1, ex1, 10, t1, st1, ct1, 10, f1, bcx1, bcx2, bcth1, bcth2, xtest,
          extest, ttest, stest, 200, res_2d);
}

TEST_CASE("test_cubic_spline_interpolator")
{
  Kokkos::print_configuration(std::cout);
  std::cout << "DefaultExecutionSpace: "
            << Kokkos::DefaultExecutionSpace::name() << std::endl;
  Kokkos::Timer timer;
  {
    Kokkos::View<double*, HostMemorySpace> res_1d("res_1d", 6);
    double gt_1d[6] = {6.7572E-04, 6.6529E-04, 6.8669E-04,
                       6.7622E-04, 6.7572E-04, 6.6529E-04};
    pspltest1(res_1d);
    for (int i = 0; i < 6; ++i) {
      REQUIRE(are_equal(res_1d(i), gt_1d[i]));
    }

    Kokkos::View<double*, HostMemorySpace> res_1d_long("res_1d_long", 6);
    pspltest1_long_periodic(res_1d_long);
    REQUIRE(res_1d_long(2) < 1.0e-3);
    REQUIRE(res_1d_long(3) < 1.0e-3);
    REQUIRE(res_1d_long(4) < 1.0e-3);
    REQUIRE(res_1d_long(5) < 1.0e-3);

    Kokkos::View<double*, HostMemorySpace> res_2d("res_2d", 4);
    double gt_2d[4] = {1.8312E-03, 6.7151E-04, 1.8312E-03, 6.7151E-04};
    pspltest2(res_2d);
    for (int i = 0; i < 4; ++i) {
      REQUIRE(are_equal(res_2d(i), gt_2d[i]));
    }
  }

  double time = timer.seconds();
  std::cout << "Total time for spline tests: " << time << " seconds"
            << std::endl;
}

namespace
{

/// Copy an n x ncol evaluate output to the host as point-major rows
std::vector<double> to_row_major(Kokkos::View<double**, TestMemorySpace> out)
{
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, out);
  std::vector<double> rows;
  for (size_t i = 0; i < h.extent(0); ++i)
    for (size_t j = 0; j < h.extent(1); ++j)
      rows.push_back(h(i, j));
  return rows;
}

/// Evaluate a 1D interpolator with a selector and return the n x ncol output
template <typename Interpolator>
std::vector<double> eval_1d(Interpolator& interp, const std::vector<LO>& sel,
                            const std::vector<double>& xq, int ncol,
                            double fill = 0.0)
{
  auto xd = CreateDeviceRank1View(xq);
  Kokkos::View<double**, TestMemorySpace> out("out", xq.size(), ncol);
  Kokkos::deep_copy(out, fill);
  interp.evaluate(CreateDeviceRank1View<LO>(sel), MakeRank1View(xd),
                  MakeRank2View(out));
  return to_row_major(out);
}

/// Evaluate a 2D interpolator with a selector and return the n x ncol output
template <typename Interpolator>
std::vector<double> eval_2d(Interpolator& interp, const std::vector<LO>& sel,
                            const std::vector<double>& xq,
                            const std::vector<double>& yq, int ncol,
                            double fill = 0.0)
{
  auto xd = CreateDeviceRank1View(xq);
  auto yd = CreateDeviceRank1View(yq);
  Kokkos::View<double**, TestMemorySpace> out("out", xq.size(), ncol);
  Kokkos::deep_copy(out, fill);
  interp.evaluate(CreateDeviceRank1View<LO>(sel), MakeRank1View(xd),
                  MakeRank1View(yd), MakeRank2View(out));
  return to_row_major(out);
}

/// Uniform grid of n points on [a, b]
std::vector<double> linspace(double a, double b, int n)
{
  std::vector<double> x(n);
  for (int i = 0; i < n; ++i)
    x[i] = a + (b - a) * i / (n - 1);
  return x;
}

/// Sample f on the tensor grid with x as the slow index
template <typename F>
std::vector<double> sample_2d(const std::vector<double>& x,
                              const std::vector<double>& y, F f)
{
  std::vector<double> v;
  for (double xi : x)
    for (double yj : y)
      v.push_back(f(xi, yj));
  return v;
}

} // namespace

TEMPLATE_TEST_CASE("two-point splines keep both endpoint derivatives",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  // f = x^3 on [1, 3]: any pair of first/second derivative end conditions
  // determines the cubic exactly
  auto x = CreateDeviceRank1View<double>({1.0, 3.0});
  auto f = CreateDeviceRank1View<double>({1.0, 27.0});
  const auto d1 = BoundaryCondition::FIRST_DERIVATIVE_MATCH;
  const auto d2 = BoundaryCondition::SECOND_DERIVATIVE_MATCH;
  auto [lo, lo_data, hi, hi_data] = GENERATE_COPY(
    std::make_tuple(d1, 3.0, d1, 27.0), std::make_tuple(d1, 3.0, d2, 18.0),
    std::make_tuple(d2, 6.0, d1, 27.0), std::make_tuple(d2, 6.0, d2, 18.0));
  TestType interp(MakeRank1View(x), MakeRank1View(f), lo, lo_data, hi, hi_data);
  auto out = eval_1d(interp, {1, 1, 1}, {1.5, 2.0, 2.5}, 3);
  for (int i = 0; i < 3; ++i) {
    const double xi = 1.5 + 0.5 * i;
    CHECK(out[3 * i] == Catch::Approx(xi * xi * xi));
    CHECK(out[3 * i + 1] == Catch::Approx(3 * xi * xi));
    CHECK(out[3 * i + 2] == Catch::Approx(6 * xi));
  }
}

TEMPLATE_TEST_CASE("non-zero boundary data on the lower side only is applied",
                   "[spline]",
                   (ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactBiCubicSplineInterpolator<double, TestMemorySpace>))
{
  // f = g(x) (y - y^2/2) has df/dy = g(x) at y = 0 and 0 at y = 1
  auto g = [](double x) { return 1 + x + x * x + x * x * x; };
  auto f = [&](double x, double y) { return g(x) * (y - 0.5 * y * y); };
  auto xh = linspace(0, 1, 5);
  auto yh = linspace(0, 1, 5);
  std::vector<double> lo_h;
  for (double xi : xh)
    lo_h.push_back(g(xi));
  auto x = CreateDeviceRank1View(xh);
  auto y = CreateDeviceRank1View(yh);
  auto v = CreateDeviceRank1View(sample_2d(xh, yh, f));
  auto lo = CreateDeviceRank1View(lo_h);
  auto hi = CreateDeviceRank1View(std::vector<double>(xh.size(), 0.0));
  auto none = CreateDeviceRank1View<double>({});
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  const auto d1 = BoundaryCondition::FIRST_DERIVATIVE_MATCH;
  TestType interp(MakeRank1View(x), MakeRank1View(y), MakeRank1View(v), nak,
                  MakeRank1View(none), nak, MakeRank1View(none), d1,
                  MakeRank1View(lo), d1, MakeRank1View(hi));
  std::vector<double> xq = {0.1, 0.45, 0.8};
  std::vector<double> yq = {0.05, 0.6, 0.95};
  auto out = eval_2d(interp, {1, 0, 0, 0, 0, 0}, xq, yq, 1);
  for (size_t i = 0; i < xq.size(); ++i)
    CHECK(out[i] == Catch::Approx(f(xq[i], yq[i])));
}

TEST_CASE("spline grids must be finite, ascending and non-degenerate",
          "[spline]")
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  auto grid = GENERATE_COPY(std::vector<double>{1.0},
                            std::vector<double>{1.0, 1.0, 1.0, 1.0},
                            std::vector<double>{3.0, 2.0, 1.0, 0.0},
                            std::vector<double>{0.0, nan, 2.0, 3.0},
                            std::vector<double>{0.0, 1.0, 2.0, inf});
  auto x = CreateDeviceRank1View(grid);
  auto f = CreateDeviceRank1View(std::vector<double>(grid.size(), 1.0));
  using Interp = CompactCubicSplineInterpolator<double, TestMemorySpace>;
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  REQUIRE_THROWS_AS(Interp(MakeRank1View(x), MakeRank1View(f), nak, 0, nak, 0),
                    pcms_error);
}

TEST_CASE("nearly uniform spline grids are rejected", "[spline]")
{
  // the solver assumes one spacing, so a 5e-4 relative deviation would be
  // silently interpolated as if uniform
  auto x = CreateDeviceRank1View<double>({0.0, 1.0, 2.0005, 3.0});
  auto f = CreateDeviceRank1View<double>({0.0, 1.0, 8.0, 27.0});
  using Interp = ExplicitCubicSplineInterpolator<double, TestMemorySpace>;
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  REQUIRE_THROWS_AS(Interp(MakeRank1View(x), MakeRank1View(f), nak, 0, nak, 0),
                    pcms_error);
}

TEMPLATE_TEST_CASE("1D queries outside the grid leave the output untouched",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  auto xh = linspace(0, 1, 6);
  auto x = CreateDeviceRank1View(xh);
  auto f = CreateDeviceRank1View(xh);
  TestType interp(MakeRank1View(x), MakeRank1View(f),
                  BoundaryCondition::NOT_A_KNOT, 0,
                  BoundaryCondition::NOT_A_KNOT, 0);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto out = eval_1d(interp, {1, 1, 0}, {nan, -1.0e6, 1.0e6}, 2);
  for (double v : out)
    CHECK(v == 0.0);
}

TEMPLATE_TEST_CASE("2D queries outside the grid leave the output untouched",
                   "[spline]",
                   (ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactBiCubicSplineInterpolator<double, TestMemorySpace>))
{
  auto xh = linspace(0, 1, 6);
  auto x = CreateDeviceRank1View(xh);
  auto v = CreateDeviceRank1View(
    sample_2d(xh, xh, [](double a, double b) { return a + b; }));
  auto none = CreateDeviceRank1View<double>({});
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  TestType interp(MakeRank1View(x), MakeRank1View(x), MakeRank1View(v), nak,
                  MakeRank1View(none), nak, MakeRank1View(none), nak,
                  MakeRank1View(none), nak, MakeRank1View(none));
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto out = eval_2d(interp, {1, 1, 0, 0, 0, 0}, {nan, 0.5, -1.0e6, 0.5},
                     {0.5, nan, 0.5, 1.0e6}, 2);
  for (double val : out)
    CHECK(val == 0.0);
}

namespace
{

/// d^k/dx^k of the cubic with coefficients c
double cubic(const double (&c)[4], double x, int k)
{
  const double fall[4][4] = {
    {1, 1, 1, 1}, {0, 1, 2, 3}, {0, 0, 2, 6}, {0, 0, 0, 6}};
  double sum = 0;
  for (int j = k; j < 4; ++j)
    sum += c[j] * fall[k][j] * std::pow(x, j - k);
  return sum;
}

double px(double x, int k)
{
  const double c[4] = {1, 2, -1, 0.5};
  return cubic(c, x, k);
}

double py(double y, int k)
{
  const double c[4] = {2, -1, 0.3, 0.2};
  return cubic(c, y, k);
}

} // namespace

TEMPLATE_TEST_CASE("1D evaluate without a selector writes every value",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  // the overload builds its own value-only selector on the device; NaN marks
  // outputs it never wrote
  auto xh = linspace(-1, 2, 6);
  std::vector<double> fh;
  for (double x : xh)
    fh.push_back(px(x, 0));
  auto x = CreateDeviceRank1View(xh);
  auto f = CreateDeviceRank1View(fh);
  TestType interp(MakeRank1View(x), MakeRank1View(f),
                  BoundaryCondition::NOT_A_KNOT, 0,
                  BoundaryCondition::NOT_A_KNOT, 0);
  const std::vector<double> xq = {-1.0, 0.3, 1.7, 2.0};
  auto xd = CreateDeviceRank1View(xq);
  Kokkos::View<double**, TestMemorySpace> out("out", xq.size(), 1);
  Kokkos::deep_copy(out, std::numeric_limits<double>::quiet_NaN());
  interp.evaluate(MakeRank1View(xd), MakeRank2View(out));
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, out);
  for (size_t i = 0; i < xq.size(); ++i) {
    INFO("point " << i);
    CHECK(h(i, 0) == Catch::Approx(px(xq[i], 0)));
  }
}

TEMPLATE_TEST_CASE("2D evaluate without a selector writes every value",
                   "[spline]",
                   (ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactBiCubicSplineInterpolator<double, TestMemorySpace>))
{
  auto xh = linspace(0, 2, 6);
  auto yh = linspace(-1, 1, 5);
  auto x = CreateDeviceRank1View(xh);
  auto y = CreateDeviceRank1View(yh);
  auto v = CreateDeviceRank1View(
    sample_2d(xh, yh, [](double a, double b) { return px(a, 0) * py(b, 0); }));
  auto none = CreateDeviceRank1View<double>({});
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  TestType interp(MakeRank1View(x), MakeRank1View(y), MakeRank1View(v), nak,
                  MakeRank1View(none), nak, MakeRank1View(none), nak,
                  MakeRank1View(none), nak, MakeRank1View(none));
  const std::vector<double> xq = {0.0, 0.7, 1.3, 2.0};
  const std::vector<double> yq = {-1.0, 0.1, -0.45, 1.0};
  auto xd = CreateDeviceRank1View(xq);
  auto yd = CreateDeviceRank1View(yq);
  Kokkos::View<double**, TestMemorySpace> out("out", xq.size(), 1);
  Kokkos::deep_copy(out, std::numeric_limits<double>::quiet_NaN());
  interp.evaluate(MakeRank1View(xd), MakeRank1View(yd), MakeRank2View(out));
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, out);
  for (size_t i = 0; i < xq.size(); ++i) {
    INFO("point " << i);
    CHECK(h(i, 0) == Catch::Approx(px(xq[i], 0) * py(yq[i], 0)).margin(1e-9));
  }
}

TEMPLATE_TEST_CASE("1D splines build and evaluate in host memory", "[spline]",
                   (ExplicitCubicSplineInterpolator<double, Kokkos::HostSpace>),
                   (CompactCubicSplineInterpolator<double, Kokkos::HostSpace>))
{
  // every kernel must run in the interpolator's own space, even when the
  // default execution space is a device
  auto xh = linspace(-1, 2, 6);
  std::vector<double> fh;
  for (double x : xh)
    fh.push_back(px(x, 0));
  std::vector<double> xq = {-1.0, 0.3, 1.7, 2.0};
  using H = Kokkos::HostSpace;
  Kokkos::View<double**, H> out("out", xq.size(), 1);
  Kokkos::deep_copy(out, std::numeric_limits<double>::quiet_NaN());
  TestType interp(make_array_view(xh), make_array_view(fh),
                  BoundaryCondition::NOT_A_KNOT, 0,
                  BoundaryCondition::NOT_A_KNOT, 0);
  interp.evaluate(make_array_view(xq), MakeRank2View(out));
  Kokkos::fence();
  for (size_t i = 0; i < xq.size(); ++i) {
    INFO("point " << i);
    CHECK(out(i, 0) == Catch::Approx(px(xq[i], 0)));
  }
}

TEMPLATE_TEST_CASE(
  "2D splines build and evaluate in host memory", "[spline]",
  (ExplicitBiCubicSplineInterpolator<double, Kokkos::HostSpace>),
  (CompactBiCubicSplineInterpolator<double, Kokkos::HostSpace>))
{
  auto xh = linspace(0, 2, 6);
  auto yh = linspace(-1, 1, 5);
  auto v =
    sample_2d(xh, yh, [](double a, double b) { return px(a, 0) * py(b, 0); });
  std::vector<double> xq = {0.0, 0.7, 1.3, 2.0};
  std::vector<double> yq = {-1.0, 0.1, -0.45, 1.0};
  using H = Kokkos::HostSpace;
  Kokkos::View<double**, H> out("out", xq.size(), 1);
  Kokkos::deep_copy(out, std::numeric_limits<double>::quiet_NaN());
  Rank1View<double, H> none(nullptr, 0);
  const auto nak = BoundaryCondition::NOT_A_KNOT;
  TestType interp(make_array_view(xh), make_array_view(yh), make_array_view(v),
                  nak, none, nak, none, nak, none, nak, none);
  interp.evaluate(make_array_view(xq), make_array_view(yq), MakeRank2View(out));
  Kokkos::fence();
  for (size_t i = 0; i < xq.size(); ++i) {
    INFO("point " << i);
    CHECK(out(i, 0) == Catch::Approx(px(xq[i], 0) * py(yq[i], 0)).margin(1e-9));
  }
}

namespace
{

/// Divided difference f[x_i, .., x_{i+k}]
double divided_difference(const std::vector<double>& x,
                          const std::vector<double>& f, int i, int k)
{
  if (k == 0)
    return f[i];
  return (divided_difference(x, f, i + 1, k - 1) -
          divided_difference(x, f, i, k - 1)) /
         (x[i + k] - x[i]);
}

/// The conditions the solver applies on n points: two points take the first
/// divided difference, three points cannot use the third one and take the
/// parabola for not-a-knot at both ends
std::pair<int, int> effective_conditions(int n, int lo, int hi)
{
  if (n == 2) {
    lo = (lo <= 0 || lo > 5) ? 5 : lo;
    hi = (hi <= 0 || hi > 5) ? 5 : hi;
  }
  if (n == 3) {
    lo = lo == 7 ? 6 : lo;
    hi = hi == 7 ? 6 : hi;
    if (lo == 0 && hi == 0)
      hi = 7;
  }
  return {lo, hi};
}

/// End data that a product of the test cubics satisfies for condition bc
double end_data(BoundaryCondition bc, double (*p)(double, int), double x)
{
  switch (bc) {
    case BoundaryCondition::FIRST_DERIVATIVE_MATCH: return p(x, 1);
    case BoundaryCondition::SECOND_DERIVATIVE_MATCH: return p(x, 2);
    default: return 0;
  }
}

} // namespace

TEMPLATE_TEST_CASE("every end condition holds on grids of 2 to 8 points",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  const double lo_data = 0.7, hi_data = -0.4;
  for (int n = 2; n <= 8; ++n) {
    for (int lo = -1; lo <= 7; ++lo) {
      for (int hi = lo == -1 ? -1 : 0; hi <= (lo == -1 ? -1 : 7); ++hi) {
        INFO("n " << n << " lo " << lo << " hi " << hi);
        auto xh = linspace(0.5, 2.0, n);
        const double h = xh[1] - xh[0];
        std::vector<double> fh;
        for (double x : xh)
          fh.push_back(std::sin(1.3 * x) + 0.2 * x * x);
        if (lo == -1)
          fh.back() = fh.front();
        auto x = CreateDeviceRank1View(xh);
        auto f = CreateDeviceRank1View(fh);
        TestType s(MakeRank1View(x), MakeRank1View(f),
                   static_cast<BoundaryCondition>(lo), lo_data,
                   static_cast<BoundaryCondition>(hi), hi_data);

        // s', s'', s''' at both ends and in the two cells next to each end
        auto mid = [&](int cell) { return xh[cell] + h / 2; };
        const std::vector<double> pts = {xh[0],
                                         mid(0),
                                         mid(std::min(1, n - 2)),
                                         mid(std::max(n - 3, 0)),
                                         mid(n - 2),
                                         xh[n - 1]};
        auto d12 = eval_1d(s, {0, 1, 1}, pts, 2);
        auto d3 = eval_1d(s, {3, 0, 0}, pts, 1);
        auto at = [&](int point, int deriv) {
          return deriv == 3 ? d3[point] : d12[2 * point + deriv - 1];
        };

        auto check_end = [&](int side, int bc) {
          const int end = side == 0 ? 0 : 5;
          const double d1 = at(end, 1), d2 = at(end, 2);
          const double d3 = at(side == 0 ? 1 : 4, 3);
          const double d3_next = at(side == 0 ? 2 : 3, 3);
          const double data = side == 0 ? lo_data : hi_data;
          auto dd = [&](int k) {
            return divided_difference(xh, fh, side == 0 ? 0 : n - 1 - k, k);
          };
          auto approx = [](double v) {
            return Catch::Approx(v).margin(1e-9).epsilon(1e-9);
          };
          switch (bc) {
            case 0: CHECK(d3 == approx(d3_next)); break;
            case 1: CHECK(d1 == approx(data)); break;
            case 2: CHECK(d2 == approx(data)); break;
            case 3: CHECK(d1 == approx(0)); break;
            case 4: CHECK(d2 == approx(0)); break;
            case 5: CHECK(d1 == approx(dd(1))); break;
            case 6: CHECK(d2 == approx(2 * dd(2))); break;
            case 7: CHECK(d3 == approx(n < 4 ? 0 : 6 * dd(3))); break;
            default: FAIL("unexpected end condition " << bc);
          }
        };
        if (lo == -1) {
          CHECK(at(0, 1) == Catch::Approx(at(5, 1)).margin(1e-9));
          CHECK(at(0, 2) == Catch::Approx(at(5, 2)).margin(1e-9));
        } else {
          const auto [elo, ehi] = effective_conditions(n, lo, hi);
          check_end(0, elo);
          check_end(1, ehi);
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE("1D splines reproduce cubics for every exact end condition",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  using BC = BoundaryCondition;
  const std::vector<BC> conditions = {
    BC::NOT_A_KNOT, BC::FIRST_DERIVATIVE_MATCH, BC::SECOND_DERIVATIVE_MATCH,
    BC::THIRD_DIV_DIFF_MATCH};
  for (int n = 4; n <= 8; ++n) {
    for (BC lo : conditions) {
      for (BC hi : conditions) {
        INFO("n " << n << " lo " << static_cast<int>(lo) << " hi "
                  << static_cast<int>(hi));
        auto xh = linspace(-1, 2, n);
        std::vector<double> fh;
        for (double x : xh)
          fh.push_back(px(x, 0));
        auto x = CreateDeviceRank1View(xh);
        auto f = CreateDeviceRank1View(fh);
        TestType s(MakeRank1View(x), MakeRank1View(f), lo,
                   end_data(lo, px, xh.front()), hi,
                   end_data(hi, px, xh.back()));
        const std::vector<double> xq = {-1.0, -0.2, 0.9, 1.6, 2.0};
        auto out = eval_1d(s, {1, 1, 1}, xq, 3);
        auto d3 = eval_1d(s, {3, 0, 0}, xq, 1);
        for (size_t q = 0; q < xq.size(); ++q) {
          for (int j = 0; j < 3; ++j)
            CHECK(out[3 * q + j] == Catch::Approx(px(xq[q], j)).margin(1e-10));
          CHECK(d3[q] == Catch::Approx(px(xq[q], 3)).margin(1e-9));
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE("periodic 1D splines wrap queries outside the period",
                   "[spline]",
                   (ExplicitCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactCubicSplineInterpolator<double, TestMemorySpace>))
{
  const double pi2 = 2 * M_PI;
  auto xh = linspace(0, pi2, 16);
  std::vector<double> fh;
  for (double a : xh)
    fh.push_back(std::sin(a) + 0.5 * std::cos(2 * a));
  fh.back() = fh.front();
  auto x = CreateDeviceRank1View(xh);
  auto f = CreateDeviceRank1View(fh);
  const auto per = BoundaryCondition::PERIODIC;
  TestType s(MakeRank1View(x), MakeRank1View(f), per, 0, per, 0);
  const std::vector<double> outside = {-0.1, pi2 + 0.1, -pi2 + 0.3,
                                       2 * pi2 + 1.0};
  const std::vector<double> inside = {pi2 - 0.1, 0.1, 0.3, 1.0};
  auto o = eval_1d(s, {1, 1, 1}, outside, 3, -999.0);
  auto i = eval_1d(s, {1, 1, 1}, inside, 3);
  for (size_t k = 0; k < o.size(); ++k)
    CHECK(o[k] == Catch::Approx(i[k]).margin(1e-10));
}

TEMPLATE_TEST_CASE("periodic 2D splines wrap queries outside the period",
                   "[spline]",
                   (ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>),
                   (CompactBiCubicSplineInterpolator<double, TestMemorySpace>))
{
  const double pi2 = 2 * M_PI;
  auto xh = linspace(0, pi2, 12);
  auto yh = linspace(0, pi2, 10);
  auto x = CreateDeviceRank1View(xh);
  auto y = CreateDeviceRank1View(yh);
  // sin and cos of the grid ends agree to round-off, so the samples are
  // periodic in both directions
  auto v = CreateDeviceRank1View(sample_2d(
    xh, yh, [](double a, double b) { return std::sin(a) * std::cos(b) + 1; }));
  auto none = CreateDeviceRank1View<double>({});
  const auto per = BoundaryCondition::PERIODIC;
  TestType s(MakeRank1View(x), MakeRank1View(y), MakeRank1View(v), per,
             MakeRank1View(none), per, MakeRank1View(none), per,
             MakeRank1View(none), per, MakeRank1View(none));
  const std::vector<double> xo = {-0.1, pi2 + 0.2, 1.0, -0.4};
  const std::vector<double> yo = {1.0, -0.3, pi2 + 0.5, pi2 + 0.1};
  const std::vector<double> xi = {pi2 - 0.1, 0.2, 1.0, pi2 - 0.4};
  const std::vector<double> yi = {1.0, pi2 - 0.3, 0.5, 0.1};
  auto o = eval_2d(s, {1, 1, 1, 1, 1, 1}, xo, yo, 6, -999.0);
  auto i = eval_2d(s, {1, 1, 1, 1, 1, 1}, xi, yi, 6);
  for (size_t k = 0; k < o.size(); ++k)
    CHECK(o[k] == Catch::Approx(i[k]).margin(1e-10));
}
