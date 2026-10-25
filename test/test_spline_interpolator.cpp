#include "mdspan/mdspan.hpp"
#include "pcms.h"
#include <Kokkos_Core.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
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
