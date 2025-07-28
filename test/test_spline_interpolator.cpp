#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "mdspan/mdspan.hpp"
#include "pcms.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <pcms/interpolator/spline_interpolator.hpp>

using namespace pcms;

using TestMemorySpace = Kokkos::DefaultExecutionSpace;

bool are_equal(double a, double b, double tolerance = 1e-7) {
  return std::abs(a - b) < tolerance;
}

void createFlatGrid(Rank1View<double, TestMemorySpace> xvec,
                    Rank1View<double, TestMemorySpace> yvec,
                    Kokkos::View<double *, TestMemorySpace> X_flat,
                    Kokkos::View<double *, TestMemorySpace> Y_flat) {
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

void tset(int nth, Kokkos::View<double *, TestMemorySpace> th,
          Kokkos::View<double *, TestMemorySpace> sth,
          Kokkos::View<double *, TestMemorySpace> cth, double thmin,
          double thmax) {
  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, nth),
      KOKKOS_LAMBDA(const int ith) {
        th[ith] = thmin + static_cast<double>(ith) * (thmax - thmin) /
                              static_cast<double>(nth - 1);
        sth[ith] = 2.0 + std::sin(th[ith]);
        cth[ith] = std::cos(th[ith]);
      });
}

void xset(int nx, Kokkos::View<double *, TestMemorySpace> x,
          Kokkos::View<double *, TestMemorySpace> ex, double xmin,
          double xmax) {
  if (nx < 2)
    return; // avoid division by zero

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, nx), KOKKOS_LAMBDA(const int ix) {
        x[ix] = xmin + static_cast<double>(ix) * (xmax - xmin) /
                           static_cast<double>(nx - 1);
        ex[ix] = std::exp(2.0 * x[ix] - 1.0);
      });
}

void ffset(int num, Kokkos::View<double *, TestMemorySpace> xf,
           Kokkos::View<double *, TestMemorySpace> tf,
           Kokkos::View<double *, TestMemorySpace> f) {
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num, num}),
      KOKKOS_LAMBDA(const int i, const int j) {
        f(i * num + j) = xf[i] * tf[j];
      });
}

void bset(Rank1View<double, TestMemorySpace> fx, int nx,
          Rank1View<double, TestMemorySpace> fth, int nth,
          Rank1View<double, TestMemorySpace> bcx1,
          Rank1View<double, TestMemorySpace> bcx2,
          Rank1View<double, TestMemorySpace> bcth1,
          Rank1View<double, TestMemorySpace> bcth2) {
  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, nth),
      KOKKOS_LAMBDA(const int ith) {
        bcx1[ith] = 2.0 * fx(0) * fth(ith);      // df/dx at x(1)
        bcx2[ith] = 2.0 * fx(nx - 1) * fth(ith); // df/dx at x(nx)
      });

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, nx), KOKKOS_LAMBDA(const int ix) {
        bcth1(ix) = fx(ix); // df/dth at th = 0 (th[0])
        bcth2(ix) = fx(ix); // df/dth at th = 2π (th[nth - 1])
      });
}

void generate_gt_2d(Rank1View<double, TestMemorySpace> gt_splinv,
                    Rank1View<double, TestMemorySpace> fxtest,
                    Rank1View<double, TestMemorySpace> fthtest, int ntest) {
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ntest, ntest}),
      KOKKOS_LAMBDA(const int j, const int i) {
        gt_splinv(j * ntest + i) = fxtest(i) * fthtest(j);
      });
}

void compare(Rank2View<double, TestMemorySpace> splinv,
             Rank1View<double, TestMemorySpace> gt_splinv,
             Rank1View<double, HostMemorySpace> res, std::string label) {
  // Compare the output of splinv with gt_splinv
  double result_fdif = 0.0;
  double result_fdifr = 0.0;
  Kokkos::parallel_reduce(
      "compare_splinv", splinv.extent(0),
      KOKKOS_LAMBDA(const int i, double &max_fdif, double &max_fdifr) {
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
                Rank1View<double, HostMemorySpace> res, std::string label) {
  // Compare the output of splinv with gt_splinv
  double result_fdif = 0.0;
  double result_fdifr = 0.0;
  Kokkos::parallel_reduce(
      "compare_splinv", splinv.extent(0),
      KOKKOS_LAMBDA(const int i, double &max_fdif, double &max_fdifr) {
        if (splinv(i, 0) != 0.0) {
          double dif = Kokkos::abs(splinv(i, 0) - gt_splinv(i));
          max_fdif = Kokkos::max(max_fdif, dif);
          max_fdifr = Kokkos::max(max_fdifr,
                                  dif / (0.5 * (splinv(i, 0) + gt_splinv(i))));
        }
      },
      Kokkos::Max<double>(result_fdif), Kokkos::Max<double>(result_fdifr));
  res(0) = result_fdif;
  res(1) = result_fdifr;
  std::cout << "2d " << label << " max absolute difference: " << res(0)
            << ", relative difference: " << res(1) << std::endl;
}

void reset_2dspan(Rank2View<double, TestMemorySpace> span) {
  Kokkos::parallel_for(
      "zero_init",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                             {span.extent(0), span.extent(1)}),
      KOKKOS_LAMBDA(int i, int j) { span(i, j) = 0.0; });
}

void dotest1(int ns, Rank1View<double, TestMemorySpace> x,
             Rank1View<double, TestMemorySpace> f,
             Rank1View<double, TestMemorySpace> fd,
             Rank2View<double, TestMemorySpace> fspl,
             Rank2View<double, TestMemorySpace> fspp,
             Rank2View<double, TestMemorySpace> fs2, int nt,
             Rank1View<double, TestMemorySpace> xt,
             Rank1View<double, TestMemorySpace> ft,
             Rank2View<double, TestMemorySpace> xpkg,
             Rank2View<double, TestMemorySpace> testa1,
             Rank2View<double, TestMemorySpace> testa2,
             Rank2View<double, TestMemorySpace> testa3,
             Rank1View<double, TestMemorySpace> wk,
             Rank1View<double, TestMemorySpace> wk2,
             Kokkos::View<double *, HostMemorySpace> res_1d) {
  using std::data;
  using std::size;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, ns), KOKKOS_LAMBDA(const int i) {
        fspl(0, i) = f(i);
        fspp(0, i) = f(i);
        fs2(0, i) = f(i);
      });

  Kokkos::View<double *, TestMemorySpace> fspl4_view("fspl4_view", 40);
  auto fspl4 = Rank2View<double, TestMemorySpace>(fspl4_view.data(), 4, 10);

  // interpolator.cspline(x, ns, fspl, 1, 1, 1, 1, wk);
  ExplicitCubicSplineInterpolator<double, TestMemorySpace>
      explicit_interpolator(x, ns, fspl, 1, 1, 1, 1, wk);

  int ict_arr[3] = {1, 0, 0};
  Kokkos::View<int *, TestMemorySpace> ict("ict", 3);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, ict.size()),
      KOKKOS_LAMBDA(const int i) { ict(i) = ict_arr[i]; });

  Kokkos::View<double *, TestMemorySpace> splinv_view("splinv_view", 3000);
  auto splinv = Rank2View<double, TestMemorySpace>(splinv_view.data(), 1000, 3);
  reset_2dspan(splinv);
  explicit_interpolator.evaluate(ict, xt, splinv);

  auto res1 = Rank1View<double, HostMemorySpace>(res_1d.data(), 2);
  compare(splinv, ft, res1, "cspline");

  reset_2dspan(splinv);
  ExplicitCubicSplineInterpolator<double, TestMemorySpace>
      explicit_interpolator_periodic(x, ns, fspp, -1, 0, -1, 0, wk);
  explicit_interpolator_periodic.evaluate(ict, xt, splinv);
  auto res2 = Rank1View<double, HostMemorySpace>(res_1d.data() + 2, 2);
  compare(splinv, ft, res2, "cspline_periodic");

  reset_2dspan(splinv);
  CompactCubicSplineInterpolator<double, TestMemorySpace> compact_interpolator(
      x, ns, fs2, fspl4, 1, 1, 1, 1, wk2);
  compact_interpolator.evaluate(ict, xt, splinv);
  auto res3 = Rank1View<double, HostMemorySpace>(res_1d.data() + 4, 2);
  compare(splinv, ft, res3, "compact_cubic");
}

void pspltest1(Kokkos::View<double *, HostMemorySpace> res_1d) {
  const double pi2 = 6.28318530718;
  const double zero = 0.0;
  int inum = 10;

  // Local arrays
  Kokkos::View<double *, TestMemorySpace> zdum_view("zdum_view", 1000);
  Kokkos::View<double *, TestMemorySpace> testa1_view("testa1_view", 1000);
  Kokkos::View<double *, TestMemorySpace> testa2_view("testa2_view", 1000);
  Kokkos::View<double *, TestMemorySpace> testa3_view("testa3_view", 1000);
  Kokkos::View<double *, TestMemorySpace> xtest_view("xtest_view", 1000);
  Kokkos::View<double *, TestMemorySpace> ftest_view("ftest_view", 1000);
  Kokkos::View<double *, TestMemorySpace> x_view("x_view", 10);
  Kokkos::View<double *, TestMemorySpace> zcos_view("zcos_view", 10);
  Kokkos::View<double *, TestMemorySpace> z2sin_view("z2sin_view", 10);
  Kokkos::View<double *, TestMemorySpace> xpkg_view("xpkg_view", 10 * 4);
  Kokkos::View<double *, TestMemorySpace> fs_view("fs_view", 40);
  Kokkos::View<double *, TestMemorySpace> fsp_view("fsp_view", 40);
  Kokkos::View<double *, TestMemorySpace> fs2_view("fs2_view", 20);

  // Prepare test data
  tset(1000, xtest_view, ftest_view, zdum_view, zero - 0.1, pi2 + 0.1);
  tset(inum, x_view, z2sin_view, zcos_view, zero, pi2);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, inum),
      KOKKOS_LAMBDA(const int ix) { zcos_view(ix) = std::cos(x_view(ix)); });

  auto zdum = Rank1View<double, TestMemorySpace>(zdum_view.data(), 10);

  auto wk2_view =
      Kokkos::View<double *, TestMemorySpace>(zdum_view.data(), 1000);
  auto wk2 = Rank1View<double, TestMemorySpace>(wk2_view.data(), 10);

  auto testa1 = Rank2View<double, TestMemorySpace>(testa1_view.data(), 1000, 3);
  auto testa2 = Rank2View<double, TestMemorySpace>(testa2_view.data(), 1000, 3);
  auto testa3 = Rank2View<double, TestMemorySpace>(testa3_view.data(), 1000, 3);
  auto x = Rank1View<double, TestMemorySpace>(x_view.data(), x_view.size());
  auto zcos =
      Rank1View<double, TestMemorySpace>(zcos_view.data(), zcos_view.size());
  auto z2sin =
      Rank1View<double, TestMemorySpace>(z2sin_view.data(), z2sin_view.size());
  auto fs = Rank2View<double, TestMemorySpace>(fs_view.data(), 4, 10);
  auto fsp = Rank2View<double, TestMemorySpace>(fsp_view.data(), 4, 10);
  auto fs2 = Rank2View<double, TestMemorySpace>(fs2_view.data(), 2, 10);
  auto xpkg = Rank2View<double, TestMemorySpace>(xpkg_view.data(), 10, 4);
  auto xtest = Rank1View<double, TestMemorySpace>(xtest_view.data(), 1000);
  auto ftest = Rank1View<double, TestMemorySpace>(ftest_view.data(), 1000);

  // Call test function
  dotest1(inum, x, z2sin, zcos, fs, fsp, fs2, 1000, xtest, ftest, xpkg, testa1,
          testa2, testa3, zdum, wk2, res_1d);
}

void dotest2(Rank1View<double, TestMemorySpace> x,
             Rank1View<double, TestMemorySpace> fx, int nx,
             Rank1View<double, TestMemorySpace> th,
             Rank1View<double, TestMemorySpace> fth,
             Rank1View<double, TestMemorySpace> dfth, int nth,
             Rank4View<double, TestMemorySpace> f,
             Rank3View<double, TestMemorySpace> fh,
             Rank2View<double, TestMemorySpace> flin,
             Rank1View<double, TestMemorySpace> bcx1,
             Rank1View<double, TestMemorySpace> bcx2,
             Rank1View<double, TestMemorySpace> bcth1,
             Rank1View<double, TestMemorySpace> bcth2,
             Rank1View<double, TestMemorySpace> xtest,
             Rank1View<double, TestMemorySpace> fxtest,
             Rank1View<double, TestMemorySpace> thtest,
             Rank1View<double, TestMemorySpace> fthtest, int ntest,
             Kokkos::View<double *, HostMemorySpace> res_2d) {
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nth, nx}),
      KOKKOS_LAMBDA(const int ith, const int ix) {
        flin(ix, ith) = f(0, 0, ix, ith);          // f
        fh(0, ix, ith) = f(0, 0, ix, ith);         // f
        fh(1, ix, ith) = 2.0 * f(0, 0, ix, ith);   // df/dx
        fh(2, ix, ith) = fx(ix) * dfth(ith);       // df/dy
        fh(3, ix, ith) = 2.0 * fx(ix) * dfth(ith); // d2f/dxdy
      });

  int nbc = 1;
  Kokkos::View<double *, TestMemorySpace> wk_vec("wk_vec", 900);
  auto wk = Rank1View<double, TestMemorySpace>(wk_vec.data(), 900);
  bset(fx, nx, fth, nth, bcx1, bcx2, bcth1, bcth2);

  ExplicitBiCubicSplineInterpolator<double, TestMemorySpace>
      explicit_interpolator(x, nx, th, nth, f, nbc, bcx1, nbc, bcx2, nbc, bcth1,
                            nbc, bcth2, wk);

  Kokkos::View<double *, TestMemorySpace> splinv_view("splinv_view", 400000);
  auto splinv =
      Rank2View<double, TestMemorySpace>(splinv_view.data(), 40000, 10);
  reset_2dspan(splinv);
  int isel_arr[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Kokkos::View<int *, TestMemorySpace> isel("isel", 10);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<TestMemorySpace>(0, isel.size()),
      KOKKOS_LAMBDA(const int i) { isel(i) = isel_arr[i]; });

  Kokkos::View<double *, TestMemorySpace> xtest_grid_view("xtest_grid_view",
                                                          40000);
  Kokkos::View<double *, TestMemorySpace> thtest_grid_view("thtest_grid_view",
                                                           40000);
  createFlatGrid(xtest, thtest, xtest_grid_view, thtest_grid_view);
  auto xtest_grid =
      Rank1View<double, TestMemorySpace>(xtest_grid_view.data(), 40000);
  auto thtest_grid =
      Rank1View<double, TestMemorySpace>(thtest_grid_view.data(), 40000);
  Kokkos::View<double *, TestMemorySpace> splinv_gt_view("splinv_view", 40000);
  auto splinv_gt =
      Rank1View<double, TestMemorySpace>(splinv_gt_view.data(), 40000);
  generate_gt_2d(splinv_gt, fxtest, fthtest, ntest);

  // TODO: evaluate function raise memeory acess violation
  explicit_interpolator.evaluate(isel, xtest_grid, thtest_grid, splinv);

  auto res1 = Rank1View<double, HostMemorySpace>(res_2d.data(), 2);
  compare_2d(splinv, splinv_gt, res1, "bcspline");

  reset_2dspan(splinv);
  CompactBiCubicSplineInterpolator<double, TestMemorySpace>
      compact_interpolator(x, nx, th, nth, fh, nbc, bcx1, nbc, bcx2, nbc, bcth1,
                           nbc, bcth2, wk);
  compact_interpolator.evaluate(isel, xtest_grid, thtest_grid, splinv);
  auto res2 = Rank1View<double, HostMemorySpace>(res_2d.data() + 2, 2);
  compare_2d(splinv, splinv_gt, res2, "mkbicub");
}

void pspltest2(Kokkos::View<double *, HostMemorySpace> res_2d) {
  const double pi2 = 6.28318530718;
  const double zero = 0.0;
  const double one = 1.0;

  Kokkos::View<double *, TestMemorySpace> x1_view("x1_view", 10);
  Kokkos::View<double *, TestMemorySpace> ex1_view("ex1_view", 10);
  Kokkos::View<double *, TestMemorySpace> t1_view("t1_view", 10);
  Kokkos::View<double *, TestMemorySpace> st1_view("st1_view", 10);
  Kokkos::View<double *, TestMemorySpace> ct1_view("ct1_view", 10);
  Kokkos::View<double *, TestMemorySpace> xtest_view("xtest_view", 200);
  Kokkos::View<double *, TestMemorySpace> extest_view("extest_view", 200);
  Kokkos::View<double *, TestMemorySpace> ttest_view("ttest_view", 200);
  Kokkos::View<double *, TestMemorySpace> stest_view("stest_view", 200);
  Kokkos::View<double *, TestMemorySpace> ctest_view("ctest_view", 200);
  Kokkos::View<double *, TestMemorySpace> f1_view("f1_view", 1600);
  Kokkos::View<double *, TestMemorySpace> f2_view("f2_view", 6400);
  Kokkos::View<double *, TestMemorySpace> f4_view("f4_view", 4 * 4 * 40 * 40);
  Kokkos::View<double *, TestMemorySpace> bcx1_view("bcx1_view", 40);
  Kokkos::View<double *, TestMemorySpace> bcx2_view("bcx2_view", 40);
  Kokkos::View<double *, TestMemorySpace> bcth1_view("bcth1_view", 40);
  Kokkos::View<double *, TestMemorySpace> bcth2_view("bcth2_view", 40);
  Kokkos::View<double *, TestMemorySpace> fh_view("fh_view", 400);
  Kokkos::View<double *, TestMemorySpace> flin_view("flin_view", 100);

  xset(10, x1_view, ex1_view, zero, one);
  // xset(20, x2_view, ex2_view, zero, one);
  // xset(40, x4_view, ex4_view, zero, one);
  xset(200, xtest_view, extest_view, zero, one);

  tset(10, t1_view, st1_view, ct1_view, zero, pi2);
  // tset(20, t2_view, st2_view, ct2_view, zero, pi2);
  // tset(40, t4_view, st4_view, ct4_view, zero, pi2);
  tset(200, ttest_view, stest_view, ctest_view, zero, pi2);

  ffset(10, ex1_view, st1_view, f1_view);
  // ffset(20, ex2_view, st2_view, f2_view);
  // ffset(40, ex4_view, st4_view, f4_view);

  auto x1 = Rank1View<double, TestMemorySpace>(x1_view.data(), x1_view.size());
  auto ex1 =
      Rank1View<double, TestMemorySpace>(ex1_view.data(), ex1_view.size());
  auto t1 = Rank1View<double, TestMemorySpace>(t1_view.data(), t1_view.size());
  auto st1 =
      Rank1View<double, TestMemorySpace>(st1_view.data(), st1_view.size());
  auto ct1 =
      Rank1View<double, TestMemorySpace>(ct1_view.data(), ct1_view.size());
  auto f1 = Rank4View<double, TestMemorySpace>(f1_view.data(),
                                               Kokkos::extents{4, 4, 10, 10});
  auto flin = Rank2View<double, TestMemorySpace>(flin_view.data(), 10, 10);

  auto fh = Rank3View<double, TestMemorySpace>(fh_view.data(), 4, 10, 10);

  auto bcx1 =
      Rank1View<double, TestMemorySpace>(bcx1_view.data(), bcx1_view.size());
  auto bcx2 =
      Rank1View<double, TestMemorySpace>(bcx2_view.data(), bcx2_view.size());
  auto bcth1 =
      Rank1View<double, TestMemorySpace>(bcth1_view.data(), bcth1_view.size());
  auto bcth2 =
      Rank1View<double, TestMemorySpace>(bcth2_view.data(), bcth2_view.size());

  auto xtest =
      Rank1View<double, TestMemorySpace>(xtest_view.data(), xtest_view.size());
  auto extest = Rank1View<double, TestMemorySpace>(extest_view.data(),
                                                   extest_view.size());
  auto ttest =
      Rank1View<double, TestMemorySpace>(ttest_view.data(), ttest_view.size());
  auto stest =
      Rank1View<double, TestMemorySpace>(stest_view.data(), stest_view.size());

  dotest2(x1, ex1, 10, t1, st1, ct1, 10, f1, fh, flin, bcx1, bcx2, bcth1, bcth2,
          xtest, extest, ttest, stest, 200, res_2d);

  // dotest2(x2,ex2,20,t2,st2,ct2,20,f2,fh,flin,
  //     bcx1,bcx2,bcth1,bcth2,
  //     xtest,extest,ttest,stest,200);

  // dotest2(x4,ex4,40,t4,st4,ct4,40,f4,fh,flin,
  //     bcx1,bcx2,bcth1,bcth2,
  //     xtest,extest,ttest,stest,200);
}

TEST_CASE("test_cubic_spline_interpolator") {
  Kokkos::print_configuration(std::cout);
  std::cout << "DefaultExecutionSpace: "
            << Kokkos::DefaultExecutionSpace::name() << std::endl;
  Kokkos::Timer timer;
  {
    Kokkos::View<double *, HostMemorySpace> res_1d("res_1d", 6);
    double gt_1d[6] =
    {6.7572E-04, 6.6529E-04, 6.8669E-04, 6.7622E-04, 6.7572E-04, 6.6529E-04};
    pspltest1(res_1d);
    for (int i = 0; i < 6; ++i) {
      REQUIRE(are_equal(res_1d(i), gt_1d[i]));
    }

    Kokkos::View<double *, HostMemorySpace> res_2d("res_2d", 4);
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