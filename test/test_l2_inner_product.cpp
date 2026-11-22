#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include "field_test_utils.hpp"
#include <cmath>

using Catch::Approx;

TEST_CASE("L2Integral is exact for polynomials of the space",
          "[l2][quadrature]")
{
  Omega_h::Library lib;
  const int dim = GENERATE(2, 3);
  auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 3)
                       : pcms::test::BuildUnitCube(lib, 2);
  const auto space = pcms::test::MakeP1Space(mesh);
  auto u = space->CreateFunction<pcms::Real>();
  auto v = space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    u, [](pcms::Real x, pcms::Real y, pcms::Real* out) { out[0] = x; });
  pcms::test::SetFieldComponents(
    v, [](pcms::Real x, pcms::Real y, pcms::Real* out) { out[0] = x + y; });
  const pcms::L2Integral product(space, space);
  // int x (x + y) = 1/3 + 1/4 on the unit square and cube
  CHECK(product.Inner(u, v) == Approx(7.0 / 12.0).epsilon(1e-13));
  CHECK(product.Inner(u, u) == Approx(1.0 / 3.0).epsilon(1e-13));
  // |x - (x + y)|^2 = y^2
  CHECK(product.Distance(u, v) == Approx(std::sqrt(1.0 / 3.0)).epsilon(1e-13));
  const pcms::L2Integral integral(space);
  CHECK(integral.Apply(v) == Approx(1.0).epsilon(1e-13));
  CHECK(integral.Norm(u) == Approx(std::sqrt(1.0 / 3.0)).epsilon(1e-13));
}

TEST_CASE("L2Integral pairs a function with a source on another mesh",
          "[l2][quadrature]")
{
  Omega_h::Library lib;
  auto coarse = pcms::test::BuildSquareMesh(lib, 2);
  auto fine = pcms::test::BuildSquareMesh(lib, 5);
  const auto coarse_space = pcms::test::MakeP1Space(coarse);
  const auto fine_space = pcms::test::MakeP1Space(fine);
  auto u = coarse_space->CreateFunction<pcms::Real>();
  auto v = fine_space->CreateFunction<pcms::Real>();
  const auto linear = [](pcms::Real x, pcms::Real y, pcms::Real* out) {
    out[0] = 1 + x - 2 * y;
  };
  pcms::test::SetFieldComponents(u, linear);
  pcms::test::SetFieldComponents(v, linear);
  const pcms::L2Integral product(coarse_space, fine_space);
  CHECK(product.Distance(u, v) < 1e-14);
  CHECK(product.Inner(u, v) == Approx(2.0 / 3.0).epsilon(1e-13)); // E[f^2]
  // a given quadrature must be on u's geometry
  CHECK_THROWS_AS(
    pcms::L2Integral(coarse_space, fine_space,
                     pcms::CreateGaussQuadrature(fine_space->GetGeometry(), 2)),
    pcms::pcms_error);
}

TEST_CASE("L2Integral::Distance of nearly equal fields keeps its accuracy",
          "[l2][quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto space = pcms::test::MakeP1Space(mesh);
  auto u = space->CreateFunction<pcms::Real>();
  auto v = space->CreateFunction<pcms::Real>();
  const pcms::Real offset = 1e4;
  const pcms::Real epsilon = 1e-6;
  pcms::test::SetFieldComponents(
    u,
    [=](pcms::Real x, pcms::Real y, pcms::Real* out) { out[0] = offset + x; });
  pcms::test::SetFieldComponents(
    v, [=](pcms::Real x, pcms::Real y, pcms::Real* out) {
      out[0] = offset + x + epsilon * y;
    });
  const pcms::L2Integral product(space, space);
  const double exact = epsilon * std::sqrt(1.0 / 3.0);
  CHECK(product.Distance(u, v) == Approx(exact).epsilon(1e-6));
  // the expansion would not: |u|^2 - 2 u.v + |v|^2 is a difference of numbers
  // near 1e8 with a result near 1e-13
  const double expanded =
    product.Inner(u, u) - 2 * product.Inner(u, v) + product.Inner(v, v);
  CHECK(std::abs(expanded - exact * exact) > 1e-3 * exact * exact);
}

TEST_CASE("L2Integral of vector fields and of P0", "[l2][quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto space = pcms::test::MakeP0Space(mesh);
  auto u = space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    u, [](pcms::Real, pcms::Real, pcms::Real* out) { out[0] = 3.0; });
  const pcms::L2Integral integral(space);
  CHECK(integral.Apply(u) == Approx(3.0).epsilon(1e-13));
  CHECK(integral.Norm(u) == Approx(3.0).epsilon(1e-13));
  const auto vector_space = pcms::CreateLagrangeSpace(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1,
    2, pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
  const pcms::L2Integral vector_integral(vector_space);
  auto w = vector_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*vector_space));
  pcms::test::SetFieldComponents(w,
                                 [](pcms::Real, pcms::Real, pcms::Real* out) {
                                   out[0] = 3.0;
                                   out[1] = 4.0;
                                 });
  CHECK(vector_integral.Norm(w) == Approx(5.0).epsilon(1e-13));
  CHECK_THROWS(vector_integral.Apply(w));
}
