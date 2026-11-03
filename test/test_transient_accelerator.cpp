#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/accelerator.hpp"

#include <cmath>
#include <vector>

using pcms::Real;
using Status = pcms::InterfaceAccelerator::UpdateStatus;

TEST_CASE("Aitken relaxation updates and resets the interface iterate",
          "[transient]")
{
  // Use an initial relaxation factor.
  pcms::AitkenRelaxation accelerator(0.5);
  std::vector<Real> next;

  accelerator.BeginWindow();
  REQUIRE(accelerator.LastStatus() == Status::Initial);

  // These vectors represent the interface iterate before and after one
  // fixed-point evaluation.
  const Real first_residual = accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);

  // The returned residual norm is ||{2,-4}||_2, and omega=0.5.
  REQUIRE(first_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(next[0] == Catch::Approx(1.0));
  REQUIRE(next[1] == Catch::Approx(-2.0));
  REQUIRE(accelerator.LastStatus() == Status::Initial);

  // A second residual lets Aitken compute omega=2/3 from iteration history.
  const Real second_residual = accelerator.Update(next, {1.5, -3.0}, next);
  REQUIRE(second_residual == Catch::Approx(std::sqrt(1.25)));
  REQUIRE(next[0] == Catch::Approx(4.0 / 3.0));
  REQUIRE(next[1] == Catch::Approx(-8.0 / 3.0));
  REQUIRE(accelerator.LastStatus() == Status::Updated);

  // The status is available through the base interface.
  const pcms::InterfaceAccelerator& base = accelerator;
  REQUIRE(base.LastStatus() == Status::Updated);

  // A new window must discard the history and reuse omega=0.5.
  accelerator.BeginWindow();
  REQUIRE(accelerator.LastStatus() == Status::Initial);
  const Real third_residual = accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);
  REQUIRE(third_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(next[0] == Catch::Approx(1.0));
  REQUIRE(next[1] == Catch::Approx(-2.0));
}

TEST_CASE("Aitken relaxation reports what happened to the relaxation factor",
          "[transient]")
{
  std::vector<Real> next;

  SECTION("the factor is kept when the residual does not change")
  {
    pcms::AitkenRelaxation accelerator(0.5);
    accelerator.BeginWindow();

    // Both iterations have the residual {2,-4}, so the change in residual is 0.
    (void)accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);
    (void)accelerator.Update({1.0, -2.0}, {3.0, -6.0}, next);

    // omega stays 0.5, so next = {1,-2} + 0.5 * {2,-4}.
    REQUIRE(accelerator.LastStatus() == Status::Fallback);
    REQUIRE(next[0] == Catch::Approx(2.0));
    REQUIRE(next[1] == Catch::Approx(-4.0));
  }

  SECTION("the factor is reset to the initial factor when it would be zero")
  {
    pcms::AitkenRelaxation accelerator(0.5);
    accelerator.BeginWindow();

    // The residuals {1,0} and {1,1} give a change {0,1} that is orthogonal to
    // the previous residual, so the new omega is 0 and is reset to 0.5.
    (void)accelerator.Update({0.0, 0.0}, {1.0, 0.0}, next);
    (void)accelerator.Update({0.0, 0.0}, {1.0, 1.0}, next);

    REQUIRE(accelerator.LastStatus() == Status::Fallback);
    REQUIRE(next[0] == Catch::Approx(0.5));
    REQUIRE(next[1] == Catch::Approx(0.5));
  }

  SECTION("the factor is clamped when it exceeds the limit")
  {
    // The unclamped omega is 2/3, which is above the limit of 0.6.
    pcms::AitkenRelaxation accelerator(0.5, 0.6);
    accelerator.BeginWindow();

    (void)accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);
    (void)accelerator.Update({1.0, -2.0}, {1.5, -3.0}, next);

    // next = {1,-2} + 0.6 * {0.5,-1}.
    REQUIRE(accelerator.LastStatus() == Status::Limited);
    REQUIRE(next[0] == Catch::Approx(1.3));
    REQUIRE(next[1] == Catch::Approx(-2.6));
  }
}
