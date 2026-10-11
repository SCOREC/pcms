#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/accelerator.hpp"

#include <cmath>
#include <vector>

namespace tr = pcms::transient;
using pcms::Real;

TEST_CASE("Aitken relaxation updates and resets the interface iterate",
          "[transient]")
{
  // Use a representative initial under-relaxation factor.
  tr::AitkenRelaxation accelerator(0.5);
  std::vector<Real> next;

  accelerator.BeginWindow();

  // These vectors represent the interface iterate before and after one
  // fixed-point evaluation.
  const Real first_residual = accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);

  // The returned residual is ||{2,-4}||_2, and the first update uses omega=0.5.
  REQUIRE(first_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(next[0] == Catch::Approx(1.0));
  REQUIRE(next[1] == Catch::Approx(-2.0));

  // A second residual lets Aitken compute omega=2/3 from iteration history.
  const Real second_residual = accelerator.Update(next, {1.5, -3.0}, next);
  REQUIRE(second_residual == Catch::Approx(std::sqrt(1.25)));
  REQUIRE(next[0] == Catch::Approx(4.0 / 3.0));
  REQUIRE(next[1] == Catch::Approx(-8.0 / 3.0));

  // A new window must discard that history and reuse omega=0.5.
  accelerator.BeginWindow();
  const Real third_residual = accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);
  REQUIRE(third_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(next[0] == Catch::Approx(1.0));
  REQUIRE(next[1] == Catch::Approx(-2.0));
}
