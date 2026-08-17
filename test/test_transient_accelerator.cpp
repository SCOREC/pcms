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
  // 0.5 is a representative initial under-relaxation factor. A production
  // input should be positive and chosen for the coupling being solved.
  tr::AitkenRelaxation accelerator(0.5);
  std::vector<Real> next;

  accelerator.BeginWindow();

  // x_in and x_out represent equally sized interface vectors before and after
  // one fixed-point evaluation. Their values are simple dummy DOFs chosen so
  // the residual and relaxed result can be calculated exactly.
  const Real residual = accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);

  // The residual is ||{2,-4}||_2 = sqrt(20).
  REQUIRE(residual == Catch::Approx(std::sqrt(20.0)));

  // With omega=0.5, x_next = x_in + omega*(x_out-x_in) = {1,-2}.
  REQUIRE(next == std::vector<Real>{1.0, -2.0});

  // Add history and then start a new window. The repeated first update must
  // again use omega=0.5 rather than a factor learned in the previous window.
  accelerator.Update(next, {1.5, -3.0}, next);
  accelerator.BeginWindow();
  accelerator.Update({0.0, 0.0}, {2.0, -4.0}, next);
  REQUIRE(next == std::vector<Real>{1.0, -2.0});
}
