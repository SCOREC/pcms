#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/timestepper.hpp"

namespace tr = pcms::transient;

TEST_CASE("FixedTimestepper always accepts its configured step",
          "[transient]")
{
  // 0.25 represents the positive fixed window size selected by the user.
  tr::FixedTimestepper timestepper(0.25);

  // The completed step (0.1) and large normalized error (100) are deliberately
  // different from the configured step. A fixed controller must ignore both.
  const auto [accepted, next_step] = timestepper.Update(0.1, 100.0);

  // The first window must use the user-configured size.
  REQUIRE(timestepper.InitialStep() == Catch::Approx(0.25));

  // A fixed timestepper never rejects a completed window based on error.
  REQUIRE(accepted);

  // Subsequent windows must retain the same configured size.
  REQUIRE(next_step == Catch::Approx(0.25));
}
