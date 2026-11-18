#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/accelerator.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"

#include <cmath>
#include <cstddef>
#include <vector>

using pcms::Real;
using Status = pcms::InterfaceAccelerator::UpdateStatus;

namespace
{

using Vector = Kokkos::View<Real*, pcms::DeviceMemorySpace>;

Vector ToDevice(const std::vector<Real>& values)
{
  Vector device("test_vector", values.size());
  pcms::CopyHostRank1ViewToDeviceView(device,
                                      pcms::make_const_array_view(values));
  return device;
}

std::vector<Real> ToHost(const Vector& device)
{
  std::vector<Real> values(device.extent(0));
  pcms::CopyRank1ViewToHost(pcms::make_array_view(values),
                            pcms::MakeRank1View(device));
  return values;
}

// Runs one Update() on device data through array views.
Real Run(pcms::InterfaceAccelerator& accelerator, const Vector& in,
         const Vector& out, Vector& next)
{
  return accelerator.Update(pcms::MakeRank1View(in), pcms::MakeRank1View(out),
                            pcms::MakeRank1View(next));
}

} // namespace

TEST_CASE("Aitken relaxation updates and resets the interface iterate",
          "[transient]")
{
  // Use an initial relaxation factor.
  pcms::AitkenRelaxation accelerator(0.5);
  Vector next("next", 2);

  accelerator.BeginWindow();
  REQUIRE(accelerator.LastStatus() == Status::Initial);

  const Real first_residual =
    Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({2.0, -4.0}), next);

  // The returned residual norm is ||{2,-4}||_2, and omega=0.5.
  auto result = ToHost(next);
  REQUIRE(first_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(result[0] == Catch::Approx(1.0));
  REQUIRE(result[1] == Catch::Approx(-2.0));
  REQUIRE(accelerator.LastStatus() == Status::Initial);

  // A second residual lets Aitken compute omega=2/3.
  const Real second_residual =
    Run(accelerator, next, ToDevice({1.5, -3.0}), next);
  result = ToHost(next);
  REQUIRE(second_residual == Catch::Approx(std::sqrt(1.25)));
  REQUIRE(result[0] == Catch::Approx(4.0 / 3.0));
  REQUIRE(result[1] == Catch::Approx(-8.0 / 3.0));
  REQUIRE(accelerator.LastStatus() == Status::Updated);

  // The status is available through the base interface.
  const pcms::InterfaceAccelerator& base = accelerator;
  REQUIRE(base.LastStatus() == Status::Updated);

  // A new window must discard the history and reuse omega=0.5.
  accelerator.BeginWindow();
  REQUIRE(accelerator.LastStatus() == Status::Initial);
  const Real third_residual =
    Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({2.0, -4.0}), next);
  result = ToHost(next);
  REQUIRE(third_residual == Catch::Approx(std::sqrt(20.0)));
  REQUIRE(result[0] == Catch::Approx(1.0));
  REQUIRE(result[1] == Catch::Approx(-2.0));
}

TEST_CASE("Aitken relaxation reports what happened to the relaxation factor",
          "[transient]")
{
  Vector next("next", 2);

  SECTION("the factor is kept when the residual does not change")
  {
    pcms::AitkenRelaxation accelerator(0.5);
    accelerator.BeginWindow();

    // Both iterations have the residual {2,-4}, so the change in residual is 0.
    (void)Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({2.0, -4.0}), next);
    (void)Run(accelerator, ToDevice({1.0, -2.0}), ToDevice({3.0, -6.0}), next);

    // omega stays 0.5, so next = {1,-2} + 0.5 * {2,-4}.
    const auto result = ToHost(next);
    REQUIRE(accelerator.LastStatus() == Status::Fallback);
    REQUIRE(result[0] == Catch::Approx(2.0));
    REQUIRE(result[1] == Catch::Approx(-4.0));
  }

  SECTION("the factor is reset to the initial factor when it would be zero")
  {
    pcms::AitkenRelaxation accelerator(0.5);
    accelerator.BeginWindow();

    // The residuals {1,0} and {1,1} give a change {0,1} that is orthogonal to
    // the previous residual, so the new omega is 0 and is reset to 0.5.
    (void)Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({1.0, 0.0}), next);
    (void)Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({1.0, 1.0}), next);

    const auto result = ToHost(next);
    REQUIRE(accelerator.LastStatus() == Status::Fallback);
    REQUIRE(result[0] == Catch::Approx(0.5));
    REQUIRE(result[1] == Catch::Approx(0.5));
  }

  SECTION("the factor is clamped when it exceeds the limit")
  {
    // The unclamped omega is 2/3, which is above the limit of 0.6.
    pcms::AitkenRelaxation accelerator(0.5, 0.6);
    accelerator.BeginWindow();

    (void)Run(accelerator, ToDevice({0.0, 0.0}), ToDevice({2.0, -4.0}), next);
    (void)Run(accelerator, ToDevice({1.0, -2.0}), ToDevice({1.5, -3.0}), next);

    // next = {1,-2} + 0.6 * {0.5,-1}.
    const auto result = ToHost(next);
    REQUIRE(accelerator.LastStatus() == Status::Limited);
    REQUIRE(result[0] == Catch::Approx(1.3));
    REQUIRE(result[1] == Catch::Approx(-2.6));
  }
}

TEST_CASE("Aitken relaxation rejects views of different sizes", "[transient]")
{
  pcms::AitkenRelaxation accelerator(0.5);
  accelerator.BeginWindow();

  const Vector in = ToDevice({0.0, 0.0});
  const Vector longer = ToDevice({0.0, 0.0, 0.0});
  Vector next("next", 2);
  Vector longer_next("longer_next", 3);

  SECTION("x_out has a different size than x_in")
  {
    REQUIRE_THROWS_AS(Run(accelerator, in, longer, next), pcms::pcms_error);
  }

  SECTION("x_next has a different size than x_in")
  {
    REQUIRE_THROWS_AS(Run(accelerator, in, in, longer_next), pcms::pcms_error);
  }

  SECTION("the size changes within a window")
  {
    (void)Run(accelerator, in, ToDevice({1.0, 0.0}), next);
    REQUIRE_THROWS_AS(Run(accelerator, longer, longer, longer_next),
                      pcms::pcms_error);

    // A new window may use a different size.
    accelerator.BeginWindow();
    REQUIRE_NOTHROW(Run(accelerator, longer, longer, longer_next));
  }
}
