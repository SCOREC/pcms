#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/application.hpp"
#include "pcms/utility/arrays.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

using pcms::Real;

namespace
{

Kokkos::View<Real*, pcms::DeviceMemorySpace> ToDevice(
  const std::vector<Real>& values)
{
  Kokkos::View<Real*, pcms::DeviceMemorySpace> device("test_vector",
                                                      values.size());
  pcms::CopyHostRank1ViewToDeviceView(device,
                                      pcms::make_const_array_view(values));
  return device;
}

std::vector<Real> ToHost(
  pcms::Rank1View<const Real, pcms::DeviceMemorySpace> view)
{
  std::vector<Real> values(view.extent(0));
  pcms::CopyRank1ViewToHost(pcms::make_array_view(values), view);
  return values;
}

// An application that adds 10 to every interface value on each AdvanceTo, so
// every step changes the state. It remembers the state before each step, up to
// a limit, and rolls back through that history.
class TestApplication final : public pcms::Application
{
public:
  explicit TestApplication(std::size_t max_steps = 8) : max_steps_(max_steps) {}

  std::string_view Name() const override { return "test-application"; }

  void AdvanceTo(Real target_time) override
  {
    Kokkos::View<Real*, pcms::DeviceMemorySpace> before("test_step",
                                                        interface_.extent(0));
    Kokkos::deep_copy(before, interface_);
    history_.push_back({time_, before});
    if (history_.size() > max_steps_)
      history_.erase(history_.begin());

    time_ = target_time;
    auto values = ToHost(pcms::MakeRank1View(interface_));
    for (Real& value : values)
      value += 10.0;
    pcms::CopyHostRank1ViewToDeviceView(interface_,
                                        pcms::make_const_array_view(values));
  }

  bool Restore(int index = -1) override
  {
    if (index >= 0 || static_cast<std::size_t>(-index) > history_.size())
      return false;
    // Go back to the state before the step, and forget the steps after it.
    const std::size_t target =
      history_.size() - static_cast<std::size_t>(-index);
    const Step step = history_[target];
    history_.resize(target);
    time_ = step.time;
    Kokkos::deep_copy(interface_, step.interface);
    return true;
  }

  pcms::Rank1View<const Real, pcms::DeviceMemorySpace> GetInterface(
    std::string_view) const override
  {
    return pcms::MakeRank1View(interface_);
  }

  void SetInterface(
    std::string_view,
    pcms::Rank1View<const Real, pcms::DeviceMemorySpace> values) override
  {
    pcms::CopyDeviceRank1ViewToDeviceView(interface_, values);
  }

  pcms::Capabilities GetCapabilities() const override
  {
    return {/*can_rollback=*/true, /*has_continuous_output=*/false,
            /*reports_qoi=*/false};
  }

  Real Time() const noexcept { return time_; }

private:
  struct Step
  {
    Real time;
    Kokkos::View<Real*, pcms::DeviceMemorySpace> interface;
  };

  std::size_t max_steps_;
  Real time_ = 0.0;
  Kokkos::View<Real*, pcms::DeviceMemorySpace> interface_{"test_interface", 2};
  std::vector<Step> history_;
};

// The interface values of an application, copied to the host.
std::vector<Real> Interface(const pcms::Application& application)
{
  return ToHost(application.GetInterface("interface"));
}

} // namespace

TEST_CASE("Application rolls back the last step by default", "[transient]")
{
  TestApplication application;
  // Use the base interface, as the coupling code will.
  pcms::Application& base = application;

  // These values represent interface DOFs in their coupling order.
  const auto initial_values = ToDevice({1.0, 2.0});
  base.SetInterface("interface", pcms::MakeRank1View(initial_values));
  base.AdvanceTo(0.5);
  base.AdvanceTo(1.0);
  REQUIRE(Interface(base)[0] == Catch::Approx(21.0));

  REQUIRE(base.Restore());

  const auto restored = Interface(base);
  REQUIRE(std::string(base.Name()) == "test-application");
  REQUIRE(application.Time() == Catch::Approx(0.5));
  REQUIRE(restored.size() == 2);
  REQUIRE(restored[0] == Catch::Approx(11.0));
  REQUIRE(restored[1] == Catch::Approx(12.0));
  REQUIRE(base.GetCapabilities().can_rollback);
  REQUIRE(base.ReportQoI().empty());
}

TEST_CASE("The index counts the steps to roll back", "[transient]")
{
  TestApplication application;
  const auto start = ToDevice({1.0, 2.0});
  application.SetInterface("interface", pcms::MakeRank1View(start));
  application.AdvanceTo(0.5);
  application.AdvanceTo(1.0);
  application.AdvanceTo(1.5);

  // -2 undoes the last two steps, which leaves the state after the first.
  REQUIRE(application.Restore(-2));
  REQUIRE(application.Time() == Catch::Approx(0.5));
  REQUIRE(Interface(application)[0] == Catch::Approx(11.0));

  // -1 undoes the one step that is left.
  REQUIRE(application.Restore(-1));
  REQUIRE(application.Time() == Catch::Approx(0.0));
  REQUIRE(Interface(application)[0] == Catch::Approx(1.0));
}

TEST_CASE("Repeating a window rolls back one step per sweep", "[transient]")
{
  // A Schwarz iteration advances every application once per sweep and rolls it
  // back before the next sweep.
  TestApplication application;
  const auto start = ToDevice({1.0, 2.0});
  application.SetInterface("interface", pcms::MakeRank1View(start));

  for (int sweep = 0; sweep < 3; ++sweep) {
    application.AdvanceTo(1.0);
    REQUIRE(Interface(application)[0] == Catch::Approx(11.0));

    REQUIRE(application.Restore());
    REQUIRE(application.Time() == Catch::Approx(0.0));
    const auto restored = Interface(application);
    REQUIRE(restored[0] == Catch::Approx(1.0));
    REQUIRE(restored[1] == Catch::Approx(2.0));
  }
}

TEST_CASE("An application that cannot honor a rollback returns false",
          "[transient]")
{
  // This application remembers only the last step.
  TestApplication application(1);
  const auto start = ToDevice({1.0, 2.0});
  application.SetInterface("interface", pcms::MakeRank1View(start));

  SECTION("nothing to roll back")
  {
    REQUIRE_FALSE(application.Restore());
    REQUIRE(Interface(application)[0] == Catch::Approx(1.0));
  }

  SECTION("more steps than it remembers")
  {
    application.AdvanceTo(0.5);
    application.AdvanceTo(1.0);

    REQUIRE_FALSE(application.Restore(-2));
    // The state did not change.
    REQUIRE(application.Time() == Catch::Approx(1.0));
    REQUIRE(Interface(application)[0] == Catch::Approx(21.0));

    // The step it does remember can still be rolled back.
    REQUIRE(application.Restore(-1));
    REQUIRE(application.Time() == Catch::Approx(0.5));
    REQUIRE(Interface(application)[0] == Catch::Approx(11.0));
  }

  SECTION("an index that is not negative")
  {
    application.AdvanceTo(0.5);
    REQUIRE_FALSE(application.Restore(0));
    REQUIRE_FALSE(application.Restore(1));
    REQUIRE(application.Time() == Catch::Approx(0.5));
  }
}
