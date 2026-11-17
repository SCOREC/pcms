#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/application.hpp"
#include "pcms/utility/arrays.h"

#include <any>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using pcms::Real;

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

std::vector<Real> ToHost(
  pcms::Rank1View<const Real, pcms::DeviceMemorySpace> view)
{
  std::vector<Real> values(view.extent(0));
  pcms::CopyRank1ViewToHost(pcms::make_array_view(values), view);
  return values;
}

class TestApplication final : public pcms::Application
{
public:
  std::string_view Name() const override { return "test-application"; }

  void AdvanceTo(Real target_time) override { time_ = target_time; }

  pcms::Checkpoint Save() const override
  {
    // The checkpoint owns its own copy of the interface values.
    Vector saved("test_checkpoint", interface_.extent(0));
    Kokkos::deep_copy(saved, interface_);
    return {time_, saved};
  }

  void Restore(const pcms::Checkpoint& checkpoint) override
  {
    time_ = checkpoint.time;
    Kokkos::deep_copy(interface_,
                      std::any_cast<const Vector&>(checkpoint.state));
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
    return {/*can_restart=*/true, /*has_dense_output=*/false,
            /*reports_qoi=*/false};
  }

  Real Time() const noexcept { return time_; }

private:
  Real time_ = 0.0;
  Vector interface_{"test_interface", 2};
};

} // namespace

TEST_CASE("Application state can be saved and restored", "[transient]")
{
  TestApplication application;
  // Use the base interface, as the coupling code will.
  pcms::Application& base = application;

  // These values represent interface DOFs in their coupling order.
  const Vector initial_values = ToDevice({1.0, 2.0});
  base.SetInterface("interface", pcms::MakeRank1View(initial_values));
  base.AdvanceTo(0.5);
  const pcms::Checkpoint checkpoint = base.Save();

  // Change both parts of the state so restore cannot pass accidentally.
  const Vector replacement_values = ToDevice({3.0, 4.0});
  base.SetInterface("interface", pcms::MakeRank1View(replacement_values));
  base.AdvanceTo(1.0);
  base.Restore(checkpoint);

  const auto restored = ToHost(base.GetInterface("interface"));

  REQUIRE(std::string(base.Name()) == "test-application");
  REQUIRE(application.Time() == Catch::Approx(0.5));
  REQUIRE(restored.size() == 2);
  REQUIRE(restored[0] == Catch::Approx(1.0));
  REQUIRE(restored[1] == Catch::Approx(2.0));
  REQUIRE(base.GetCapabilities().can_restart);
  REQUIRE(base.ReportQoI().empty());
}

TEST_CASE("A checkpoint can be restored more than once", "[transient]")
{
  // A Schwarz iteration restores the same window-start checkpoint every sweep.
  TestApplication application;
  const Vector start = ToDevice({1.0, 2.0});
  application.SetInterface("interface", pcms::MakeRank1View(start));
  const pcms::Checkpoint checkpoint = application.Save();

  for (int sweep = 0; sweep < 2; ++sweep) {
    const Vector changed = ToDevice({5.0 + sweep, 6.0 + sweep});
    application.SetInterface("interface", pcms::MakeRank1View(changed));
    application.AdvanceTo(1.0 + sweep);

    application.Restore(checkpoint);
    const auto restored = ToHost(application.GetInterface("interface"));
    REQUIRE(application.Time() == Catch::Approx(0.0));
    REQUIRE(restored[0] == Catch::Approx(1.0));
    REQUIRE(restored[1] == Catch::Approx(2.0));
  }
}
