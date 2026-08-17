#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/transient/participant.hpp"

#include <any>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace tr = pcms::transient;
using pcms::Real;

namespace
{

class TestParticipant final : public tr::Participant
{
public:
  std::string_view Name() const override { return "test-participant"; }

  void AdvanceTo(Real target_time) override { time_ = target_time; }

  tr::Checkpoint Save() const override
  {
    const auto values = interface_.View();
    return {time_, std::vector<Real>(values.begin(), values.end())};
  }

  void Restore(const tr::Checkpoint& checkpoint) override
  {
    time_ = checkpoint.time;
    interface_ =
      tr::InterfaceState(std::any_cast<std::vector<Real>>(checkpoint.state));
  }

  std::span<const Real> GetInterface(std::string_view) const override
  {
    return interface_.View();
  }

  void SetInterface(std::string_view, std::span<const Real> values) override
  {
    interface_ =
      tr::InterfaceState(std::vector<Real>(values.begin(), values.end()));
  }

  tr::Capabilities GetCapabilities() const override
  {
    return {/*can_restart=*/true, /*has_dense_output=*/false,
            /*reports_qoi=*/false};
  }

  Real Time() const noexcept { return time_; }

private:
  Real time_ = 0.0;
  tr::InterfaceState interface_{2, 0.0};
};

} // namespace

TEST_CASE("Participant state can be saved and restored", "[transient]")
{
  TestParticipant participant;

  // The vector represents FieldData-owned overlap DOFs in the ordering expected
  // by the participant. The numbers are arbitrary but distinct so that the test
  // can detect an ordering or restore error. SetInterface receives only a view
  // and must consume or copy it during the call.
  const std::vector<Real> initial_values{1.0, 2.0};
  participant.SetInterface("interface", std::span<const Real>(initial_values));

  // 0.5 represents the start of a coupling window. A real solver would advance
  // its complete numerical state to this time before creating the checkpoint.
  participant.AdvanceTo(0.5);
  const auto checkpoint = participant.Save();

  // Change both time and interface values so Restore cannot pass by leaving the
  // current state untouched.
  const std::vector<Real> replacement_values{3.0, 4.0};
  participant.SetInterface("interface",
                           std::span<const Real>(replacement_values));
  participant.AdvanceTo(1.0);
  participant.Restore(checkpoint);

  const auto restored = participant.GetInterface("interface");

  // The participant name must remain stable for coupling configuration.
  REQUIRE(std::string(participant.Name()) == "test-participant");

  // Restore must return to the saved window time, not the later time 1.0.
  REQUIRE(participant.Time() == Catch::Approx(0.5));

  // GetInterface returns a borrowed view. It must expose the restored shape and
  // original DOFs rather than the replacement values {3.0, 4.0}.
  REQUIRE(restored.size() == 2);
  REQUIRE(restored[0] == Catch::Approx(1.0));
  REQUIRE(restored[1] == Catch::Approx(2.0));

  // The transient driver may repeat a window only when restart is supported.
  REQUIRE(participant.GetCapabilities().can_restart);
}
