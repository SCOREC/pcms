#include "pcms/transient/timestepper.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms
{

namespace
{

/// Always accepts and keeps a constant time step.
class FixedTimestepper final : public Timestepper
{
public:
  explicit FixedTimestepper(Real dt) : dt_(dt) {}

  [[nodiscard]] Real InitialStep() const override { return dt_; }

  /// Ignores \p dt and \p err; always accepts and returns the fixed step.
  [[nodiscard]] Step UpdateStep(Real /*dt*/, Real /*err*/) override
  {
    return {true, dt_};
  }

private:
  Real dt_;
};

} // namespace

std::unique_ptr<Timestepper> CreateFixedTimestepper(Real dt)
{
  if (!(dt > 0)) {
    throw pcms_error("CreateFixedTimestepper: the step must be positive");
  }
  return std::make_unique<FixedTimestepper>(dt);
}

} // namespace pcms
