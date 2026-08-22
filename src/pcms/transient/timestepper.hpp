#ifndef PCMS_TRANSIENT_TIMESTEPPER_HPP
#define PCMS_TRANSIENT_TIMESTEPPER_HPP

#include "pcms/utility/types.h"
#include <utility>

namespace pcms::transient
{

// Decides whether to accept a completed time window and selects the next step.
class Timestepper
{
public:
  // Return the positive step size used for the first coupling window.
  [[nodiscard]] virtual Real InitialStep() const = 0;

  // Given the current time step that was actually completed and its normalized error,
  // return whether to accept it and the positive step size to try next.
  virtual std::pair<bool, Real> Update(Real dt, Real err) = 0;

  virtual ~Timestepper() = default;
};

// Always accepts and keeps a constant time step.
class FixedTimestepper : public Timestepper
{
public:
  // Store the positive step size used for every window.
  explicit FixedTimestepper(Real dt);

  // Return the configured fixed step.
  Real InitialStep() const override;

  // Keep the configured fixed time step.
  std::pair<bool, Real> Update(Real dt, Real err) override;

private:
  Real dt_;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_TIMESTEPPER_HPP
