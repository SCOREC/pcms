#ifndef PCMS_TRANSIENT_TIMESTEPPER_HPP
#define PCMS_TRANSIENT_TIMESTEPPER_HPP

#include "pcms/utility/types.hpp"
#include <memory>

namespace pcms
{

/// Result of Timestepper::UpdateStep.
struct Step
{
  /// Whether the completed time window is accepted.
  bool accepted;
  /// Positive step size to try next.
  Real next_step;
};

/// Decides whether to accept a completed time window and selects the next step.
class Timestepper
{
public:
  /// \return The positive step size used for the first coupling window.
  [[nodiscard]] virtual Real InitialStep() const = 0;

  /// Decide whether to accept a completed time window and pick the next step.
  /// \param dt  The time step that was actually completed.
  /// \param err Normalized error of the completed step.
  /// \return Whether the step is accepted and the positive step size to try
  /// next.
  [[nodiscard]] virtual Step UpdateStep(Real dt, Real err) = 0;

  virtual ~Timestepper() = default;
};

/// A timestepper that always accepts and keeps the constant step `dt`.
/// @throws pcms_error unless `dt` is positive
[[nodiscard]] std::unique_ptr<Timestepper> CreateFixedTimestepper(Real dt);

} // namespace pcms

#endif // PCMS_TRANSIENT_TIMESTEPPER_HPP
