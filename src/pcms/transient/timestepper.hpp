#ifndef PCMS_TRANSIENT_TIMESTEPPER_HPP
#define PCMS_TRANSIENT_TIMESTEPPER_HPP

#include "pcms/utility/types.h"
#include <utility>

namespace pcms::transient
{

// Layer 2 — the window-Δt feedback controller (plan §3.3d). FEEDBACK, not
// predictive: the window is advanced at some Δt and measured, then Update
// decides accept/reject and the next Δt from the normalized coupling error.
//
// Step-size CONTRACT (the reason there is no Propose(t)): the controller learns
// Δt ONLY through Update's `dt` argument, which the orchestrator guarantees is
// the step ACTUALLY integrated (post-clamping to t_end / output points), paired
// with the error measured for exactly that step. The controller must build all
// of its history — previous errors AND previous step sizes — from those
// (dt, err) pairs, never from a value it proposed. InitialStep() only seeds the
// very first window and is read once. This makes the proposed-vs-actual Δt
// desync (which bites multi-step PI/PID controllers on a clamped final window)
// unrepresentable, and maps 1:1 onto SUNDIALS SUNAdaptController, whose
// primitive is EstimateStep(h, p, dsm)->hnew with no separate "propose".
//
// Accept/reject history: within one window Update may be called several times
// (reject → reject → accept) for the SAME physical window. Only the ACCEPTED
// (dt, err) pair advances the smooth multi-step history; the controller owns the
// accept decision, so it can manage this. Reset() clears history at a modelled
// discontinuity or across reuse (maps to SUNAdaptController_Reset).
//
// Collective control (plan §6): Update must consume the GLOBALLY-reduced error
// so every replica of the controller computes an identical dt_next and the
// participants stay in lockstep. Here (serial) the error is already global.
class Timestepper
{
public:
  // First-window Δt seed; read exactly once by the orchestrator.
  [[nodiscard]] virtual Real InitialStep() const = 0;
  // {accept?, dt_next}. `dt` is the actually-integrated step; `err` is
  // normalized to tol (accept iff ≤ 1).
  virtual std::pair<bool, Real> Update(Real dt, Real err) = 0;
  virtual ~Timestepper() = default;
};

// Ignores the error signal: always accepts, always returns the same Δt. The
// plan's baseline / debugging / CI stepper (§3.3d) and the driver for Step 1.
class FixedTimestepper : public Timestepper
{
public:
  explicit FixedTimestepper(Real dt) : dt_(dt) {}
  Real InitialStep() const override { return dt_; }
  std::pair<bool, Real> Update(Real, Real) override { return {true, dt_}; }

private:
  Real dt_;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_TIMESTEPPER_HPP

