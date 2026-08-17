#ifndef PCMS_TRANSIENT_PARTICIPANT_HPP
#define PCMS_TRANSIENT_PARTICIPANT_HPP

#include "pcms/utility/types.h"
#include <any>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

// Transient coupled-simulation layer (docs/transient_coupled_simulation_plan.md).
//
// The core contracts remain transport- and mesh-independent. Redev, ADIOS2,
// Omega_h, and field-transfer implementations live in participant_adapter/.
namespace pcms::transient
{

// InterfaceState carries the interface transmission DOFs (the overlap Schwarz
// unknown, plan §7) for a single named interface. For the in-process demo this
// is a host vector; §10.2 records the intended evolution: wrap these DOFs as a
// SUNDIALS N_Vector over the (possibly distributed) interface layout so KINSOL
// and the collective reductions (§6) operate in parallel / on device without a
// host copy per sweep.
class InterfaceState
{
public:
  InterfaceState() = default;
  explicit InterfaceState(std::size_t n, Real value = 0.0) : dofs_(n, value) {}
  explicit InterfaceState(std::vector<Real> dofs) : dofs_(std::move(dofs)) {}

  [[nodiscard]] std::size_t Size() const noexcept { return dofs_.size(); }
  [[nodiscard]] Real& operator[](std::size_t i) noexcept { return dofs_[i]; }
  [[nodiscard]] Real operator[](std::size_t i) const noexcept
  {
    return dofs_[i];
  }
  [[nodiscard]] std::span<const Real> View() const noexcept { return dofs_; }
  [[nodiscard]] std::vector<Real>& Data() noexcept { return dofs_; }
  [[nodiscard]] const std::vector<Real>& Data() const noexcept { return dofs_; }

private:
  std::vector<Real> dofs_;
};

// A code becomes couplable by producing a Checkpoint of its start-of-window
// state and restoring from it (plan §3.1). The payload is type-erased: each
// adapter stores whatever it needs (for the demo solver, its nodal field plus
// boundary values). Save/Restore re-entrancy is the one hard requirement
// Schwarz imposes — a window is re-integrated once per Schwarz iteration.
struct Checkpoint
{
  Real time = 0.0;
  std::any state; // adapter-defined payload
};

// What information a participant exposes, so the coupler advertises only the
// guarantee the assembled configuration can actually back (plan §7 honesty
// discipline; consumed by TransientOverlapCoupled::NegotiateCapabilities).
struct Capabilities
{
  bool can_restart = false;      // ⇒ Schwarz + step-doubling are possible
  bool has_dense_output = false; // ⇒ cheap overlap-residual coupling-error est.
  bool reports_qoi = false;      // ⇒ global-conservation monitoring / step-doubling QoI
};

// Layer 0 — the only thing a code implements. A thin, non-invasive adapter may
// wrap a local solver directly or be exposed to an MPMD driver through the
// optional transient participant-adapter component. The transient layer never
// invents a new application field type; it drives the interface transmission
// data the adapter already owns.
class Participant
{
public:
  [[nodiscard]] virtual std::string_view Name() const = 0;

  // Advance the code from its current time to t_target. The code subcycles
  // internally with whatever integrator it owns (black box) using the interface
  // BCs currently set on it. Repeated calls after Restore must be reproducible.
  virtual void AdvanceTo(Real t_target) = 0;

  virtual Checkpoint Save() const = 0;
  virtual void Restore(const Checkpoint&) = 0;

  // Produce this participant's transmission data for a named interface (what a
  // neighbour consumes). For a mesh code this is a field sampled/restricted onto
  // the interface; here it is the subdomain solution evaluated at the overlap.
  [[nodiscard]] virtual InterfaceState GetInterface(
    std::string_view name) const = 0;

  // Impose a received interface iterate as this participant's transmission BC
  // for the next AdvanceTo.
  virtual void SetInterface(std::string_view name, const InterfaceState&) = 0;

  // Optional interior scalars (plan Rung 1): used by the step-doubling estimator
  // as the compared quantity of interest and, eventually, to monitor global
  // conservation. Empty ⇒ estimators fall back to the interface state.
  [[nodiscard]] virtual std::span<const Real> ReportQoI() const { return {}; }

  [[nodiscard]] virtual Capabilities GetCapabilities() const = 0;

  virtual ~Participant() = default;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_PARTICIPANT_HPP

