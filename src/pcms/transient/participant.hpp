#ifndef PCMS_TRANSIENT_PARTICIPANT_HPP
#define PCMS_TRANSIENT_PARTICIPANT_HPP

#include "pcms/utility/types.h"
#include <any>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

// Transient simulation interfaces.
namespace pcms::transient
{
// Stores interface values retained by the Schwarz fixed-point iteration.
// Schwarz must preserve x_in while a participant advances or restores and may
// overwrite its FieldData. Aitken also needs stable current and previous
// residuals to compute the next iterate. A view into FieldData cannot provide
// that lifetime, so these iteration values are copied into owned storage.
class InterfaceState
{
public:
  // Create an empty state when no interface values are available yet.
  InterfaceState();

  // Create n interface DOFs initialized to value.
  explicit InterfaceState(std::size_t n, Real value = 0.0);

  // Take ownership of an existing vector of interface DOFs.
  explicit InterfaceState(std::vector<Real> dofs);

  // Return the number of interface DOFs.
  [[nodiscard]] std::size_t Size() const noexcept;

  // Access one DOF for modification.
  [[nodiscard]] Real& operator[](std::size_t i) noexcept;

  // Read one DOF without allowing modification.
  [[nodiscard]] Real operator[](std::size_t i) const noexcept;

  // Borrow all retained DOFs without transferring ownership.
  [[nodiscard]] std::span<const Real> View() const noexcept;

private:
  std::vector<Real> dofs_;
};

// Restart state saved at the beginning of a coupling window.
struct Checkpoint
{
  Real time = 0.0;
  std::any state; // Participant-defined payload.
};

// Optional features supported by a participant.
struct Capabilities
{
  bool can_restart = false;
  bool has_dense_output = false;
  bool reports_qoi = false;
};

// Interface implemented by each solver participating in the simulation.
class Participant
{
public:
  // Return the stable name used to identify this participant.
  [[nodiscard]] virtual std::string_view Name() const = 0;

  // Advance the solver from its current state to t_target.
  virtual void AdvanceTo(Real t_target) = 0;

  // Capture all state required to repeat the current coupling window.
  virtual Checkpoint Save() const = 0;

  // Restore a checkpoint previously created by this participant.
  virtual void Restore(const Checkpoint&) = 0;

  // Borrow values produced on the requested configured interface. The returned
  // view must remain valid until the next operation that changes participant
  // state.
  [[nodiscard]] virtual std::span<const Real> GetInterface(
    std::string_view name) const = 0;

  // Apply borrowed values with the size and ordering expected by the named
  // interface. Implementations must consume or copy them before returning.
  virtual void SetInterface(std::string_view name,
                            std::span<const Real> values) = 0;

  // Return optional scalar quantities used for error or conservation checks.
  // The returned storage must remain valid until the participant is advanced.
  [[nodiscard]] virtual std::span<const Real> ReportQoI() const;

  // Report the features this implementation can reliably provide.
  [[nodiscard]] virtual Capabilities GetCapabilities() const = 0;

  virtual ~Participant() = default;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_PARTICIPANT_HPP
