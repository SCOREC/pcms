#ifndef PCMS_TRANSIENT_APPLICATION_HPP
#define PCMS_TRANSIENT_APPLICATION_HPP

#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <span>
#include <string_view>

namespace pcms
{

/**
 * @brief Describes the optional transient features an application provides.
 */
struct Capabilities
{
  /// The application can roll back to an earlier step.
  bool can_rollback = false;
  /// The application can evaluate its solution between time steps.
  bool has_continuous_output = false;
  /// The application reports quantities for error and conservation checks.
  bool reports_qoi = false;
};

/**
 * @brief Abstract base class that an application adapter implements to take
 * part in a transient coupled simulation.
 */
class Application
{
public:
  /**
   * @brief Return the stable name used in the coupling configuration.
   */
  [[nodiscard]] virtual std::string_view Name() const = 0;

  /**
   * @brief Advance from the current state to the requested time.
   *
   * @param target_time time to advance to.
   */
  virtual void AdvanceTo(Real target_time) = 0;

  /**
   * @brief Roll the application back by a number of steps.
   *
   * A step is one call to AdvanceTo(). The index counts back from the current
   * state: -1 undoes the last AdvanceTo(), -2 undoes the last two, and so on.
   *
   * @param index Negative number of steps to undo. The default undoes the last
   * step.
   * @return True if the application rolled back by the requested number of
   * steps. False if the adapter cannot roll back, in which case its state is
   * unchanged.
   */
  [[nodiscard]] virtual bool Restore(int index = -1) = 0;

  /**
   * @brief Borrow the values produced on a named coupling interface.
   *
   * The view stays valid until the next call to AdvanceTo(), SetInterface() or
   * Restore().
   *
   * @param name Name of the coupling interface.
   * @return A view of the interface values in device memory.
   */
  [[nodiscard]] virtual Rank1View<const Real, DeviceMemorySpace> GetInterface(
    std::string_view name) const = 0;

  /**
   * @brief Apply values to a named coupling interface.
   *
   * Implementations must use or copy them before returning.
   *
   * @param name Name of the coupling interface.
   * @param values Interface values in device memory.
   */
  virtual void SetInterface(
    std::string_view name, Rank1View<const Real, DeviceMemorySpace> values) = 0;

  /**
   * @brief Return quantities used for error and conservation checks.
   *
   * @return The reported quantities, or an empty span by default.
   */
  [[nodiscard]] virtual std::span<const Real> ReportQoI() const { return {}; }

  /**
   * @brief Report which  transient features this application provides.
   */
  [[nodiscard]] virtual Capabilities GetCapabilities() const = 0;

  virtual ~Application() = default;
};

} // namespace pcms

#endif // PCMS_TRANSIENT_APPLICATION_HPP
