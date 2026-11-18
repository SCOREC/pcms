#ifndef PCMS_TRANSIENT_ACCELERATOR_HPP
#define PCMS_TRANSIENT_ACCELERATOR_HPP

#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

namespace pcms
{

/**
 * @brief Interface for accelerating the fixed-point iteration that couples
 * two applications across an interface.
 *
 * Coupling interface values are updated by an accelerator during a fixed-point
 * iteration. Each coupling window is started with BeginWindow() and then
 * advanced one iteration at a time with Update(). How the last Update()
 * produced the next state is reported by LastStatus().
 *
 * All data is passed as views in device memory and is used in place.
 */
class InterfaceAccelerator
{
public:
  /**
   * @brief Describes how the last Update() produced the next interface state.
   *
   * Exactly one value is reported.
   */
  enum class UpdateStatus
  {
    /// The accelerator used its initial parameters.
    Initial,
    /// The accelerator applied its normal update.
    Updated,
    /// The accelerator kept or restored its fallback parameters.
    Fallback,
    /// The normal update was applied but limited to a permitted range.
    Limited
  };

  /**
   * @brief Reset the iteration history before solving a new coupling window.
   */
  virtual void BeginWindow() = 0;

  /**
   * @brief Compute the next interface state and return the residual norm.
   *
   * The residual of the iteration is @c x_out - @c x_in.
   *
   * @param x_in Interface state that was passed to the application.
   * @param x_out Interface state produced by the application. Must have the
   * same size as @c x_in.
   * @param x_next Receives the next interface state. Must have the same size as
   * @c x_in. It may alias @c x_in.
   * @return The L2 norm of the residual.
   * @throws pcms_error If the sizes differ, or if the size changed since the
   * first Update() of the window.
   */
  [[nodiscard]] virtual Real Update(
    Rank1View<const Real, DeviceMemorySpace> x_in,
    Rank1View<const Real, DeviceMemorySpace> x_out,
    Rank1View<Real, DeviceMemorySpace> x_next) = 0;

  /**
   * @brief Report how the last Update() produced the next interface state.
   *
   * @return The status of the last Update(), or UpdateStatus::Initial if there
   * has been no Update() since construction or BeginWindow().
   */
  [[nodiscard]] virtual UpdateStatus LastStatus() const noexcept = 0;

  virtual ~InterfaceAccelerator() = default;
};

/**
 * @brief Aitken relaxation of the interface state.
 *
 * Computes the next state as @c x_in + omega * (x_out - x_in). The relaxation
 * factor omega is the initial factor on the first iteration of a window. On
 * later iterations it is updated from the residuals of two consecutive
 * iterations, r_k and r_{k-1}:
 *
 * @code
 * d = r_k - r_{k-1}
 * omega_{k+1} = -omega_k * dot(r_{k-1}, d) / dot(d, d)
 * @endcode
 *
 * If the new factor is not finite or is exactly zero, it is reset to the
 * initial factor. The factor is then clamped to the range [-omega_max,
 * omega_max].
 *
 * The residual of the previous iteration is kept in device memory.
 */
class AitkenRelaxation : public InterfaceAccelerator
{
public:
  /**
   * @brief Construct an Aitken relaxation accelerator.
   *
   * @param omega0 The initial relaxation factor, also used when the factor is
   * reset.
   * @param omega_max The limit on the magnitude of the relaxation factor.
   */
  explicit AitkenRelaxation(Real omega0 = 0.5, Real omega_max = 1e6);

  /**
   * @brief Reset the relaxation factor to the initial factor and discard the
   * stored residual.
   */
  void BeginWindow() override;

  /**
   * @copydoc InterfaceAccelerator::Update
   */
  [[nodiscard]] Real Update(Rank1View<const Real, DeviceMemorySpace> x_in,
                            Rank1View<const Real, DeviceMemorySpace> x_out,
                            Rank1View<Real, DeviceMemorySpace> x_next) override;

  /**
   * @copydoc InterfaceAccelerator::LastStatus
   */
  [[nodiscard]] UpdateStatus LastStatus() const noexcept override
  {
    return last_status_;
  }

private:
  Real omega0_;    ///< Initial relaxation factor.
  Real omega_;     ///< Current relaxation factor.
  Real omega_max_; ///< Limit on the magnitude of omega_.

  /// Residual of the previous iteration, in device memory.
  Kokkos::View<Real*, DeviceMemorySpace> previous_residual_;
  bool have_previous_residual_ = false; ///< Whether previous_residual_ is set.

  /// Outcome of the last Update().
  UpdateStatus last_status_ = UpdateStatus::Initial;
};

} // namespace pcms

#endif // PCMS_TRANSIENT_ACCELERATOR_HPP
