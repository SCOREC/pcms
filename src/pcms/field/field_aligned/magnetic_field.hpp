#ifndef PCMS_FIELD_FIELD_ALIGNED_MAGNETIC_FIELD_HPP
#define PCMS_FIELD_FIELD_ALIGNED_MAGNETIC_FIELD_HPP

#include "pcms/utility/eqdsk.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

/**
 * @brief Axisymmetric magnetic field from EQDSK equilibrium data
 *
 * Coordinates are right-handed cylindrical (R, phi, Z), with R and Z in meters
 * and phi in radians. With psi the poloidal flux in Wb/rad and I(psi) the
 * poloidal current function in T*m, the field is
 *   B = s grad(phi) x grad(psi) + I(psi) grad(phi)
 * and a field line, with phi as the independent variable, satisfies
 *   dR/dphi =  s * R * psi_Z / I(psi)
 *   dZ/dphi = -s * R * psi_R / I(psi)
 *
 * The poloidal sign s selects the convention of the equilibrium:
 * s = +1 for B = grad(phi) x grad(psi) + I grad(phi), and s = -1 for
 * B = grad(psi) x grad(phi) + I grad(phi) (STOMMS, EFIT). It has no default
 * because codes differ and a wrong sign reverses the poloidal motion.
 *
 * I keeps the sign stored in the EQDSK data; its magnitude is never taken,
 * because the sign sets the toroidal direction of B.
 *
 * psi(R, Z) is a not-a-knot bicubic spline on the EQDSK grid and I(psi) is a
 * not-a-knot cubic spline on psi_grid. A strictly descending psi_grid is
 * reversed together with I_psi first. I is evaluated at psi clamped to the
 * psi_grid range, so it is constant outside that range (vacuum field).
 */
class AxisymmetricMagneticField
{
public:
  /**
   * @brief Builds the psi(R, Z) and I(psi) splines
   *
   * @param eqdsk equilibrium data; it is copied, so it need not outlive this
   *        object
   * @param poloidal_sign s in the field-line equations, +1 or -1
   *
   * @throws pcms_error if poloidal_sign is not +1 or -1, the R-Z grid has
   *         fewer than 2 points or a non-positive spacing in either direction,
   *         or psi_grid is not uniform and strictly ascending or descending
   */
  AxisymmetricMagneticField(const EQDSKData& eqdsk, int poloidal_sign);

  /**
   * @brief Evaluates psi, its gradient and I at a batch of points in one call
   *
   * A point with R outside [GetRMin(), GetRMax()] or Z outside
   * [GetZMin(), GetZMax()] gets inside = false and is not extrapolated; its
   * psi, psi_R, psi_Z and I are set to zero.
   *
   * @param R major radius of each point, in m (size n)
   * @param Z height of each point, in m (size n)
   * @param psi receives psi in Wb/rad (size n)
   * @param psi_R receives d(psi)/dR (size n)
   * @param psi_Z receives d(psi)/dZ (size n)
   * @param I receives I(psi) in T*m, with its stored sign (size n)
   * @param inside receives whether the point lies in the R-Z grid box (size n)
   *
   * @throws pcms_error if the views do not all have the same size
   */
  void Evaluate(Kokkos::View<const Real*, DeviceMemorySpace> R,
                Kokkos::View<const Real*, DeviceMemorySpace> Z,
                Kokkos::View<Real*, DeviceMemorySpace> psi,
                Kokkos::View<Real*, DeviceMemorySpace> psi_R,
                Kokkos::View<Real*, DeviceMemorySpace> psi_Z,
                Kokkos::View<Real*, DeviceMemorySpace> I,
                Kokkos::View<bool*, DeviceMemorySpace> inside) const;

  /// Smallest R of the EQDSK grid, in m
  [[nodiscard]] Real GetRMin() const noexcept { return r_min_; }
  /// Largest R of the EQDSK grid, in m
  [[nodiscard]] Real GetRMax() const noexcept { return r_max_; }
  /// Smallest Z of the EQDSK grid, in m
  [[nodiscard]] Real GetZMin() const noexcept { return z_min_; }
  /// Largest Z of the EQDSK grid, in m
  [[nodiscard]] Real GetZMax() const noexcept { return z_max_; }
  /// Smallest psi_grid value; I(psi) below it equals I at this value
  [[nodiscard]] Real GetPsiMin() const noexcept { return psi_min_; }
  /// Largest psi_grid value; I(psi) above it equals I at this value
  [[nodiscard]] Real GetPsiMax() const noexcept { return psi_max_; }
  /// The poloidal sign s, +1 or -1
  [[nodiscard]] int GetPoloidalSign() const noexcept { return poloidal_sign_; }

private:
  // Defined in the .cpp so this header does not pull in the spline library.
  struct Impl;
  std::shared_ptr<Impl> impl_;
  Real r_min_;
  Real r_max_;
  Real z_min_;
  Real z_max_;
  Real psi_min_;
  Real psi_max_;
  int poloidal_sign_;
};

} // namespace pcms

#endif // PCMS_FIELD_FIELD_ALIGNED_MAGNETIC_FIELD_HPP
