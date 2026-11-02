#ifndef PCMS_FIELD_FIELD_ALIGNED_FIELD_LINE_TRACER_HPP
#define PCMS_FIELD_FIELD_ALIGNED_FIELD_LINE_TRACER_HPP

#include "pcms/field/field_aligned/magnetic_field.hpp"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <vector>

namespace pcms
{

/**
 * @brief Positions of a batch of field lines at a list of toroidal angles
 */
struct TracedPoints
{
  /// Toroidal angle of each station, in rad (size S); stations(0) is the start
  Kokkos::View<Real*, DeviceMemorySpace> stations;
  /// positions(k, i, 0) = R and positions(k, i, 1) = Z, in m, of field line i
  /// at station k (size S x n x 2)
  Kokkos::View<Real***, Kokkos::LayoutRight, DeviceMemorySpace> positions;
  /// false once field line i has left the R-Z grid box; its later positions
  /// repeat its last valid position (size n)
  Kokkos::View<bool*, DeviceMemorySpace> valid;
};

/**
 * @brief Traces a batch of field lines of an axisymmetric field in phi
 *
 * Integrates, with phi as the independent variable,
 *   dR/dphi =  s * R * psi_Z / I(psi)
 *   dZ/dphi = -s * R * psi_R / I(psi)
 * with s = field.GetPoloidalSign() and I keeping its stored sign, using
 * classic fourth-order Runge-Kutta. R is re-evaluated at every stage; holding
 * it fixed within a step would be only first-order accurate. Between
 * consecutive stations it takes ceil(|phi_{k+1} - phi_k| / max_dphi) equal
 * substeps, so every station is reached exactly, and each stage evaluates the
 * field for all points in one batched call.
 *
 * A field line becomes invalid when a stage or the end of a step leaves the
 * R-Z grid box, or I(psi) is zero there; it keeps the position from the start
 * of that step. A start outside the box is invalid from the start. Invalid
 * points never throw.
 *
 * With project_to_psi, every valid point is moved back to the psi of its
 * starting point after each station by three Newton iterations
 *   x -= (psi - psi0) * grad(psi) / |grad(psi)|^2.
 * A zero gradient leaves the point unchanged; a projection that leaves the box
 * marks the point invalid and keeps its unprojected position.
 *
 * @param field the axisymmetric magnetic field
 * @param R0 starting R of each field line at phi = stations[0], in m (size n)
 * @param Z0 starting Z of each field line at phi = stations[0], in m (size n)
 * @param stations toroidal angles to record, in rad, strictly increasing or
 *        strictly decreasing; stations[0] is the starting phi
 * @param max_dphi largest substep in |phi|, in rad, > 0
 * @param project_to_psi project each point back to its starting psi at every
 *        station
 * @return the stations, the positions at every station and the validity of
 *         each field line
 *
 * @throws pcms_error if R0 and Z0 differ in size, stations is empty or not
 *         strictly monotonic, or max_dphi is not positive
 */
TracedPoints TraceFieldLines(const AxisymmetricMagneticField& field,
                             Kokkos::View<const Real*, DeviceMemorySpace> R0,
                             Kokkos::View<const Real*, DeviceMemorySpace> Z0,
                             const std::vector<Real>& stations, Real max_dphi,
                             bool project_to_psi = false);

/**
 * @brief TraceFieldLines with the stations in a device view
 *
 * The stations are copied to the host once, since the step sizes are chosen
 * there; see the std::vector overload for the rest.
 */
TracedPoints TraceFieldLines(
  const AxisymmetricMagneticField& field,
  Kokkos::View<const Real*, DeviceMemorySpace> R0,
  Kokkos::View<const Real*, DeviceMemorySpace> Z0,
  Kokkos::View<const Real*, DeviceMemorySpace> stations, Real max_dphi,
  bool project_to_psi = false);

/**
 * @brief Toroidal stations of wedge layers with cubic geometry along the field
 *
 * Layer k = 0, ..., num_layers - 1 contributes
 * phi0 + k * dphi + {0, dphi / 3, 2 * dphi / 3}, and the list ends with
 * phi0 + num_layers * dphi: 3 * num_layers + 1 values in total.
 *
 * @param phi0 toroidal angle of the first plane, in rad
 * @param dphi toroidal width of one layer, in rad
 * @param num_layers number of wedge layers, >= 0
 * @return the stations in order
 *
 * @throws pcms_error if num_layers is negative
 */
std::vector<Real> WedgeStations(Real phi0, Real dphi, int num_layers);

} // namespace pcms

#endif // PCMS_FIELD_FIELD_ALIGNED_FIELD_LINE_TRACER_HPP
