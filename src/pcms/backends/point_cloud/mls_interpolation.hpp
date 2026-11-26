#ifndef PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP
#define PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP

#include "pcms/field/mls_options.hpp"
#include "pcms/geometry/supports.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

#include <Kokkos_Core.hpp>

namespace pcms::detail
{

/**
 * @brief Performs Moving Least Squares (MLS) interpolation at target points.
 *
 * This function computes interpolated values at a set of target coordinates
 * using Moving Least Squares (MLS) based on the provided source values and
 * coordinates. It supports different radial basis functions (RBFs), polynomial
 * degrees, spatial dimensions, optional regularization, and solver tolerance.
 *
 * @param source_values        A flat array of source data values. Length should
 *                             match the number of source points
 * (`num_sources`).
 * @param source_coordinates   A flat array of source point coordinates of size
 *                             `num_sources * dim`.
 * @param target_coordinates   A flat array of target point coordinates of size
 *                             `num_targets * dim`.
 * @param support              A structure containing neighbor information for
 *                             each target point in Compressed Sparse Row (CSR)
 * format.
 * @param dim                  Spatial dimension of the coordinate data (e.g., 2
 * or 3).
 * @param degree               Degree of the polynomial basis used in MLS.
 * @param bf                   The radial basis function (RBF) used for
 * weighting (e.g., Gaussian, C4).
 * @param lambda               Tikhonov regularization parameter (default: 0.0).
 *                             A small positive value improves numerical
 * stability for ill-conditioned local systems.
 *                             - Well-conditioned systems: 1e-8 to 1e-6
 *                             - Mildly ill-conditioned: 1e-5 to 1e-3
 *                             - Highly ill-conditioned: 1e-2 to 1e-1
 * @param tol                  Optional solver tolerance (default: 1e-6).
 * Singular values below this threshold are ignored.
 * @param decay_factor         Controls the spread of the RBF influence width
 * (default: 5.0).
 *                             - 1–3  → wide influence (smooth interpolation)
 *                             - 3–8  → moderate locality (balanced,
 * recommended)
 *                             - 8–20 → highly local behavior (sharp detail,
 * possible conditioning issues)
 *
 * @return A device view containing interpolated values at each target
 * point. Length matches the number of target points. Rows without supports
 * are zero; the point evaluator applies its missing-point policy to them.
 *
 * @note
 * - All input arrays must reside in device memory (e.g., Kokkos device views).
 * - Coordinate arrays must be correctly sized as `num_points * dim`.
 * - Choose `lambda` and `decay_factor` carefully to balance accuracy and
 * stability: a larger decay factor often requires a slightly higher `lambda`.
 */
Kokkos::View<Real*, DeviceMemorySpace> mls_interpolation(
  Kokkos::View<const Real*, DeviceMemorySpace> source_values,
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, RadialBasisFunction bf,
  double lambda = 0, double tol = 1e-6, double decay_factor = 5.0);

/**
 * @brief The coefficient rows of the MLS interpolation at the target points:
 * the value at target i is the sum over its supports k (the CSR row of
 * `support`) of coefficients[k] * source_values[support.supports_idx[k]]. The
 * rows do not depend on the source values, so an evaluator that evaluates the
 * same targets repeatedly computes them once. Arguments as for
 * mls_interpolation.
 *
 * @return one coefficient per support entry, aligned with
 * support.supports_idx; targets without supports have no entries
 */
Kokkos::View<Real*, DeviceMemorySpace> mls_coefficients(
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, RadialBasisFunction bf,
  double lambda = 0, double tol = 1e-6, double decay_factor = 5.0);

} // namespace pcms::detail
#endif // PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP
