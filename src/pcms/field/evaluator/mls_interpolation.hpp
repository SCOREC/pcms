#ifndef PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP
#define PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP

#include <Omega_h_array.hpp>
#include "pcms/localization/mls_options.hpp"
#include "pcms/localization/support_results.hpp"

namespace pcms
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
 * @return A Write<Real> array containing interpolated values at each target
 * point. Length matches the number of target points.
 *
 * @note
 * - All input arrays must reside in device memory (e.g., Kokkos device views).
 * - Coordinate arrays must be correctly sized as `num_points * dim`.
 * - Choose `lambda` and `decay_factor` carefully to balance accuracy and
 * stability: a larger decay factor often requires a slightly higher `lambda`.
 */
Omega_h::Write<Omega_h::Real> mls_interpolation(
  const Omega_h::Reals source_values, const Omega_h::Reals source_coordinates,
  const Omega_h::Reals target_coordinates, const SupportResults& support,
  const Omega_h::LO& dim, const Omega_h::LO& degree, RadialBasisFunction bf,
  double lambda = 0, double tol = 1e-6, double decay_factor = 5.0);

} // namespace pcms
#endif // PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_HPP
