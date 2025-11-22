#ifndef MLS_RBF_OPTIONS_HPP
#define MLS_RBF_OPTIONS_HPP

#include "mls_interpolation_impl.hpp"

namespace pcms
{
/**
 * @brief Enumeration of supported radial basis functions (RBFs) for MLS
 * interpolation.
 *
 * Specifies the kernel type used for computing weights in Moving Least Squares
 * (MLS) interpolation. Each RBF defines a different spatial weighting behavior,
 * influencing locality, smoothness, and conditioning of the interpolant.
 *
 * ### Values:
 *
 * - **RBF_GAUSSIAN (0)**
 *   \f$ \phi(r) = \exp\left(-\left(\frac{a r}{r_c}\right)^2\right) \f$
 *   Infinitely smooth and compactly supported.
 *   `a` is the shape (decay) parameter — smaller values give wider influence.
 *   Good for smooth interpolants; may cause ill-conditioning for large `a`.
 *
 * - **RBF_C4 (1)**
 *   Compactly supported, \f$ C^4\f$-continuous polynomial kernel.
 *   Does not require a shape parameter.
 *   Efficient, stable, and good for local approximations.
 *
 * - **RBF_CONST (2)**
 *   Uniform weights within support:
 *   \f$ \phi(r) = 1 \text{ for } r < r_c \f$, otherwise 0.
 *   Equivalent to averaging.
 *
 * - **NO_OP (3)**
 *   Disables internal RBF weighting logic.
 *   Always returns \f$ \phi = 1 \f$ regardless of distance.
 *   Useful for methods like SPR patching
 *
 * - **RBF_MULTIQUADRIC (4)**
 *   \f$ \phi(r) = \sqrt{1 + \left(\frac{a r}{r_c}\right)^2} \f$
 *   Smooth and locally supported.
 *   `a` controls the width; higher values increase conditioning issues.
 *
 * - **RBF_INVMULTIQUADRIC (5)**
 *   \f$ \phi(r) = \frac{1}{\sqrt{1 + \left(\frac{a r}{r_c}\right)^2}} \f$
 *   Compactly supported and generally more stable than MQ.
 *   `a` adjusts the decay; larger `a` leads to faster fall-off.
 *
 * - **RBF_THINPLATESPLINE (6)**
 *   \f$ \phi(r) = \left(\frac{a r}{r_c}\right)^2 \log\left(\frac{a
 * r}{r_c}\right) \f$ Classic 2D surface fitting kernel. No explicit shape
 * parameter is required, but `a` can scale smoothness.
 *
 * - **RBF_CUBIC (7)**
 *   \f$ \phi(r) = \left(\frac{a r}{r_c}\right)^3 \f$
 *   Smooth with compact support.
 *   Simple and effective for local interpolation.
 *
 * ---
 *
 * @note
 * - The shape parameter `a` (also called decay or spread factor) is required
 * for Gaussian, MQ, IMQ, ThinPlateSpline, and Cubic kernels.
 * - All kernels are compactly supported within the cutoff radius \f$ r_c \f$.
 * - A good starting value for `a` is 3–5. Avoid `a = 0` (leads to constant
 * behavior).
 * - In cases of ill-conditioning, use regularization (`lambda > 0`).
 *
 * @see mls_interpolation()
 */
enum class RadialBasisFunction : LO
{
  RBF_GAUSSIAN = 0,
  RBF_C4,
  RBF_CONST,
  NO_OP,
  RBF_MULTIQUADRIC,
  RBF_INVMULTIQUADRIC,
  RBF_THINPLATESPLINE,
  RBF_CUBIC
};

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
#endif
