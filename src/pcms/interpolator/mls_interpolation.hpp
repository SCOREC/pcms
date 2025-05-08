#ifndef MLS_RBF_OPTIONS_HPP
#define MLS_RBF_OPTIONS_HPP

#include "mls_interpolation_impl.hpp"

namespace pcms
{

/**
 * @brief Enumeration of supported radial basis functions (RBF) for MLS
 *interpolation.
 *
 * This enum specifies the type of radial basis function used to weight source
 *points in the Moving Least Squares (MLS) interpolation process.
 *
 * Members:
 * - RBF_GAUSSIAN:
 *     A Gaussian RBF: `exp(- a^2 * r^2)`, smooth and commonly used. Good for
 *localized support. `a` is a spreading/decay factor
 *
 * - RBF_C4:
 *     A compactly supported C4-continuous function. Useful for bounded support
 *and efficiency.
 *
 * - RBF_CONST:
 *     Constant basis function. Effectively uses uniform weights.
 *
 * - NO_OP:
 *     No operation. Disables RBF weighting — typically used when weights
 *     are externally defined or not needed.
 *
 * @note These are intended to be passed into function `mls_interpolation` to
 *control weighting behavior.
 */
enum class RadialBasisFunction : LO
{
  RBF_GAUSSIAN = 0,
  RBF_C4,
  RBF_CONST,
  NO_OP

};

/**
 * @brief Performs Moving Least Squares (MLS) interpolation at target points.
 *
 * This function computes interpolated values at a set of target coordinates
 * using Moving Least Squares (MLS) based on the provided source values and
 * coordinates. It supports different radial basis functions (RBFs), polynomial
 * degrees, dimension, optional regularization and tolerance
 *
 * @param source_values        A flat array of source data values. Length should
 * be `num_sources`.
 * @param source_coordinates   A flat array of source point coordinates. Length
 * should be `num_sources * dim`.
 * @param target_coordinates   A flat array of target point coordinates. Length
 * should be `num_targets * dim`.
 * @param support              A data structure holding neighbor information for
 * 								each target (in
 * CSR format).
 * @param dim                  Dimension
 * @param degree               Polynomial degree
 * @param bf                   The radial basis function used for weighting
 * (e.g., Gaussian, C4).
 * @param lambda               Optional regularization parameter (default is
 * 0.0). Helps with stability in ill-conditioned systems.
 * @param tol                  Optional solver tolerance (default is 1e-6).
 * Small singular values below this are discarded.
 *
 * @return A Write<Real> array containing the interpolated values at each target
 * point.
 *
 * @note
 * - All input arrays are expected to reside in device memory (e.g., Kokkos
 * device views).
 * - Ensure consistency in dimensions: coordinate arrays must be sized as
 * `num_points * dim`.
 * - The result array length will match the number of target points.
 */
Write<Real> mls_interpolation(const Omega_h::Reals source_values,
                              const Omega_h::Reals source_coordinates,
                              const Omega_h::Reals target_coordinates,
                              const SupportResults& support,
                              const Omega_h::LO& dim, const Omega_h::LO& degree,
                              RadialBasisFunction bf, double lambda = 0,
                              double tol = 1e-6, double decay_factor = 5.0);
} // namespace pcms
#endif
