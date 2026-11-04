#ifndef PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H
#define PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H

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
enum class RadialBasisFunction
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


/// Configuration for Moving Least Squares evaluation.
/// Passed to PolynomialReconstructionFunctionSpace::Create and
/// PointCloudEvaluatorFactory.
struct MLSOptions
{
  double radius = 0.5;
  unsigned min_req_supports = 10;
  unsigned degree = 3;
  bool adapt_radius = true;
  double lambda = 0.0;
  double tol = 1e-6;
  double decay_factor = 5.0;
  RadialBasisFunction basis = RadialBasisFunction::RBF_GAUSSIAN;
};

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H
