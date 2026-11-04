#ifndef PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H
#define PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H

#include "pcms/field/evaluator/mls_interpolation.hpp"

namespace pcms
{

/// Configuration of the Moving Least Squares reconstruction. The supports it
/// is computed on come from the geometry's search (NeighbourhoodOptions).
/// Passed to CreatePolynomialReconstructionSpace and
/// PointCloudEvaluatorFactory.
struct MLSOptions
{
  unsigned degree = 3;
  double lambda = 0.0;
  double tol = 1e-6;
  double decay_factor = 5.0;
  RadialBasisFunction basis = RadialBasisFunction::RBF_GAUSSIAN;
};

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H
