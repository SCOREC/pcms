#ifndef PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H
#define PCMS_FIELD_EVALUATOR_MLS_OPTIONS_H

#include "pcms/field/evaluator/mls_interpolation.hpp"

namespace pcms
{

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
