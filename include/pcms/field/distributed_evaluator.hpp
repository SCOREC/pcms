#ifndef PCMS_FIELD_DISTRIBUTED_EVALUATOR_HPP
#define PCMS_FIELD_DISTRIBUTED_EVALUATOR_HPP

#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/distributed_geometry.hpp"

#include <memory>

namespace pcms
{

/**
 * Evaluates `local_source` over the whole decomposition of `geometry`: bind
 * the query points on `geometry` (collective) and evaluate fields of the local
 * source as usual; evaluation is collective too. The fields' ghost values must
 * be synchronized.
 * \param local_source an evaluator source on geometry's local geometry, e.g.
 *        a space created on it
 * \throws pcms_error if either is null or `local_source` is not on the local
 *         geometry
 */
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory>
CreateDistributedEvaluator(
  std::shared_ptr<const DistributedGeometry> geometry,
  std::shared_ptr<const PointEvaluatorFactory> local_source);

} // namespace pcms

#endif // PCMS_FIELD_DISTRIBUTED_EVALUATOR_HPP
