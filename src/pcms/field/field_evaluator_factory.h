#ifndef PCMS_FIELD_EVALUATOR_FACTORY_H
#define PCMS_FIELD_EVALUATOR_FACTORY_H

#include "field_layout.h"
#include "pcms/field/coordinate_view.hpp"
#include "evaluation_request.h"
#include "out_of_bounds_policy.h"
#include "point_evaluator.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/discretization/discretization.h"
#include <memory>

namespace pcms
{

// FieldEvaluatorFactory<T> owns backend-specific knowledge needed to evaluate a
// field (mesh geometry, model architecture, analytic parameters, etc.). It is
// reusable across any number of query point sets and FieldData objects and is
// created once per field configuration.
//
// Contract:
//   - A factory instance is tied to one backend/layout configuration.
//   - CreatePointEvaluator(...) either returns a fully usable evaluator or
//     throws pcms_error immediately; unsupported operations must not be
//     deferred to PointEvaluator::Evaluate().
//   - A PointEvaluator created by this factory may be reused with multiple
//     FieldData objects only when they are compatible with this factory's
//     layout/backend contract.
//   - OutOfBoundsPolicy is fixed at evaluator creation time.
//
// FieldEvaluatorFactory<T> is a standalone abstract type rather than being
// merged into concrete field factories (LagrangeFunctionSpace, etc.). This is
// intentional: field families such as analytic functions and ML surrogates can
// be evaluated at arbitrary points but have no meaningful field construction
// operation. Concrete field factories implement both this interface and their
// own field creation interface.
//
// Concrete field factories store a FieldEvaluatorFactory<T> by composition and
// expose CreatePointEvaluator(request) as a convenience that delegates to the
// internal factory.
//
// Example (sketch of a concrete field factory):
//
//   class LagrangeFunctionSpace {
//   public:
//     std::unique_ptr<FieldData<Real>>   CreateFieldReal() const;
//     std::shared_ptr<const FieldLayout> GetLayout() const;
//
//     // convenience — delegates to the internal FieldEvaluatorFactory:
//     std::unique_ptr<PointEvaluator<Real>> CreatePointEvaluator(
//       const EvaluationRequest& request) const;
//
//   private:
//     std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory_;
//   };
template <typename T>
class FieldEvaluatorFactory
{
public:
  // Returns the layout this factory was constructed for. Used by callers to
  // verify FieldData compatibility before calling Evaluate.
  virtual const FieldLayout& GetLayout() const = 0;

  virtual bool HasDOFHolderCoordinates() const = 0;

  virtual CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const = 0;

  // Whether this factory supports OutOfBoundsMode::NearestBoundary.
  // Backends that return false will throw a descriptive error if
  // NearestBoundary is requested via OutOfBoundsPolicy.
  virtual bool SupportsNearestBoundary() const = 0;

  // Performs all localization work for the given evaluation request and
  // returns a PointEvaluator that caches the resolved evaluation state. The
  // PointEvaluator may be reused across Evaluate calls as long as the query
  // points and underlying geometry/layout remain unchanged. Implementations may
  // discard construction-time request metadata once the cached evaluator state
  // has been built. On misuse or unsupported capability, implementations
  // should throw pcms_error.
  virtual std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    const EvaluationRequest& request) const = 0;

  virtual ~FieldEvaluatorFactory() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_FACTORY_H
