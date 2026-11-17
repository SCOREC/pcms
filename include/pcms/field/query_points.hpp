#ifndef PCMS_FIELD_QUERY_POINTS_HPP
#define PCMS_FIELD_QUERY_POINTS_HPP

#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/differential_operator.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <map>
#include <memory>
#include <tuple>
#include <utility>
#include <variant>

namespace pcms
{

class CoordinateMap; // pcms/field/coordinate_map.hpp

/// @brief Evaluation at points that move between calls: one bound source on a
///        geometry, rebound in place by Locate, and the evaluators of the
///        functions evaluated on it, created on first use after each Locate.
///
/// Locate releases the evaluators before it rebinds, so Geometry::Rebind
/// reuses the bound source's storage unless the caller still holds it
/// through Bound(). Evaluators of an Interpolator hold their own bound
/// source and are not affected. Not thread-safe.
class QueryPoints
{
public:
  /// @param geometry geometry the points are located on
  /// @param policy out-of-bounds policy of the evaluation
  /// @param map map from the query coordinate system to the geometry's, or
  ///        null, as for Geometry::Bind
  /// @throws pcms_error if `geometry` is null
  explicit QueryPoints(std::shared_ptr<const Geometry> geometry,
                       OutOfBoundsPolicy policy = {},
                       std::shared_ptr<const CoordinateMap> map = nullptr);

  /// @brief Binds the points on the first call and rebinds them after,
  ///        reusing the storage when they do not outgrow it; drops the
  ///        evaluators.
  /// @param coords query points, as for Geometry::Bind; may be released after
  /// @throws pcms_error as Geometry::Bind; no points are located after a throw
  void Locate(CoordinateView<DeviceMemorySpace> coords);

  /// @brief `op` of `function` at the located points.
  ///
  /// The evaluator for (function's space, op) is created on first use after
  /// Locate and reused.
  /// @param out output, shaped and tagged as PointEvaluator::Evaluate requires
  /// @throws pcms_error if Locate was not called (or threw), the space is not
  ///         on this geometry, or the space does not support `op`; otherwise
  ///         as PointEvaluator::Evaluate
  template <typename T>
  void Evaluate(const Function<T>& function,
                ValueView<T, DeviceMemorySpace> out,
                DifferentialOperator op = DifferentialOperator::Value)
  {
    if (bound_ == nullptr) {
      throw pcms_error("QueryPoints::Evaluate: no points are located; call "
                       "Locate first");
    }
    Key key{function.GetSpacePtr(), op, detail::TypeEnumFromType<T>()};
    auto it = evaluators_.find(key);
    if (it == evaluators_.end()) {
      auto evaluator =
        function.GetSpace().template CreatePointEvaluator<T>(bound_, op);
      it = evaluators_.emplace(std::move(key), std::move(evaluator)).first;
    }
    std::get<std::unique_ptr<PointEvaluator<T>>>(it->second)
      ->Evaluate(function, out);
  }

  /// @return the bound source of the located points (FoundPoints(),
  ///         MissingPoints()), or null before Locate. A copy held past the
  ///         next Locate keeps these points; that Locate then binds anew
  ///         instead of reusing the storage.
  [[nodiscard]] const std::shared_ptr<const BoundSource>& Bound() const noexcept
  {
    return bound_;
  }

private:
  // the space is held so that its address is not reused by another space
  using Key = std::tuple<std::shared_ptr<const FunctionSpace>,
                         DifferentialOperator, Type>;

  std::shared_ptr<const Geometry> geometry_;
  OutOfBoundsPolicy policy_;
  std::shared_ptr<const CoordinateMap> map_;
  std::shared_ptr<const BoundSource> bound_;
  std::map<Key, PointEvaluatorVariant> evaluators_;
};

/// @brief `op` of `function` at `points`, in one call.
///
/// Binds `points` on every call: for repeated evaluation at moving points use
/// QueryPoints, and for fixed points bind once and keep the evaluators.
/// @param points query points, as for Geometry::Bind
/// @param out output, shaped and tagged as PointEvaluator::Evaluate requires
/// @param policy out-of-bounds policy
/// @param map map from the query coordinate system to the geometry's, or null
/// @throws pcms_error if the function's space does not evaluate at points, or
///         as QueryPoints::Locate and QueryPoints::Evaluate
template <typename T>
void Evaluate(const Function<T>& function,
              CoordinateView<DeviceMemorySpace> points,
              ValueView<T, DeviceMemorySpace> out,
              DifferentialOperator op = DifferentialOperator::Value,
              OutOfBoundsPolicy policy = {},
              std::shared_ptr<const CoordinateMap> map = nullptr)
{
  QueryPoints query(function.GetSpace().GetGeometry(), policy, std::move(map));
  query.Locate(std::move(points));
  query.Evaluate(function, std::move(out), op);
}

} // namespace pcms

#endif // PCMS_FIELD_QUERY_POINTS_HPP
