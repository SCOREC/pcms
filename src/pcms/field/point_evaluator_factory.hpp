#ifndef PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP
#define PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/field_factory.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/types.h"
#include <memory>
#include <string>
#include <variant>

namespace pcms
{

class PointEvaluatorFactory
{
public:
  /// Coordinate system that query points handed to CreatePointEvaluator must
  /// be expressed in.
  [[nodiscard]] virtual const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const = 0;

  /// Creates a point evaluator bound to the request's query points.
  /// @param request query coordinates, out-of-bounds policy, and optional
  ///        provenance; the coordinates must be in GetCoordinateSystem()
  template <typename T>
  [[nodiscard]] std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    const EvaluationRequest& request) const;

  virtual ~PointEvaluatorFactory() noexcept = default;

protected:
  virtual PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const = 0;
};

template <typename T>
std::unique_ptr<PointEvaluator<T>> PointEvaluatorFactory::CreatePointEvaluator(
  const EvaluationRequest& request) const
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  if (!SameCoordinateSystem(request.coords.GetCoordinateSystem(),
                            GetCoordinateSystem())) {
    throw pcms_error(
      "CreatePointEvaluator: query coordinate system '" +
      std::string(request.coords.GetCoordinateSystem()->Kind()) +
      "' is not the factory's coordinate system '" +
      std::string(GetCoordinateSystem()->Kind()) +
      "'; map the points with a CoordinateMap or evaluate through a "
      "TransformedSource");
  }
  return std::get<std::unique_ptr<PointEvaluator<T>>>(
    CreatePointEvaluatorImpl(TypeEnumFromType<T>(), request));
}

} // namespace pcms

#endif // PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP
