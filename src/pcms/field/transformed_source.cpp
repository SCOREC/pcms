#include "pcms/field/transformed_source.hpp"
#include "pcms/field/function_space.h"
#include <string>

namespace pcms
{

namespace detail
{

void FillRows(Rank2View<Real, DeviceMemorySpace> values,
              const Kokkos::View<const LO*, DeviceMemorySpace>& rows,
              Real fill_value)
{
  const auto num_cols = static_cast<int>(values.extent(1));
  Kokkos::parallel_for(
    "transformed_source_fill_rows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(rows.extent(0))),
    KOKKOS_LAMBDA(const LO k) {
      for (int c = 0; c < num_cols; ++c) {
        values(rows(k), c) = fill_value;
      }
    });
}

} // namespace detail

TransformedSource::TransformedSource(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const CoordinateMap> to_source)
  : source_(std::move(source)), to_source_(std::move(to_source))
{
  if (source_ == nullptr) {
    throw pcms_error("TransformedSource: source must not be null");
  }
  if (to_source_ == nullptr) {
    throw pcms_error("TransformedSource: to_source must not be null");
  }
  if (!SameCoordinateSystem(to_source_->GetTargetCoordinateSystem(),
                            source_->GetCoordinateSystem())) {
    throw pcms_error(
      "TransformedSource: to_source maps into '" +
      std::string(to_source_->GetTargetCoordinateSystem()->Kind()) +
      "' but the source's coordinate system is '" +
      std::string(source_->GetCoordinateSystem()->Kind()) +
      "'; expected a map ending in '" +
      std::string(source_->GetCoordinateSystem()->Kind()) + "'");
  }
  system_ = to_source_->GetSourceCoordinateSystem();
}

template <typename T>
std::unique_ptr<PointEvaluator<T>> TransformedSource::Create(
  const EvaluationRequest& request) const
{
  PCMS_FUNCTION_TIMER;
  auto bound = to_source_->Bind(request.coords);
  auto inner = source_->CreatePointEvaluator<T>(
    EvaluationRequest::FromCoordinates(bound->MappedPoints(), request.policy));
  return std::make_unique<TransformedPointEvaluator<T>>(
    std::move(bound), std::move(inner), request.policy);
}

PointEvaluatorVariant TransformedSource::CreatePointEvaluatorImpl(
  Type value_type, const EvaluationRequest& request) const
{
  switch (value_type) {
    case Type::Real: return Create<Real>(request);
    case Type::LO: return Create<LO>(request);
    case Type::GO: return Create<GO>(request);
    case Type::Int8: return Create<int8_t>(request);
    case Type::Float: return Create<float>(request);
  }
  throw pcms_error("TransformedSource::CreatePointEvaluator: unsupported "
                   "value type");
}

} // namespace pcms
