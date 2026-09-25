#include "pcms/field/transformed_source.hpp"
#include <string>

namespace pcms
{

namespace detail
{

void ThrowIfAnyPointInvalid(const PointStatusView& status, const char* context)
{
  if (status.size() == 0) {
    return;
  }
  LO invalid = 0;
  Kokkos::parallel_reduce(
    "transformed_source_count_invalid",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(status.size())),
    KOKKOS_LAMBDA(const LO i, LO& count) {
      if (status(i) != PointStatus::Valid) {
        ++count;
      }
    },
    invalid);
  if (invalid > 0) {
    throw pcms_error(std::string(context) + ": " + std::to_string(invalid) +
                     " query points could not be mapped into the source "
                     "coordinate system");
  }
}

void RestoreFillRows(
  const Kokkos::View<Real**, DeviceMemorySpace>& pre_transform,
  Rank2View<Real, DeviceMemorySpace> values, Real fill_value)
{
  const auto num_cols = static_cast<int>(pre_transform.extent(1));
  Kokkos::parallel_for(
    "transformed_source_restore_fill_rows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(pre_transform.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      bool filled = true;
      for (int c = 0; c < num_cols; ++c) {
        const Real value = pre_transform(i, c);
        filled =
          filled && ((value == fill_value) ||
                     (Kokkos::isnan(value) && Kokkos::isnan(fill_value)));
      }
      if (filled) {
        for (int c = 0; c < num_cols; ++c) {
          values(i, c) = fill_value;
        }
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
  detail::ThrowIfAnyPointInvalid(bound->Status(),
                                 "TransformedSource::CreatePointEvaluator");
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
