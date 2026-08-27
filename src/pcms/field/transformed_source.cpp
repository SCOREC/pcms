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
  const PointEvaluatorFactory& source,
  std::shared_ptr<const CoordinateMap> to_source, ValueBasis source_basis,
  ValueBasis target_basis)
  : source_(&source),
    to_source_(std::move(to_source)),
    source_basis_(std::move(source_basis)),
    target_basis_(std::move(target_basis))
{
  if (to_source_ == nullptr) {
    throw pcms_error("TransformedSource: to_source must not be null");
  }
  if (!SameCoordinateSystem(to_source_->GetTargetCoordinateSystem(),
                            source.GetCoordinateSystem())) {
    throw pcms_error(
      "TransformedSource: to_source maps into '" +
      std::string(to_source_->GetTargetCoordinateSystem()->Kind()) +
      "' but the source's coordinate system is '" +
      std::string(source.GetCoordinateSystem()->Kind()) + "'");
  }
  system_ = to_source_->GetSourceCoordinateSystem();
  detail::ValidateValueBasis(source_basis_);
  detail::ValidateValueBasis(target_basis_);
  if (source_basis_.Rank() != target_basis_.Rank()) {
    throw pcms_error("TransformedSource: source and target bases have "
                     "different ranks");
  }
  if (source_basis_.Rank() == 0) {
    return;
  }
  if (!(source_basis_.variance == target_basis_.variance)) {
    throw pcms_error("TransformedSource: source and target bases have "
                     "different variances; index raising and lowering is "
                     "not supported");
  }
  if (source_basis_.Rank() > 1) {
    throw pcms_error("TransformedSource: rank-" +
                     std::to_string(source_basis_.Rank()) +
                     " transformations are not implemented");
  }
  if (!SameCoordinateSystem(source_basis_.system,
                            to_source_->GetTargetCoordinateSystem())) {
    throw pcms_error("TransformedSource: source_basis must be in the "
                     "source's coordinate system");
  }
  if (!SameCoordinateSystem(target_basis_.system, system_)) {
    throw pcms_error("TransformedSource: target_basis must be in the "
                     "to_source map's source coordinate system");
  }
}

template <typename T>
std::unique_ptr<PointEvaluator<T>> TransformedSource::Create(
  const EvaluationRequest& request) const
{
  PCMS_FUNCTION_TIMER;
  MappedPoints mapped = to_source_->Map(request.coords);
  detail::ThrowIfAnyPointInvalid(mapped.status,
                                 "TransformedSource::CreatePointEvaluator");
  if constexpr (!std::is_same_v<T, Real>) {
    if (source_basis_.Rank() > 0) {
      throw pcms_error("TransformedSource::CreatePointEvaluator: basis "
                       "transformations require T == Real");
    }
  }
  auto inner = source_->CreatePointEvaluator<T>(
    EvaluationRequest::FromCoordinates(mapped.View(), request.policy));
  std::unique_ptr<BoundBasisTransformation> law;
  if (source_basis_.Rank() > 0) {
    if constexpr (std::is_same_v<T, Real>) {
      law = to_source_->MakeBasisTransformation(request.coords, mapped);
      if (law == nullptr) {
        throw pcms_error("TransformedSource::CreatePointEvaluator: the "
                         "coordinate map provides no basis transformation");
      }
      if (!SameValueBasis(law->GetSourceBasis(), source_basis_) ||
          !SameValueBasis(law->GetTargetBasis(), target_basis_)) {
        throw pcms_error("TransformedSource::CreatePointEvaluator: the "
                         "coordinate map's basis transformation does not "
                         "bridge source_basis to target_basis");
      }
    }
  }
  return std::make_unique<TransformedPointEvaluator<T>>(
    std::move(inner), std::move(mapped), std::move(law), source_basis_,
    target_basis_, request.policy);
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
