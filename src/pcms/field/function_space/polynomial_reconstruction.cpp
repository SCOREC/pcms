#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/layout/omega_h_entity.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/field/evaluator/point_cloud.h"
#include "pcms/field/data/simple.h"
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/localization/point_cloud_localization.h"
#include "pcms/localization/mesh_localization.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <string>

namespace pcms
{

PolynomialReconstructionFunctionSpace::PolynomialReconstructionFunctionSpace(
  Key, std::shared_ptr<const FieldLayout> layout,
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept
  : layout_(std::move(layout)), evaluator_factory_(std::move(evaluator_factory))
{
}

std::shared_ptr<PolynomialReconstructionFunctionSpace>
PolynomialReconstructionFunctionSpace::Create(
  CoordinateView<HostMemorySpace> coords, MLSOptions options)
{
  const auto values = coords.GetValues();
  const int dim = static_cast<int>(values.extent(1));
  Kokkos::View<const Real**, Kokkos::HostSpace> host_view(
    values.data_handle(), values.extent(0), values.extent(1));
  auto device_view = Kokkos::View<Real**>("device_view", host_view.extent(0),
                                          host_view.extent(1));
  DeepCopyMismatchLayouts(device_view, host_view);
  auto pc_layout = std::make_shared<PointCloudLayout>(
    dim, device_view, coords.GetCoordinateSystem());
  auto localization =
    std::make_shared<PointCloudLocalizationFactory>(pc_layout, options);
  auto eval_factory = std::make_shared<PointCloudEvaluatorFactory>(
    pc_layout, localization, options);
  return std::make_shared<PolynomialReconstructionFunctionSpace>(
    Key{}, pc_layout, std::move(eval_factory));
}

std::shared_ptr<PolynomialReconstructionFunctionSpace>
PolynomialReconstructionFunctionSpace::FromMesh(
  Omega_h::Mesh& mesh, int source_entity_dim,
  std::shared_ptr<const CoordinateSystem> coordinate_system, MLSOptions options)
{
  if (source_entity_dim < 0 || source_entity_dim > mesh.dim()) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace::FromMesh: source_entity_dim is "
      "out of range");
  }

  auto mesh_layout = std::make_shared<OmegaHEntityLayout>(
    mesh, source_entity_dim, 1, std::move(coordinate_system));
  // MLS assumes unweigthed Euclidean in distance calculation
  if (!HasIdentityMetric(*mesh_layout->GetCoordinateSystem())) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace::FromMesh: MLS requires a "
      "coordinate system whose metric is the identity; got '" +
      std::string(mesh_layout->GetCoordinateSystem()->Kind()) + "'");
  }
  auto localization = std::make_shared<AdjacencyLocalizationFactory>(
    mesh, source_entity_dim, options);
  auto eval_factory = std::make_shared<PointCloudEvaluatorFactory>(
    mesh_layout, localization, options);
  return std::make_shared<PolynomialReconstructionFunctionSpace>(
    Key{}, mesh_layout, std::move(eval_factory));
}

std::shared_ptr<const FieldLayout>
PolynomialReconstructionFunctionSpace::GetLayout() const noexcept
{
  return layout_;
}

FieldVariant PolynomialReconstructionFunctionSpace::CreateFieldImpl(
  Type storage_type, ValueBasis basis) const
{
  return apply_to_type(storage_type, [&](auto tag) -> FieldVariant {
    using T = typename decltype(tag)::type;
    if constexpr (!std::is_same_v<T, double>) {
      throw pcms_error(
        "PolynomialReconstructionFunctionSpace: only double (Real) is "
        "supported");
    } else {
      return WrapField<double>(
        layout_,
        std::make_unique<SimpleFieldData<double>>(layout_, std::move(basis)));
    }
  });
}

FieldVariant PolynomialReconstructionFunctionSpace::CreateFieldImpl(
  FieldDataVariant data) const
{
  if (!std::holds_alternative<std::unique_ptr<FieldData<double>>>(data)) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace: only double (Real) is "
      "supported");
  }
  auto fd = std::move(std::get<std::unique_ptr<FieldData<double>>>(data));
  PCMS_ALWAYS_ASSERT(fd != nullptr);
  if (dynamic_cast<const SimpleFieldData<double>*>(fd.get()) == nullptr) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace::CreateField: requires "
      "SimpleFieldData<double>");
  }
  if (fd->GetDOFHolderDataHost().size() !=
      detail::ExpectedFlatFieldDataSize(*layout_)) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace::CreateField: field data size "
      "does not match layout");
  }
  return WrapField<double>(layout_, std::move(fd));
}

PointEvaluatorVariant
PolynomialReconstructionFunctionSpace::CreatePointEvaluatorImpl(
  Type value_type, const EvaluationRequest& request) const
{
  if (value_type != Type::Real) {
    throw pcms_error(
      "PolynomialReconstructionFunctionSpace: point evaluation only supports "
      "double (Real)");
  }
  PCMS_ALWAYS_ASSERT(evaluator_factory_ != nullptr);
  return evaluator_factory_->CreatePointEvaluator(request);
}

} // namespace pcms
