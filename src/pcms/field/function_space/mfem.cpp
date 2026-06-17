#include "pcms/field/function_space/mfem.h"

#include "pcms/field/data/mfem.h"
#include "pcms/utility/assert.h"

#include <memory>

namespace pcms
{

MFEMFunctionSpace::MFEMFunctionSpace(std::shared_ptr<const MFEMLayout> layout,
                                     mfem::ParFiniteElementSpace& pfes,
                                     mfem::ParGridFunction& gf,
                                     CoordinateSystem coordinate_system) noexcept
  : layout_(std::move(layout)),
    pfes_(pfes),
    gf_(gf),
    coordinate_system_(coordinate_system)
{
}

MFEMFunctionSpace MFEMFunctionSpace::FromMesh(
  mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
  mfem::ParGridFunction& gf, CoordinateSystem coordinate_system)
{
  auto layout =
    std::make_shared<const MFEMLayout>(pmesh, pfes, coordinate_system);
  return MFEMFunctionSpace(std::move(layout), pfes, gf, coordinate_system);
}

std::shared_ptr<const FieldLayout> MFEMFunctionSpace::GetLayout() const noexcept
{
  return layout_;
}

CoordinateSystem MFEMFunctionSpace::GetCoordinateSystem() const noexcept
{
  return coordinate_system_;
}

FieldVariant MFEMFunctionSpace::CreateFieldImpl(Type value_type,
                                                FieldMetadata metadata) const
{
  if (value_type != Type::Real) {
    throw pcms_error("MFEMFunctionSpace only supports double (Real) fields");
  }
  auto data = std::make_unique<MFEMVertexFieldData>(pfes_, gf_, metadata);
  return WrapField<Real>(layout_, std::move(data), /*evaluator_factory=*/nullptr);
}

FieldVariant MFEMFunctionSpace::CreateFieldImpl(FieldDataVariant data) const
{
  auto* real_data = std::get_if<std::unique_ptr<FieldData<Real>>>(&data);
  if (real_data == nullptr || *real_data == nullptr) {
    throw pcms_error(
      "MFEMFunctionSpace only supports double (Real) field data");
  }
  return WrapField<Real>(layout_, std::move(*real_data),
                         /*evaluator_factory=*/nullptr);
}

PointEvaluatorVariant MFEMFunctionSpace::CreatePointEvaluatorImpl(
  Type, const EvaluationRequest&) const
{
  throw pcms_error("MFEMFunctionSpace is communication-only; point evaluation "
                   "is not supported");
}

} // namespace pcms
