#ifndef PCMS_FUNCTION_SPACE_MFEM_H
#define PCMS_FUNCTION_SPACE_MFEM_H

#include "pcms/field/coordinate_system.h"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/mfem.h"

#include <mfem.hpp>

#include <memory>

namespace pcms
{

// Communication-only function space for an MFEM order-1 H1 (vertex) scalar
// field. It exists to wrap the live mfem::ParGridFunction into a Field that the
// coupler can send/receive; it does not support point evaluation.
//
// CreateField<Real>() returns a Field whose data is bound to the grid function
// passed at construction, so coupled get/set operate on the solver's state.
class MFEMFunctionSpace : public FunctionSpace
{
public:
  [[nodiscard]] static MFEMFunctionSpace FromMesh(
    mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
    mfem::ParGridFunction& gf, CoordinateSystem coordinate_system);

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override;

  [[nodiscard]] CoordinateSystem GetCoordinateSystem() const noexcept override;

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(
    Type value_type, FieldMetadata metadata) const override;

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override;

  [[nodiscard]] PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  MFEMFunctionSpace(std::shared_ptr<const MFEMLayout> layout,
                    mfem::ParFiniteElementSpace& pfes,
                    mfem::ParGridFunction& gf,
                    CoordinateSystem coordinate_system) noexcept;

  std::shared_ptr<const MFEMLayout> layout_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;
  CoordinateSystem coordinate_system_;
};

} // namespace pcms

#endif // PCMS_FUNCTION_SPACE_MFEM_H
