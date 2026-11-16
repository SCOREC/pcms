#ifndef PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
#define PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H

#include "pcms/field/field.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/field_data.h"
#include "pcms/field/value_view.hpp"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"

#include <Omega_h_mesh.hpp>
#include <memory>

namespace pcms
{

template <typename T>
class FieldEvaluatorFactory;

class PolynomialReconstructionFunctionSpace : public FunctionSpace
{
  // Passkey: only members (the Create/FromMesh factories) can create a Key, so
  // only they can construct a space, while std::make_shared can call the ctor.
  struct Key
  {
    explicit Key() = default;
  };

public:
  PolynomialReconstructionFunctionSpace(
    Key, std::shared_ptr<const FieldLayout> layout,
    std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept;

  [[nodiscard]] static std::shared_ptr<PolynomialReconstructionFunctionSpace>
  Create(CoordinateView<HostMemorySpace> coords, MLSOptions options = {});

  [[nodiscard]] static std::shared_ptr<PolynomialReconstructionFunctionSpace>
  FromMesh(Omega_h::Mesh& mesh, int source_entity_dim,
           std::shared_ptr<const CoordinateSystem> coordinate_system,
           MLSOptions options = {});

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override;

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(Type storage_type,
                                             ValueBasis basis) const override;

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override;

  [[nodiscard]] PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory_;
};

} // namespace pcms

#endif // PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
