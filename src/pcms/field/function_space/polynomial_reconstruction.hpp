#ifndef PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
#define PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H

#include "pcms/field/field.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"

#include <Omega_h_mesh.hpp>
#include <memory>

namespace pcms::field
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

} // namespace pcms::field

#endif // PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
