#ifndef PCMS_LAGRANGE_FIELD_FACTORY_H
#define PCMS_LAGRANGE_FIELD_FACTORY_H

#include <Omega_h_mesh.hpp>
#include "pcms/configuration.h"
#include "pcms/field/field.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/function_space.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_metadata.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/uniform_grid.h"

#include <functional>
#include <memory>
#include <string>
#include <utility>

namespace pcms
{

class LagrangeFunctionSpace : public FunctionSpace
{
  // Passkey: only members (the From* factories) can create a Key, so only they
  // can construct a space, while std::make_shared can still call the ctor.
  struct Key
  {
    explicit Key() = default;
  };

public:
  LagrangeFunctionSpace(
    Key, std::shared_ptr<const FieldLayout> layout,
    std::function<FieldDataVariant(Type, FieldMetadata)> create_field_data_fn,
    std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept;

  enum class Backend
  {
    MeshFields,
    OmegaH
  };

#ifdef PCMS_ENABLE_MESHFIELDS
  static constexpr Backend DefaultBackend = Backend::MeshFields;
#else
  static constexpr Backend DefaultBackend = Backend::OmegaH;
#endif

  // Unstructured mesh — dispatches to MeshFields or native Omega_h backend
  [[nodiscard]] static std::shared_ptr<LagrangeFunctionSpace> FromMesh(
    Omega_h::Mesh& mesh, int order, int num_components,
    std::shared_ptr<const CoordinateSystem> mesh_coordinate_system,
    std::string global_id_name = "global", Backend backend = DefaultBackend,
    std::string layout_name = "");

  [[nodiscard]] static std::shared_ptr<LagrangeFunctionSpace> FromMesh(
    Omega_h::Mesh& mesh, int order, int num_components,
    std::shared_ptr<const CoordinateSystem> mesh_coordinate_system,
    Omega_h::Read<Omega_h::I8> owned_mask,
    std::string global_id_name = "global", Backend backend = DefaultBackend,
    std::string layout_name = "");

  // Structured uniform grid — order-1 H1-conforming nodal field on a regular
  // grid
  [[nodiscard]] static std::shared_ptr<LagrangeFunctionSpace> FromUniformGrid(
    const UniformGrid<2>& grid, int num_components,
    std::shared_ptr<const CoordinateSystem> mesh_coordinate_system,
    int order = 1, std::string layout_name = "");

  [[nodiscard]] static std::shared_ptr<LagrangeFunctionSpace> FromUniformGrid(
    const UniformGrid<3>& grid, int num_components,
    std::shared_ptr<const CoordinateSystem> mesh_coordinate_system,
    int order = 1, std::string layout_name = "");

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override;

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(Type value_type, FieldMetadata metadata) const override;

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override;

  [[nodiscard]] PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::function<FieldDataVariant(Type, FieldMetadata)> create_field_data_fn_;
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory_;
};

} // namespace pcms

#endif // PCMS_LAGRANGE_FIELD_FACTORY_H
