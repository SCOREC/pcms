#include "pcms/field/function_space/lagrange.h"
#include "pcms/configuration.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"
#ifdef PCMS_ENABLE_MESHFIELDS
#include "pcms/field/layout/mesh_fields.h"
#include "pcms/field/evaluator/mesh_fields.h"
#include "pcms/field/data/mesh_fields.h"
#endif
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/evaluator/omega_h_lagrange.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/field/evaluator/uniform_grid.h"
#include "pcms/field/data/simple.h"
#include "pcms/utility/uniform_grid.h"

#include <stdexcept>

namespace pcms
{

namespace
{

template <typename T>
void ValidateLagrangeWrappedFieldData(const FieldLayout& layout,
                                      const FieldData<T>& data)
{
#ifdef PCMS_ENABLE_MESHFIELDS
  if (dynamic_cast<const MeshFieldsAdapterLayout*>(&layout) != nullptr) {
    if (dynamic_cast<const MeshFieldsFieldData<T>*>(&data) == nullptr) {
      throw pcms_error(
        "LagrangeFunctionSpace::CreateField: MeshFields layout requires "
        "MeshFieldsFieldData");
    }
  } else
#endif
  {
    if (dynamic_cast<const SimpleFieldData<T>*>(&data) == nullptr) {
      throw pcms_error(
        "LagrangeFunctionSpace::CreateField: this backend requires "
        "SimpleFieldData");
    }
  }

  if (data.GetDOFHolderDataHost().size() !=
      detail::ExpectedFlatFieldDataSize(layout)) {
    throw pcms_error(
      "LagrangeFunctionSpace::CreateField: field data size does not match "
      "layout");
  }
}

} // namespace

LagrangeFunctionSpace::LagrangeFunctionSpace(
  Key, std::shared_ptr<const FieldLayout> layout,
  std::function<FieldDataVariant(Type, FieldMetadata)> create_field_data_fn,
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept
  : layout_(std::move(layout)),
    create_field_data_fn_(std::move(create_field_data_fn)),
    evaluator_factory_(std::move(evaluator_factory))
{
}

std::shared_ptr<LagrangeFunctionSpace> LagrangeFunctionSpace::FromMesh(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system, std::string global_id_name,
  Backend backend, std::string layout_name)
{
  // https://github.com/SCOREC/meshFields/issues/88
  if (backend == Backend::MeshFields && order == 0) {
    // MeshFields does not support order-0 fields; fall back to OmegaH.
    backend = Backend::OmegaH;
  }
  if (backend == Backend::MeshFields) {
#ifdef PCMS_ENABLE_MESHFIELDS
    std::array<int, 4> nodes_per_dim{};
    switch (order) {
      case 1: nodes_per_dim = {1, 0, 0, 0}; break;
      case 2: nodes_per_dim = {1, 1, 0, 0}; break;
      default: throw pcms_error("Unimplemented Lagrange order");
    }
    auto mesh_layout = std::make_shared<MeshFieldsAdapterLayout>(
      mesh, nodes_per_dim, num_components, coordinate_system,
      std::move(global_id_name));
    mesh_layout->SetName(layout_name);
    auto eval_factory =
      std::make_shared<MeshFieldsEvaluatorFactory<Real>>(mesh_layout);
    return std::make_shared<LagrangeFunctionSpace>(
      Key{}, mesh_layout,
      [mesh_layout](Type t, FieldMetadata metadata) -> FieldDataVariant {
        if (t == Type::Float) {
          if constexpr (std::is_same_v<MeshField::Real4, float> ||
                        std::is_same_v<MeshField::Real8, float>) {
            return std::make_unique<MeshFieldsFieldData<float>>(mesh_layout,
                                                                metadata);
          }
          throw pcms_error(
            "LagrangeFunctionSpace: MeshFields backend does not support "
            "float in this build");
        }
        if (t == Type::Real) {
          return std::make_unique<MeshFieldsFieldData<double>>(mesh_layout,
                                                               metadata);
        }
        throw pcms_error(
          "LagrangeFunctionSpace: MeshFields backend only supports the "
          "MeshFields scalar types enabled in this build");
      },
      std::move(eval_factory));
#else
    throw pcms_error(
      "LagrangeFunctionSpace::FromMesh: MeshFields backend requested but "
      "PCMS_ENABLE_MESHFIELDS is not set");
#endif
  }
  if (order > 1) {
    throw pcms_error(
      "LagrangeFunctionSpace::FromMesh: OmegaH backend only supports "
      "Lagrange orders 0 and 1");
  }
  auto layout = std::make_shared<OmegaHLagrangeLayout>(
    mesh, order, num_components, coordinate_system, std::move(global_id_name));
  layout->SetName(std::move(layout_name));
  auto eval_factory =
    std::make_shared<OmegaHLagrangeEvaluatorFactory<Real>>(layout);
  return std::make_shared<LagrangeFunctionSpace>(
    Key{}, layout,
    [layout](Type t, FieldMetadata metadata) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("LagrangeFunctionSpace: int8_t is not supported");
        } else {
          return std::make_unique<SimpleFieldData<T>>(layout, metadata);
        }
      });
    },
    std::move(eval_factory));
}

std::shared_ptr<LagrangeFunctionSpace> LagrangeFunctionSpace::FromMesh(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system, Omega_h::Read<Omega_h::I8> owned_mask,
  std::string global_id_name, Backend backend, std::string layout_name)
{
  if (backend == Backend::MeshFields) {
    throw pcms_error(
      "LagrangeFunctionSpace::FromMesh: owned-mask construction is only "
      "supported by the OmegaH backend");
  }
  if (order > 1) {
    throw pcms_error(
      "LagrangeFunctionSpace::FromMesh: OmegaH backend only supports "
      "Lagrange orders 0 and 1");
  }
  auto layout = std::make_shared<OmegaHLagrangeLayout>(
    mesh, order, num_components, coordinate_system, std::move(owned_mask),
    std::move(global_id_name));
  layout->SetName(std::move(layout_name));
  auto eval_factory =
    std::make_shared<OmegaHLagrangeEvaluatorFactory<Real>>(layout);
  return std::make_shared<LagrangeFunctionSpace>(
    Key{}, layout,
    [layout](Type t, FieldMetadata metadata) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("LagrangeFunctionSpace: int8_t is not supported");
        } else {
          return std::make_unique<SimpleFieldData<T>>(layout, metadata);
        }
      });
    },
    std::move(eval_factory));
}

std::shared_ptr<LagrangeFunctionSpace> LagrangeFunctionSpace::FromUniformGrid(
  const UniformGrid<2>& grid, int num_components,
  CoordinateSystem coordinate_system, int order, std::string layout_name)
{
  if (order != 0 && order != 1) {
    throw std::invalid_argument("LagrangeFunctionSpace::FromUniformGrid: only "
                                "orders 0 and 1 are supported");
  }
  auto ug_layout = std::make_shared<UniformGridFieldLayout<2>>(
    grid, num_components, coordinate_system, order);
  ug_layout->SetName(std::move(layout_name));
  auto eval_factory =
    std::make_shared<UniformGridEvaluatorFactory<2>>(ug_layout);
  return std::make_shared<LagrangeFunctionSpace>(
    Key{}, ug_layout,
    [ug_layout](Type t, FieldMetadata metadata) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("LagrangeFunctionSpace: int8_t is not supported");
        } else {
          return std::make_unique<SimpleFieldData<T>>(ug_layout, metadata);
        }
      });
    },
    std::move(eval_factory));
}

std::shared_ptr<LagrangeFunctionSpace> LagrangeFunctionSpace::FromUniformGrid(
  const UniformGrid<3>& grid, int num_components,
  CoordinateSystem coordinate_system, int order, std::string layout_name)
{
  if (order != 0 && order != 1) {
    throw std::invalid_argument("LagrangeFunctionSpace::FromUniformGrid: only "
                                "orders 0 and 1 are supported");
  }
  auto ug_layout = std::make_shared<UniformGridFieldLayout<3>>(
    grid, num_components, coordinate_system, order);
  ug_layout->SetName(std::move(layout_name));
  auto eval_factory =
    std::make_shared<UniformGridEvaluatorFactory<3>>(ug_layout);
  return std::make_shared<LagrangeFunctionSpace>(
    Key{}, ug_layout,
    [ug_layout](Type t, FieldMetadata metadata) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("LagrangeFunctionSpace: int8_t is not supported");
        } else {
          return std::make_unique<SimpleFieldData<T>>(ug_layout, metadata);
        }
      });
    },
    std::move(eval_factory));
}

std::shared_ptr<const FieldLayout> LagrangeFunctionSpace::GetLayout()
  const noexcept
{
  return layout_;
}

CoordinateSystem LagrangeFunctionSpace::GetCoordinateSystem() const noexcept
{
  return evaluator_factory_->GetCoordinateSystem();
}

FieldVariant LagrangeFunctionSpace::CreateFieldImpl(
  Type value_type, FieldMetadata metadata) const
{
  auto field_data = create_field_data_fn_(value_type, metadata);
  return std::visit(
    [this](auto&& fd) -> FieldVariant {
      using FD = std::decay_t<decltype(fd)>;
      using T = typename FD::element_type::value_type;
      return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
    },
    std::move(field_data));
}

FieldVariant LagrangeFunctionSpace::CreateFieldImpl(FieldDataVariant data) const
{
  return std::visit(
    [this](auto&& fd) -> FieldVariant {
      using FD = std::decay_t<decltype(fd)>;
      using T = typename FD::element_type::value_type;
      PCMS_ALWAYS_ASSERT(fd != nullptr);
      if constexpr (std::is_same_v<T, int8_t>) {
        throw pcms_error(
          "LagrangeFunctionSpace: int8_t is not a supported field type");
      } else {
        ValidateLagrangeWrappedFieldData(*layout_, *fd);
        return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
      }
    },
    std::move(data));
}

PointEvaluatorVariant LagrangeFunctionSpace::CreatePointEvaluatorImpl(
  Type value_type, const EvaluationRequest& request) const
{
  if (value_type != Type::Real) {
    throw pcms_error(
      "LagrangeFunctionSpace: point evaluation only supports double (Real)");
  }
  if (!evaluator_factory_) {
    throw pcms_error(
      "LagrangeFunctionSpace::CreatePointEvaluatorImpl: evaluator construction "
      "is not available for this backend");
  }
  return evaluator_factory_->CreatePointEvaluator(request);
}

} // namespace pcms
