#include "pcms/field/function_space/lagrange.h"
#include "pcms/configuration.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"
#ifdef PCMS_ENABLE_MESHFIELDS
#include "pcms/field/layout/mesh_fields.h"
#include "pcms/field/evaluator/mesh_fields.h"
#include "pcms/field/data/mesh_fields.h"
#endif
#include "pcms/field/evaluator/omega_h_geometry.hpp"
#include "pcms/field/evaluator/uniform_grid_geometry.hpp"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/evaluator/omega_h_lagrange.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/field/evaluator/uniform_grid.h"
#include "pcms/field/data/omega_h.h"
#include "pcms/field/data/simple.h"
#include "pcms/utility/uniform_grid.h"

#include <array>
#include <stdexcept>
#include <string>

namespace pcms
{

namespace
{

/// Checks that `backend` supports Lagrange order `order` (and an owned mask
/// if `owned_mask`) before anything is built.
/// @return the backend to use: MeshFields falls back to OmegaH for order 0
LagrangeBackend CheckOrderAndBackend(int order, LagrangeBackend backend,
                                     bool owned_mask)
{
  // https://github.com/SCOREC/meshFields/issues/88
  if (backend == LagrangeBackend::MeshFields && order == 0 && !owned_mask) {
    // MeshFields does not support order-0 fields; fall back to OmegaH.
    backend = LagrangeBackend::OmegaH;
  }
  if (backend == LagrangeBackend::MeshFields) {
    if (owned_mask) {
      throw pcms_error("CreateLagrangeSpace: the geometry's owned mask is "
                       "only supported by the OmegaH backend");
    }
#ifdef PCMS_ENABLE_MESHFIELDS
    if (order != 1 && order != 2) {
      throw pcms_error("CreateLagrangeSpace: the MeshFields backend only "
                       "supports Lagrange orders 1 and 2");
    }
#else
    throw pcms_error("CreateLagrangeSpace: MeshFields backend requested but "
                     "PCMS_ENABLE_MESHFIELDS is not set");
#endif
    return backend;
  }
  if (order > 1) {
    throw pcms_error("CreateLagrangeSpace: the OmegaH backend only supports "
                     "Lagrange orders 0 and 1");
  }
  return backend;
}

/// A Lagrange space on the mesh of an Omega_h geometry; dispatches to
/// MeshFields or the native Omega_h backend.
std::shared_ptr<FunctionSpace> CreateOmegaHLagrangeSpace(
  std::shared_ptr<const OmegaHGeometry> geometry, int order, int num_components,
  const LagrangeOptions& options)
{
  const auto backend =
    CheckOrderAndBackend(order, options.backend,
                         /*owned_mask=*/geometry->GetOwnedMask().size() > 0);
  if (backend == LagrangeBackend::MeshFields) {
#ifdef PCMS_ENABLE_MESHFIELDS
    const std::array<int, 4> nodes_per_dim = order == 1
                                               ? std::array<int, 4>{1, 0, 0, 0}
                                               : std::array<int, 4>{1, 1, 0, 0};
    auto mesh_layout = std::make_shared<MeshFieldsAdapterLayout>(
      geometry, nodes_per_dim, num_components);
    mesh_layout->SetName(options.layout_name);
    auto eval_factory =
      std::make_shared<MeshFieldsEvaluatorFactory<Real>>(geometry, mesh_layout);
    return std::make_shared<FunctionSpace>(
      mesh_layout,
      [mesh_layout](Type t, ValueBasis basis) -> FieldDataVariant {
        if (t == Type::Float) {
          if constexpr (std::is_same_v<MeshField::Real4, float> ||
                        std::is_same_v<MeshField::Real8, float>) {
            return std::make_unique<MeshFieldsFieldData<float>>(
              mesh_layout, std::move(basis));
          }
          throw pcms_error("CreateLagrangeSpace: MeshFields backend does not "
                           "support float in this build");
        }
        if (t == Type::Real) {
          return std::make_unique<MeshFieldsFieldData<double>>(
            mesh_layout, std::move(basis));
        }
        throw pcms_error("CreateLagrangeSpace: MeshFields backend only "
                         "supports the MeshFields scalar types enabled in "
                         "this build");
      },
      std::move(eval_factory));
#endif
  }
  auto layout =
    std::make_shared<OmegaHLagrangeLayout>(geometry, order, num_components);
  layout->SetName(options.layout_name);
  auto eval_factory = std::make_shared<OmegaHLagrangeEvaluatorFactory<Real>>(
    std::move(geometry), layout);
  return std::make_shared<FunctionSpace>(
    layout,
    [layout](Type t, ValueBasis basis) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("CreateLagrangeSpace: int8_t is not supported");
        } else {
          return std::make_unique<OmegaHFieldData<T>>(layout, basis);
        }
      });
    },
    std::move(eval_factory));
}

/// A Lagrange space of order 0 or 1 on the grid of a uniform-grid geometry.
template <unsigned Dim>
std::shared_ptr<FunctionSpace> CreateUniformGridLagrangeSpace(
  std::shared_ptr<const UniformGridGeometry<Dim>> geometry, int order,
  int num_components, const LagrangeOptions& options)
{
  if (order != 0 && order != 1) {
    throw pcms_error("CreateLagrangeSpace: uniform-grid geometries only "
                     "support Lagrange orders 0 and 1");
  }
  auto ug_layout = std::make_shared<UniformGridFieldLayout<Dim>>(
    geometry, num_components, order);
  ug_layout->SetName(options.layout_name);
  auto eval_factory = std::make_shared<UniformGridEvaluatorFactory<Dim>>(
    std::move(geometry), ug_layout);
  return std::make_shared<FunctionSpace>(
    ug_layout,
    [ug_layout](Type t, ValueBasis basis) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (std::is_same_v<T, int8_t>) {
          throw pcms_error("CreateLagrangeSpace: int8_t is not supported");
        } else {
          return std::make_unique<SimpleFieldData<T>>(ug_layout, basis);
        }
      });
    },
    std::move(eval_factory));
}

} // namespace

std::shared_ptr<FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Geometry> geometry, int order, int num_components,
  LagrangeOptions options)
{
  if (geometry == nullptr) {
    throw pcms_error("CreateLagrangeSpace: the geometry must not be null");
  }
  if (auto omega_h =
        std::dynamic_pointer_cast<const OmegaHGeometry>(geometry)) {
    return CreateOmegaHLagrangeSpace(std::move(omega_h), order, num_components,
                                     options);
  }
  if (auto grid =
        std::dynamic_pointer_cast<const UniformGridGeometry<2>>(geometry)) {
    return CreateUniformGridLagrangeSpace<2>(std::move(grid), order,
                                             num_components, options);
  }
  if (auto grid =
        std::dynamic_pointer_cast<const UniformGridGeometry<3>>(geometry)) {
    return CreateUniformGridLagrangeSpace<3>(std::move(grid), order,
                                             num_components, options);
  }
  throw pcms_error("CreateLagrangeSpace: only Omega_h and uniform-grid "
                   "geometries are supported");
}

} // namespace pcms
