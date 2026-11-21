#include <algorithm>
#include "pcms/field/function_spaces.hpp"
#include "pcms/configuration.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/common.hpp"
#ifdef PCMS_ENABLE_OMEGA_H
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
#include "pcms/backends/omega_h_meshfields/layout.hpp"
#include "pcms/backends/omega_h_meshfields/evaluator.hpp"
#include "pcms/backends/omega_h_meshfields/field_data.hpp"
#endif
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/backends/omega_h/lagrange_evaluator.hpp"
#include "pcms/backends/omega_h/field_data.hpp"
#endif
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/backends/uniform_grid/evaluator.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace pcms
{
#ifdef PCMS_ENABLE_OMEGA_H
using omega_h::OmegaHFieldData;
using omega_h::OmegaHGeometry;
using omega_h::OmegaHLagrangeEvaluatorFactory;
using omega_h::OmegaHLagrangeLayout;
#endif
using uniform_grid::UniformGridEvaluatorFactory;
using uniform_grid::UniformGridFieldLayout;
using uniform_grid::UniformGridGeometry;
using uniform_grid::VisitUniformGrid;
} // namespace pcms

namespace pcms
{

namespace
{

/// @brief A Lagrange space of order p; its quadrature defaults to the
///        geometry's rule of degree max(1, 2p).
class LagrangeFunctionSpace : public FunctionSpace
{
public:
  LagrangeFunctionSpace(std::shared_ptr<const FieldLayout> layout,
                        std::shared_ptr<const PointEvaluatorFactory> evaluator,
                        std::shared_ptr<const Geometry> geometry, int order)
    : FunctionSpace(std::move(layout), std::move(evaluator),
                    std::move(geometry)),
      order_(order)
  {
  }

  [[nodiscard]] Quadrature CreateQuadrature() const override
  {
    return CreateGaussQuadrature(GetGeometry(), std::max(1, 2 * order_));
  }

private:
  int order_;
};

#ifdef PCMS_ENABLE_OMEGA_H
/// Throws unless every alternative of `data` is a `Data<T>`.
template <template <typename> class Data>
void CheckFieldDataType(const FieldDataVariant& data, const char* message)
{
  std::visit(
    [&](const auto& fd) {
      using T = typename std::decay_t<decltype(fd)>::element_type::value_type;
      if (dynamic_cast<const Data<T>*>(fd.get()) == nullptr) {
        throw pcms_error(message);
      }
    },
    data);
}

#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
/// The MeshFields flavour of the space: its fields are MeshFieldsFieldData.
class MeshFieldsLagrangeFunctionSpace final : public LagrangeFunctionSpace
{
public:
  MeshFieldsLagrangeFunctionSpace(
    std::shared_ptr<const omega_h_meshfields::MeshFieldsAdapterLayout> layout,
    std::shared_ptr<const PointEvaluatorFactory> evaluator,
    std::shared_ptr<const Geometry> geometry, int order)
    : LagrangeFunctionSpace(layout, std::move(evaluator), std::move(geometry),
                            order),
      layout_(std::move(layout))
  {
  }

protected:
  FieldDataVariant CreateFieldData(Type t, ValueBasis basis) const override
  {
    if (t == Type::Float) {
      if constexpr (std::is_same_v<MeshField::Real4, float> ||
                    std::is_same_v<MeshField::Real8, float>) {
        return std::make_unique<omega_h_meshfields::MeshFieldsFieldData<float>>(
          layout_, std::move(basis));
      }
      throw pcms_error("CreateLagrangeSpace: MeshFields backend does not "
                       "support float in this build");
    }
    if (t == Type::Real) {
      return std::make_unique<omega_h_meshfields::MeshFieldsFieldData<double>>(
        layout_, std::move(basis));
    }
    throw pcms_error("CreateLagrangeSpace: MeshFields backend only "
                     "supports the MeshFields scalar types enabled in "
                     "this build");
  }

  void ValidateFieldData(const FieldDataVariant& data) const override
  {
    CheckFieldDataType<omega_h_meshfields::MeshFieldsFieldData>(
      data, "FunctionSpace::CreateField: MeshFields layout requires "
            "MeshFieldsFieldData");
    FunctionSpace::ValidateFieldData(data);
  }

private:
  std::shared_ptr<const omega_h_meshfields::MeshFieldsAdapterLayout> layout_;
};
#endif

/// The native flavour of the space: its fields are OmegaHFieldData.
class NativeOmegaHLagrangeFunctionSpace final : public LagrangeFunctionSpace
{
public:
  NativeOmegaHLagrangeFunctionSpace(
    std::shared_ptr<const OmegaHLagrangeLayout> layout,
    std::shared_ptr<const PointEvaluatorFactory> evaluator,
    std::shared_ptr<const Geometry> geometry, int order)
    : LagrangeFunctionSpace(layout, std::move(evaluator), std::move(geometry),
                            order),
      layout_(std::move(layout))
  {
  }

protected:
  FieldDataVariant CreateFieldData(Type t, ValueBasis basis) const override
  {
    return detail::apply_to_type(t, [&](auto tag) -> FieldDataVariant {
      using T = typename decltype(tag)::type;
      if constexpr (std::is_same_v<T, int8_t>) {
        throw pcms_error("CreateLagrangeSpace: int8_t is not supported");
      } else {
        return std::make_unique<OmegaHFieldData<T>>(layout_, basis);
      }
    });
  }

  void ValidateFieldData(const FieldDataVariant& data) const override
  {
    CheckFieldDataType<OmegaHFieldData>(
      data, "FunctionSpace::CreateField: OmegaH layout requires "
            "OmegaHFieldData");
    FunctionSpace::ValidateFieldData(data);
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
};

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
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
    if (order != 1 && order != 2) {
      throw pcms_error("CreateLagrangeSpace: the MeshFields backend only "
                       "supports Lagrange orders 1 and 2");
    }
#else
    throw pcms_error("CreateLagrangeSpace: MeshFields backend requested but "
                     "PCMS_ENABLE_OMEGA_H_MESHFIELDS is not set");
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
std::shared_ptr<const FunctionSpace> CreateOmegaHLagrangeSpace(
  std::shared_ptr<const OmegaHGeometry> geometry, int order, int num_components,
  const LagrangeOptions& options)
{
  const auto backend =
    CheckOrderAndBackend(order, options.backend,
                         /*owned_mask=*/geometry->GetOwnedMask().size() > 0);
  if (backend == LagrangeBackend::MeshFields) {
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
    const std::array<int, 4> nodes_per_dim = order == 1
                                               ? std::array<int, 4>{1, 0, 0, 0}
                                               : std::array<int, 4>{1, 1, 0, 0};
    auto mesh_layout =
      std::make_shared<omega_h_meshfields::MeshFieldsAdapterLayout>(
        geometry, nodes_per_dim, num_components, options.layout_name);
    auto eval_factory =
      std::make_shared<omega_h_meshfields::MeshFieldsEvaluatorFactory<Real>>(
        geometry, mesh_layout);
    return std::make_shared<MeshFieldsLagrangeFunctionSpace>(
      mesh_layout, eval_factory, eval_factory->GetGeometry(), order);
#endif
  }
  auto layout = std::make_shared<OmegaHLagrangeLayout>(
    geometry, order, num_components, options.layout_name);
  auto eval_factory = std::make_shared<OmegaHLagrangeEvaluatorFactory<Real>>(
    std::move(geometry), layout);
  return std::make_shared<NativeOmegaHLagrangeFunctionSpace>(
    layout, eval_factory, eval_factory->GetGeometry(), order);
}

#endif // PCMS_ENABLE_OMEGA_H

/// A Lagrange space of order 0 or 1 on the grid of a uniform-grid geometry.
template <unsigned Dim>
std::shared_ptr<const FunctionSpace> CreateUniformGridLagrangeSpace(
  std::shared_ptr<const UniformGridGeometry<Dim>> geometry, int order,
  int num_components, const LagrangeOptions& options)
{
  if (order != 0 && order != 1) {
    throw pcms_error("CreateLagrangeSpace: uniform-grid geometries only "
                     "support Lagrange orders 0 and 1");
  }
  auto ug_layout = std::make_shared<UniformGridFieldLayout<Dim>>(
    geometry, num_components, order, options.layout_name);
  auto eval_factory = std::make_shared<UniformGridEvaluatorFactory<Dim>>(
    std::move(geometry), ug_layout);
  return std::make_shared<LagrangeFunctionSpace>(
    ug_layout, eval_factory, eval_factory->GetGeometry(), order);
}

} // namespace

std::shared_ptr<const FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Geometry> geometry, int order, int num_components,
  LagrangeOptions options)
{
  if (geometry == nullptr) {
    throw pcms_error("CreateLagrangeSpace: the geometry must not be null");
  }
#ifdef PCMS_ENABLE_OMEGA_H
  if (auto omega_h =
        std::dynamic_pointer_cast<const OmegaHGeometry>(geometry)) {
    return CreateOmegaHLagrangeSpace(std::move(omega_h), order, num_components,
                                     options);
  }
#endif
  if (auto space = VisitUniformGrid(geometry, [&](auto grid) {
        return CreateUniformGridLagrangeSpace(std::move(grid), order,
                                              num_components, options);
      })) {
    return space;
  }
  throw pcms_error("CreateLagrangeSpace: only Omega_h and uniform-grid "
                   "geometries are supported");
}

} // namespace pcms
