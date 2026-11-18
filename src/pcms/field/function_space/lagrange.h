#ifndef PCMS_LAGRANGE_FIELD_FACTORY_H
#define PCMS_LAGRANGE_FIELD_FACTORY_H

#include "pcms/configuration.h"
#include "pcms/field/function_space.h"
#include "pcms/geometry/geometry.hpp"

#include <memory>
#include <string>

namespace pcms
{

enum class LagrangeBackend
{
  MeshFields,
  OmegaH
};

#ifdef PCMS_ENABLE_MESHFIELDS
inline constexpr LagrangeBackend DefaultLagrangeBackend =
  LagrangeBackend::MeshFields;
#else
inline constexpr LagrangeBackend DefaultLagrangeBackend =
  LagrangeBackend::OmegaH;
#endif

struct LagrangeOptions
{
  LagrangeBackend backend = DefaultLagrangeBackend;
  /// name of the space's layout; "" is the default
  std::string layout_name;
};

/// A Lagrange space of order `order` with `num_components` components on
/// `geometry`: an Omega_h geometry (backend OmegaH: orders 0 and 1; MeshFields:
/// orders 1 and 2) or a uniform-grid geometry (orders 0 and 1, backend
/// ignored). Spaces on one geometry accept each other's bound sources.
/// @throws pcms_error for another geometry kind, a null geometry, or an order
///         the backend does not support
std::shared_ptr<FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Geometry> geometry, int order, int num_components,
  LagrangeOptions options = {});

} // namespace pcms

#endif // PCMS_LAGRANGE_FIELD_FACTORY_H
