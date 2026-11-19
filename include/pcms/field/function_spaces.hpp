#ifndef PCMS_FIELD_FUNCTION_SPACES_HPP
#define PCMS_FIELD_FUNCTION_SPACES_HPP

// Factories for the function spaces PCMS provides. Each takes the Geometry the
// space is defined on (uniform_grid::CreateGeometry,
// point_cloud::CreateGeometry, omega_h::CreateGeometry) and returns the
// FunctionSpace.

#include "pcms/configuration.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/utility/coordinate_system.hpp"

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
/// ignored). Spaces on one geometry accept each other's BoundPoints.
/// @throws pcms_error for another geometry kind, a null geometry, or an order
///         the backend does not support
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Geometry> geometry, int order, int num_components,
  LagrangeOptions options = {});

/// A bicubic spline space on a 2D uniform-grid geometry. Spaces on one
/// geometry accept each other's BoundPoints.
/// @throws pcms_error if `geometry` is null or not a 2D uniform-grid geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Geometry> geometry);

/// An MLS reconstruction space on a point-cloud geometry (its sites are the
/// DOF holders); `options` are the reconstruction parameters, the support
/// search is the geometry's. Spaces on one geometry accept each other's bound
/// sources.
/// @throws pcms_error if `geometry` is null or not a point-cloud geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace>
CreatePolynomialReconstructionSpace(std::shared_ptr<const Geometry> geometry,
                                    MLSOptions options = {});

/// Layout without degrees of freedom, for evaluation-only spaces that
/// implement FunctionSpace outside PCMS (e.g. fields of an external code).
std::shared_ptr<const FieldLayout> CreateEmptyLayout(
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  int num_components = 1);

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACES_HPP
