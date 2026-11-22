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
#include <vector>

namespace pcms
{

enum class LagrangeBackend
{
  MeshFields,
  OmegaH
};

#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
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
/// ignored; an axis identified over the grid's extent has no last node).
/// Spaces on one geometry accept each other's BoundPoints.
/// @throws pcms_error for another geometry kind, a null geometry, or an order
///         the backend does not support
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Geometry> geometry, int order, int num_components,
  LagrangeOptions options = {});

/// @brief A cubic (1D) or bicubic (2D) spline space on a uniform-grid
///        geometry; each component is splined on its own. On a 1D grid
///        identified over its extent, the periodic cubic spline of the node
///        values (with EvaluateTranspose and GetMatrix).
/// @throws pcms_error if `geometry` is null or not a 1D or 2D uniform-grid
///         geometry, has a 2D grid with an identified axis, or has fewer than
///         3 cells on an identified axis
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Geometry> geometry, int num_components = 1);

/// @name 1D spaces on a uniform grid identified over its extent
/// Real, with Value and Gradient evaluation, EvaluateTranspose and
/// GetMatrix. Hermite holder 2 i is the value and 2 i + 1 the derivative of
/// node i; Fourier holders have no positions.
/// @throws pcms_error if `geometry` is not such a grid
/// @{

/// Cubic Hermite: the value and derivative of each node.
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreatePeriodicHermiteSpace(
  std::shared_ptr<const Geometry> geometry, int num_components = 1);

/// Real Fourier series: for each m of `modes` the DOFs cos(m stride theta)
/// and, for m > 0, sin(m stride theta), in the order of the list, theta =
/// 2 pi (x - lower) / length; `stride` is nfp for fields periodic in a field
/// period.
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateFourierSpace(
  std::shared_ptr<const Geometry> geometry, std::vector<int> modes,
  int stride = 1, int num_components = 1);
/// @}

/// The product of two scalar spaces: the basis phi_i(a) psi_k(b), DOF (i, k)
/// with global id gid_a * (global DOFs of b) + gid_b.
/// Evaluation (Value, and Gradient when both factors have it),
/// EvaluateTranspose and GetMatrix are the Kronecker products of the factors'
/// rows, so both factors' evaluators must provide GetMatrix (for Gradient also
/// of the partials); the natural quadrature is the tensor product of the
/// factors'. Serial: the factors' layouts are on one rank.
/// @param geometry the product geometry to create the space on, or null for
///        CreateProductGeometry(a's geometry, b's geometry); to use the one
///        made with a twist-shift identification, or to share one between
///        spaces
/// @throws pcms_error if a factor is null or not scalar, `geometry` is not the
///         product geometry of the factors' geometries, or a factor's layout
///         is distributed
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateProductSpace(
  std::shared_ptr<const FunctionSpace> a,
  std::shared_ptr<const FunctionSpace> b,
  std::shared_ptr<const Geometry> geometry = nullptr);

#ifdef PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION
/// An MLS reconstruction space on a point-cloud geometry (its sites are the
/// DOF holders); `options` are the reconstruction parameters, the support
/// search is the geometry's. Spaces on one geometry accept each other's bound
/// sources.
/// @throws pcms_error if `geometry` is null or not a point-cloud geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace>
CreatePolynomialReconstructionSpace(std::shared_ptr<const Geometry> geometry,
                                    MLSOptions options = {});
#endif

/// Layout without degrees of freedom, for evaluation-only spaces that
/// implement FunctionSpace outside PCMS (e.g. fields of an external code).
std::shared_ptr<const FieldLayout> CreateEmptyLayout(
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  int num_components = 1);

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACES_HPP
