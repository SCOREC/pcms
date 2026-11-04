#ifndef PCMS_FIELD_FUNCTION_SPACES_HPP
#define PCMS_FIELD_FUNCTION_SPACES_HPP

// Factories for the function spaces PCMS provides. Each takes a
// Discretization (see pcms/discretization/discretization.hpp and the backend
// headers) and returns the FunctionSpace interface.

#include "pcms/discretization/discretization.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.hpp"

#include <memory>
#include <string>

namespace pcms
{

/// Which implementation evaluates Lagrange fields on an Omega_h
/// discretization. Default is MeshFields when PCMS is built with MeshFields,
/// Omega_h otherwise.
enum class LagrangeImplementation
{
  Default,
  MeshFields,
  OmegaH
};

struct LagrangeOptions
{
  LagrangeImplementation implementation = LagrangeImplementation::Default;
  /// Name of the layout; it identifies the layout in the coupler.
  std::string layout_name;
};

/// Nodal Lagrange space of the given order with num_components components on
/// an Omega_h discretization (orders 0-2, depending on the implementation)
/// or a 2D/3D uniform-grid discretization (order 1). Throws pcms_error for
/// other discretizations.
std::shared_ptr<FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Discretization> discretization, int order,
  int num_components, std::shared_ptr<const CoordinateSystem> coordinate_system,
  LagrangeOptions options = {});

/// Bicubic spline space on a 2D uniform-grid discretization.
std::shared_ptr<FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Discretization> discretization,
  std::shared_ptr<const CoordinateSystem> coordinate_system);

/// Moving-least-squares polynomial reconstruction from values at points.
std::shared_ptr<FunctionSpace> CreatePolynomialReconstructionSpace(
  CoordinateView<HostMemorySpace> points, MLSOptions options = {});

/// As above, with the points being the entities of dimension
/// source_entity_dim of an Omega_h discretization.
std::shared_ptr<FunctionSpace> CreatePolynomialReconstructionSpace(
  std::shared_ptr<const Discretization> discretization, int source_entity_dim,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  MLSOptions options = {});

/// Layout without degrees of freedom, for evaluation-only spaces that
/// implement FunctionSpace outside PCMS (e.g. fields of an external code).
std::shared_ptr<const FieldLayout> CreateEmptyLayout(
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  int num_components = 1);

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACES_HPP
