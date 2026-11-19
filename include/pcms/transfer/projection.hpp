#ifndef PCMS_TRANSFER_PROJECTION_HPP
#define PCMS_TRANSFER_PROJECTION_HPP

// Factories for the Galerkin projections. They are always declared; in a
// build without the implementing backends (PETSc and MeshFields) they throw
// pcms_error.

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"

#include <cstdint>
#include <memory>

namespace pcms
{

struct ConservativeProjectionOptions
{
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
};

/// Conservative (Galerkin L2) projection from source_space to target_space
/// using exact mesh intersection. Both spaces must be scalar Lagrange spaces
/// on Omega_h discretizations. The operator holds both spaces.
/// @throws pcms_error if a space is null
std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(
  std::shared_ptr<const FunctionSpace> source_space,
  std::shared_ptr<const FunctionSpace> target_space,
  ConservativeProjectionOptions options = {});

struct ControlVariateProjectionOptions
{
  int samples_per_element = 1;
  std::uint64_t seed = 8675309;
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
  /// map from the target's coordinate system to the source's, or null if they
  /// are the same; the source is evaluated at the mapped points
  std::shared_ptr<const CoordinateMap> to_source;
};

/// Galerkin projection of any evaluable source onto target_space, with the
/// right-hand side estimated by Monte Carlo sampling with a control variate.
/// target_space must be a scalar Lagrange space on an Omega_h
/// discretization. The operator holds `source` and `target_space`.
/// @throws pcms_error if `source` or `target_space` is null
std::unique_ptr<TransferOperator<Real>> CreateControlVariateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space,
  ControlVariateProjectionOptions options = {});

} // namespace pcms

#endif // PCMS_TRANSFER_PROJECTION_HPP
