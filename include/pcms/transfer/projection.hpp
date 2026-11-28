#ifndef PCMS_TRANSFER_PROJECTION_HPP
#define PCMS_TRANSFER_PROJECTION_HPP

// Factories for the Galerkin projections. They are always built;
// CreateConservativeProjection throws pcms_error in a build without
// PCMS_ENABLE_OMEGA_H_MESHFIELDS.

#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"

#include <cstdint>
#include <memory>

namespace pcms
{

/// Galerkin L2 projection of any evaluable source onto target_space: solves
/// the mass system M x = b with b = E^T W f, where E evaluates the target
/// space at the points of a quadrature, W holds its weights and f is the
/// source at the points. The right-hand side uses a Gauss rule on the target
/// space's geometry exact for its natural quadrature's degree (at least 2),
/// the mass matrix the target space's natural quadrature. A field with several
/// components is projected jointly (the components do not couple). The source
/// is evaluated at the target's points in the target's coordinate system: a
/// source in another frame needs a target geometry embedded in that frame.
/// The operator holds `source` and `target_space`.
/// @throws pcms_error if `source` or `target_space` is null, the source has no
///         geometry, or the target space's evaluator does not support the
///         transpose or has no natural quadrature
std::unique_ptr<TransferOperator<Real>> CreateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space);

struct ConservativeProjectionOptions
{
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
};

/// Conservative (Galerkin L2) projection from source_space to target_space
/// using exact mesh intersection. Both spaces must be scalar Lagrange spaces
/// on Omega_h discretizations. The operator holds both spaces.
/// @throws pcms_error if a space is null, or in a build without
///         PCMS_ENABLE_OMEGA_H_MESHFIELDS
std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(
  std::shared_ptr<const FunctionSpace> source_space,
  std::shared_ptr<const FunctionSpace> target_space,
  ConservativeProjectionOptions options = {});

struct ControlVariateProjectionOptions
{
  int samples_per_element = 1;
  std::uint64_t seed = DefaultQuadratureSeed;
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
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
