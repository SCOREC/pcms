#ifndef PCMS_TRANSFER_PROJECTION_HPP
#define PCMS_TRANSFER_PROJECTION_HPP

// Factories for the Galerkin projections. They are always declared; in a
// build without the implementing backends (PETSc and MeshFields) they throw
// pcms_error.

#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/monte_carlo_sampling.hpp"
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
/// on Omega_h discretizations.
std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  ConservativeProjectionOptions options = {});

struct ControlVariateProjectionOptions
{
  int samples_per_element = 1;
  MonteCarloSampling sampling = MonteCarloSampling::UniformRandom;
  std::uint64_t seed = 8675309;
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
};

/// Galerkin projection of any evaluable source onto target_space, with the
/// right-hand side estimated by Monte Carlo sampling with a control variate.
/// target_space must be a scalar Lagrange space on an Omega_h
/// discretization.
std::unique_ptr<TransferOperator<Real>> CreateControlVariateProjection(
  const PointEvaluatorFactory& source, const FunctionSpace& target_space,
  ControlVariateProjectionOptions options = {});

} // namespace pcms

#endif // PCMS_TRANSFER_PROJECTION_HPP
