#ifndef PCMS_TRANSFER_PROJECTION_HPP
#define PCMS_TRANSFER_PROJECTION_HPP

// Factories for the Galerkin projections. They are always built;
// CreateConservativeProjection throws pcms_error in a build without
// PCMS_ENABLE_OMEGA_H_MESHFIELDS, and the Assembled solve in one without
// PCMS_ENABLE_PETSC.

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

namespace pcms
{

/// How the mass system of a projection is solved.
struct LinearSolveOptions
{
  enum class Method
  {
    /// Jacobi-preconditioned conjugate gradients on all components jointly,
    /// written on the evaluation and the vector operations of the kinds: no
    /// matrix, no host round trips
    MatrixFree,
    /// the matrix E^T W E is assembled from the evaluation matrix into PETSc
    /// (serial) and solved by a KSP: a direct factorization unless the
    /// options of `petsc_options_prefix` say otherwise. Needs a build with
    /// PCMS_ENABLE_PETSC
    Assembled
  };
  Method method = Method::MatrixFree;
  /// MatrixFree: relative residual to stop at
  Real relative_tolerance = 1e-12;
  /// MatrixFree: iteration limit; 0 is the number of unknowns
  int max_iterations = 0;
  /// Assembled: the prefix of the PETSc options of the KSP and its PC
  std::string petsc_options_prefix;
};

struct ProjectionOptions
{
  /// the quadrature of the right-hand side, on the target space's geometry
  /// (e.g. CreateGaussQuadrature of the sum of the degrees of the source and
  /// the target space, CreateMonteCarloQuadrature, or the intersection of two
  /// meshes, which carries the points bound on the source's geometry); unset,
  /// the Gauss rule of degree 2
  std::optional<Quadrature> quadrature;
  /// project only the residual of the source from its interpolation onto the
  /// target space (a control variate), which is what the Monte Carlo
  /// right-hand side then estimates: far less variance, and exact for
  /// fields of the target space
  bool control_variate = false;
  /// Consistent: the mass matrix E^T W E; Lumped: its row sums, a diagonal
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
  LinearSolveOptions solve;
  /// map from the target's coordinate system to the source's, or null if they
  /// are the same; the source is evaluated at the mapped points
  std::shared_ptr<const CoordinateMap> to_source;
};

/// Galerkin L2 projection of any evaluable source onto target_space: solves
/// the mass system M x = b with b = E^T W f, where E evaluates the target
/// space at the points of a quadrature, W holds its weights and f is the
/// source at the points (through PointEvaluator::EvaluateTranspose). The
/// mass matrix M = E^T W E takes the target space's natural quadrature
/// (FunctionSpace::CreateQuadrature), and is exact for it. A field with several
/// components is projected jointly (the components do not couple). The
/// operator holds `source` and `target_space`.
/// @throws pcms_error if `source` or `target_space` is null, the source has no
///         geometry, or the target space's evaluator does not support the
///         transpose or has no natural quadrature, or `options.solve` asks for
///         the Assembled solve in a build without PCMS_ENABLE_PETSC
std::unique_ptr<TransferOperator<Real>> CreateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space,
  ProjectionOptions options = {});

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
