#ifndef PCMS_TRANSFER_DETAIL_PROJECTION_HPP
#define PCMS_TRANSFER_DETAIL_PROJECTION_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <optional>
#include <string>

namespace pcms::detail
{

/// @brief Options of the mass-system solve of a projection
struct LinearSolveOptions
{
  /// @brief Solution method
  enum class Method
  {
    /// Assembled when PCMS has PETSc, the layout is on one rank and the
    /// evaluator has evaluation rows; MatrixFree otherwise
    Automatic,
    /// Conjugate gradients preconditioned by the lumped mass, matrix free
    MatrixFree,
    /// PETSc KSP on the assembled matrix (conjugate gradients preconditioned
    /// by its diagonal, unless the PETSc options say otherwise); needs
    /// PCMS_ENABLE_PETSC
    Assembled
  };
  /// Solution method
  Method method = Method::Automatic;
  /// Relative residual tolerance
  Real relative_tolerance = 1e-12;
  /// Iteration limit; 0 means the number of unknowns
  int max_iterations = 0;
  /// Assembled: PETSc options prefix of the KSP
  std::string petsc_options_prefix;
};

struct ProjectionOptions
{
  /// the quadrature of the right-hand side, on the target space's geometry
  /// (e.g. CreateGaussQuadrature of the sum of the degrees of the source and
  /// the target space, CreateMonteCarloQuadrature, or the intersection of two
  /// meshes, whose `other` names the cells of the source's geometry the points
  /// lie in, to evaluate the source without a search); unset, as
  /// pcms::CreateProjection
  std::optional<Quadrature> quadrature;
  /// project only the residual of the source from its interpolation onto the
  /// target space (a control variate), which is what the Monte Carlo
  /// right-hand side then estimates: far less variance, and exact for
  /// fields of the target space; for a target whose DOFs are not values at
  /// points, the control variate is the projection with the target's natural
  /// quadrature. Unset: on when the quadrature is not polynomial-exact
  /// (Quadrature::degree < 0, e.g. Monte Carlo), off otherwise
  std::optional<bool> control_variate;
  /// Consistent: the mass matrix E^T W E; Lumped: its row sums, a diagonal
  MassMatrixType mass_matrix_type = MassMatrixType::Consistent;
  LinearSolveOptions solve;
  /// map from the target's coordinate system to the source's, or null if they
  /// are the same; the source is evaluated at the mapped points
  std::shared_ptr<const CoordinateMap> to_source;
};

/// @brief pcms::CreateProjection with the choices of `options`.
/// @throws pcms_error as pcms::CreateProjection, or if `options.solve` asks
///         for the Assembled solve in a build without PCMS_ENABLE_PETSC
std::unique_ptr<TransferOperator<Real>> CreateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space, ProjectionOptions options);

} // namespace pcms::detail

#endif // PCMS_TRANSFER_DETAIL_PROJECTION_HPP
