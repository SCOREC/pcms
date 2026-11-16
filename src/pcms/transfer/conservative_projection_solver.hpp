#ifndef PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP
#define PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <petscksp.h>

#include <pcms/field/field.h>
#include <pcms/field/point_evaluator.h>
#include <pcms/transfer/bilinear_form_integrator.hpp>
#include <pcms/transfer/linear_form_integrator.hpp>
#include <pcms/utility/arrays.h>

namespace pcms
{

// Cached Galerkin projection solver for M*x = f.
//
// The mass matrix is obtained from mass_integrator (which owns and assembled
// it) and the KSP is set up (preconditioned/factored) at construction. The RHS
// integrator is also bound at construction: its integration points fix the
// size of the cached sample buffer, which is allocated once and reused. Each
// Solve() call only fills that buffer, assembles a fresh RHS vector, and calls
// KSPSolve, reusing the cached factorization.
//
// Binding one RHS integrator per solver matches every current consumer (each
// builds M and the integrator together and reuses them in lockstep) and lets
// the sample buffer be sized exactly once. Reusing a factored M with a
// different integrator is intentionally not supported; if that need arises,
// add a setter that rebinds the integrator and resizes the buffer.
//
// Non-copyable and non-movable because it owns PETSc handles.
class GalerkinProjectionSolver
{
public:
  // Retrieves the assembled matrix from mass_integrator, sets up the KSP, and
  // sizes the sample buffer from rhs_integrator's integration points. PETSc
  // reference-counts the matrix, so mass_integrator need not outlive this
  // solver, but rhs_integrator must (it is assembled into on every Solve).
  GalerkinProjectionSolver(const BilinearFormIntegrator& mass_integrator,
                           LinearFormIntegrator& rhs_integrator);
  ~GalerkinProjectionSolver();

  GalerkinProjectionSolver(const GalerkinProjectionSolver&) = delete;
  GalerkinProjectionSolver& operator=(const GalerkinProjectionSolver&) = delete;

  // Evaluates source_field at the bound integrator's integration points via
  // evaluator into the cached buffer, assembles the RHS, solves M*x = f using
  // the cached KSP, and returns nodal values on the target mesh.
  //
  // May be called repeatedly with different source fields, provided the target
  // mesh (and therefore M) is unchanged.
  Omega_h::Reals Solve(const PointEvaluator<Real>& evaluator,
                       const Field<Real>& source_field) const;

  // Assembles the RHS from values already sampled at the bound integrator's
  // integration points (shape [num_points][num_components]) and solves
  // M*x = f using the cached KSP. Used by workflows that pre-process the
  // sampled values (e.g. control-variate residuals) before assembly.
  Omega_h::Reals Solve(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) const;

private:
  PetscInt nverts_;
  KSP ksp_ = nullptr;
  LinearFormIntegrator* rhs_integrator_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> sampled_values_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP
