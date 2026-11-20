#ifndef PCMS_TRANSFER_LINEAR_SOLVER_HPP
#define PCMS_TRANSFER_LINEAR_SOLVER_HPP

#include "pcms/field/field.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/transfer/projection.hpp"
#include <memory>
#include <optional>

namespace pcms::transfer::detail
{

/// out = numerator / denominator, entry by entry over the owned values of
/// fields of one layout (a diagonal solve: the Jacobi preconditioner, the
/// lumped mass).
void Divide(const Field<Real>& numerator, const Field<Real>& denominator,
            Field<Real>& out);

/// Solves A x = b for the mass operator A of a BilinearForm: b a functional,
/// x coefficients, both fields of the form's layout.
///
/// MatrixFree: Jacobi-preconditioned conjugate gradients on all components
/// jointly, on the vector operations of the fields (field_ops.hpp) and the
/// form's Apply; the pairings are Apply(functional, function). Without the
/// diagonal (an evaluator with no matrix) the preconditioner is the identity.
/// Assembled: a PETSc matrix from the form's COO entries and a KSP, serial;
/// without PCMS_ENABLE_PETSC the constructor throws.
class LinearSolver
{
public:
  /// Holds a reference to `form`, which must outlive the solver.
  LinearSolver(const BilinearForm& form, LinearSolveOptions options);
  ~LinearSolver();
  LinearSolver(const LinearSolver&) = delete;
  LinearSolver& operator=(const LinearSolver&) = delete;

  /// @return the iterations taken (0 for Assembled)
  /// @throws pcms_error if the conjugate gradients do not converge
  int Solve(const Field<Real>& functional, Field<Real>& solution) const;

private:
  struct Petsc; // the PETSc objects of the Assembled method
  /// The vectors of the conjugate gradients, overwritten on every Solve.
  struct CgScratch
  {
    Function<Real> residual, direction, image, preconditioned;
  };
  const BilinearForm& form_;
  LinearSolveOptions options_;
  std::unique_ptr<Petsc> petsc_;
  // MatrixFree scratch
  mutable std::optional<Function<Real>> diagonal_;
  std::unique_ptr<CgScratch> cg_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> zeros_;
};

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_LINEAR_SOLVER_HPP
