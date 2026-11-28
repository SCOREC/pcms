#ifndef PCMS_TRANSFER_PETSC_SOLVER_HPP
#define PCMS_TRANSFER_PETSC_SOLVER_HPP

#include "pcms/transfer/linear_solver.hpp"
#include <memory>

namespace pcms::transfer::detail
{

/// @return a PETSc KSP on the assembled matrix of a square `form` (serial):
///         conjugate gradients preconditioned by its diagonal, to the
///         tolerance and iteration limit of `options`, unless the options of
///         its PETSc prefix say otherwise. Holds a reference to `form`.
/// @throws pcms_error if the layout is distributed or the test evaluator has
///         no evaluation rows
std::unique_ptr<LinearSolver> CreatePetscSolver(
  const BilinearForm& form, const pcms::detail::LinearSolveOptions& options);

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_PETSC_SOLVER_HPP
