#ifndef PCMS_TRANSFER_LINEAR_SOLVER_HPP
#define PCMS_TRANSFER_LINEAR_SOLVER_HPP

#include "pcms/field/field.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/detail/projection.hpp"
#include <memory>
#include <optional>
#include <string>

namespace pcms::transfer::detail
{

/// @throws pcms_error naming `caller` unless the field data of the functions
///         of `space` provides the FieldData coefficient operations
void RequireCoefficientOperations(
  const std::shared_ptr<const FunctionSpace>& space, const std::string& caller);

/// @brief Solves B x = b for a square BilinearForm B: b a functional, x
///        coefficients, both fields of the form's test layout.
class LinearSolver
{
public:
  virtual ~LinearSolver() noexcept = default;

  /// @brief solution = B^-1 functional
  virtual void Solve(const Field<Real>& functional,
                     Field<Real>& solution) const = 0;
};

/// @brief Conjugate gradients on all components jointly, on the form's Apply
///        and the coefficient operations; B must be symmetric positive
///        definite. The preconditioner is the lumped form B 1, or the
///        identity where that is not positive.
class CgSolver final : public LinearSolver
{
public:
  /// Holds a reference to `form`, which must outlive the solver.
  CgSolver(const BilinearForm& form, Real relative_tolerance,
           int max_iterations);

  /// @throws pcms_error if the iterations do not converge
  void Solve(const Field<Real>& functional,
             Field<Real>& solution) const override;

  /// @return the iterations of the last Solve
  [[nodiscard]] int Iterations() const noexcept { return iterations_; }

private:
  /// The vectors in one value basis and the lumped form.
  struct Scratch
  {
    Function<Real> residual, direction, image, preconditioned;
    std::optional<Function<Real>> lumped;
  };

  /// @return the scratch in `basis`, made anew when the basis changes
  Scratch& InBasis(const ValueBasis& basis) const;

  const BilinearForm& form_;
  Real relative_tolerance_;
  int max_iterations_;
  mutable std::unique_ptr<Scratch> scratch_;
  mutable int iterations_ = 0;
};

/// @return the solver of `form`: x = b / (B 1) for a Lumped `type`, else the
///         method of `options`. Holds a reference to `form`, which must
///         outlive the solver.
/// @throws pcms_error if the form is not square, a lumped form is not
///         positive, or the method is not available
std::unique_ptr<LinearSolver> CreateLinearSolver(
  const BilinearForm& form, MassMatrixType type,
  const pcms::detail::LinearSolveOptions& options);

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_LINEAR_SOLVER_HPP
