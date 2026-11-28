#include "pcms/transfer/linear_solver.hpp"
#include "pcms/configuration.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/point_rows.hpp"
#ifdef PCMS_ENABLE_PETSC
#include "pcms/transfer/petsc_solver.hpp"
#endif
#include "pcms/utility/assert.hpp"
#include <utility>

namespace pcms::transfer::detail
{

namespace
{

constexpr Real RowSumTolerance = 1e-12;

// whether every owned coefficient, on every rank, exceeds a round-off
// fraction of the largest magnitude (collective)
bool AllPositive(const Field<Real>& values)
{
  return pcms::detail::Min(values) >
         RowSumTolerance * pcms::detail::MaxNorm(values);
}

Function<Real> CreateInBasis(const FunctionSpace& space,
                             const ValueBasis& basis)
{
  return space.CreateFunction<Real>("", basis.variance, basis.system,
                                    basis.scaling);
}

// lumped = B 1 in `basis`; whether it is positive
bool Lump(const BilinearForm& form, const ValueBasis& basis,
          Function<Real>& lumped)
{
  auto ones = CreateInBasis(*form.TestSpace(), basis);
  pcms::detail::Const(Real(1), ones);
  form.Apply(ones, lumped);
  return AllPositive(lumped);
}

// x = b / (B 1)
class LumpedSolver final : public LinearSolver
{
public:
  explicit LumpedSolver(const BilinearForm& form) : form_(form)
  {
    InBasis(
      form_.TestSpace()->CreateFunction<Real>().GetData().GetValueBasis());
  }

  void Solve(const Field<Real>& functional,
             Field<Real>& solution) const override
  {
    pcms::detail::Div(functional, InBasis(functional.GetData().GetValueBasis()),
                      solution);
  }

private:
  const Function<Real>& InBasis(const ValueBasis& basis) const
  {
    if (!lumped_ || lumped_->GetData().GetValueBasis() != basis) {
      lumped_.emplace(CreateInBasis(*form_.TestSpace(), basis));
      if (!Lump(form_, basis, *lumped_)) {
        lumped_.reset();
        throw pcms_error("LinearSolver: a row sum of the lumped form is not "
                         "positive; lumping needs a non-negative basis");
      }
    }
    return *lumped_;
  }

  const BilinearForm& form_;
  mutable std::optional<Function<Real>> lumped_;
};

// whether the assembled solve applies to `form`
bool CanAssemble([[maybe_unused]] const BilinearForm& form)
{
#ifdef PCMS_ENABLE_PETSC
  return !form.TestSpace()->GetLayout()->IsDistributed() &&
         pcms::detail::AsRowsSource(form.TestEvaluator()) != nullptr;
#else
  return false;
#endif
}

} // namespace

void RequireCoefficientOperations(
  const std::shared_ptr<const FunctionSpace>& space, const std::string& caller)
{
  auto probe = space->CreateFunction<Real>();
  try {
    pcms::detail::Const(Real(0), probe);
  } catch (const pcms_error&) {
    throw pcms_error(caller + ": the field data of the target space provides "
                              "no coefficient operations");
  }
}

CgSolver::CgSolver(const BilinearForm& form, Real relative_tolerance,
                   int max_iterations)
  : form_(form),
    relative_tolerance_(relative_tolerance),
    max_iterations_(max_iterations)
{
}

CgSolver::Scratch& CgSolver::InBasis(const ValueBasis& basis) const
{
  if (scratch_ == nullptr ||
      scratch_->residual.GetData().GetValueBasis() != basis) {
    const auto& space = *form_.TestSpace();
    scratch_ = std::make_unique<Scratch>(Scratch{
      CreateInBasis(space, basis), CreateInBasis(space, basis),
      CreateInBasis(space, basis), CreateInBasis(space, basis), std::nullopt});
    auto& lumped = scratch_->lumped.emplace(CreateInBasis(space, basis));
    if (!Lump(form_, basis, lumped)) {
      scratch_->lumped.reset();
    }
  }
  return *scratch_;
}

void CgSolver::Solve(const Field<Real>& functional, Field<Real>& solution) const
{
  auto& scratch = InBasis(solution.GetData().GetValueBasis());
  auto& r = scratch.residual;
  auto& p = scratch.direction;
  auto& ap = scratch.image;
  auto& z = scratch.preconditioned;
  const auto precondition = [&]() {
    if (scratch.lumped) {
      pcms::detail::Div(r, *scratch.lumped, z);
    } else {
      pcms::detail::Scale(Real(1), r, z);
    }
  };
  iterations_ = 0;
  pcms::detail::Const(Real(0), solution);
  pcms::detail::Scale(Real(1), functional, r);
  const Real target = relative_tolerance_ * pcms::detail::Norm2(functional);
  if (pcms::detail::Norm2(r) <= target) {
    return;
  }
  precondition();
  pcms::detail::Scale(Real(1), z, p);
  Real rz = pcms::detail::Apply(r, z);
  const int max_iterations =
    max_iterations_ > 0
      ? max_iterations_
      : static_cast<int>(form_.TestSpace()->GetLayout()->GlobalSize());
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    form_.Apply(p, ap);
    const Real curvature = pcms::detail::Apply(ap, p);
    if (!(curvature > 0)) {
      break;
    }
    const Real alpha = rz / curvature;
    pcms::detail::LinearSum(alpha, p, Real(1), solution, solution);
    pcms::detail::LinearSum(-alpha, ap, Real(1), r, r);
    if (pcms::detail::Norm2(r) <= target) {
      iterations_ = iteration;
      return;
    }
    precondition();
    const Real rz_next = pcms::detail::Apply(r, z);
    pcms::detail::LinearSum(Real(1), z, rz_next / rz, p, p);
    rz = rz_next;
  }
  throw pcms_error("CgSolver: the conjugate gradients did not converge");
}

std::unique_ptr<LinearSolver> CreateLinearSolver(
  const BilinearForm& form, MassMatrixType type,
  const pcms::detail::LinearSolveOptions& options)
{
  if (!form.IsSquare()) {
    throw pcms_error("CreateLinearSolver: the form is not square");
  }
  RequireCoefficientOperations(form.TestSpace(), "CreateLinearSolver");
  if (type == MassMatrixType::Lumped) {
    return std::make_unique<LumpedSolver>(form);
  }
  if (options.method == pcms::detail::LinearSolveOptions::Method::Assembled ||
      (options.method == pcms::detail::LinearSolveOptions::Method::Automatic &&
       CanAssemble(form))) {
#ifdef PCMS_ENABLE_PETSC
    return CreatePetscSolver(form, options);
#else
    throw pcms_error("CreateLinearSolver: the assembled solve needs PETSc; "
                     "configure PCMS with -DPCMS_ENABLE_PETSC=ON or use the "
                     "MatrixFree method");
#endif
  }
  return std::make_unique<CgSolver>(form, options.relative_tolerance,
                                    options.max_iterations);
}

} // namespace pcms::transfer::detail
