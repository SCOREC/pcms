#ifndef PCMS_TRANSFER_FORMS_HPP
#define PCMS_TRANSFER_FORMS_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms::transfer::detail
{

// The internal forms of the Galerkin projections, written on the evaluation of
// a space at the points of a quadrature (E) and the quadrature weights (W);
// they know nothing of the kind of the space. Their interfaces are meant to
// be published later.

/// @brief The quadrature of a BilinearForm and the map to its trial.
struct FormOptions
{
  /// the quadrature, on the test space's geometry; unset: the test space's
  /// GetQuadrature()
  std::optional<Quadrature> quadrature;
  /// map from the test space's geometry to the trial's, or null
  std::shared_ptr<const CoordinateMap> to_trial;
};

/// @brief The linear form E^T W of a space: samples at the points of a
///        quadrature to a functional.
class LinearForm
{
public:
  /// @throws pcms_error if the quadrature is not on the space's geometry
  LinearForm(std::shared_ptr<const FunctionSpace> space, Quadrature quadrature);

  /// @brief functional = E^T W samples
  /// @throws pcms_error if the evaluator does not support the transpose
  void Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                Field<Real>& functional) const;

  /// @return W samples, in a scratch the form owns (valid until the next
  ///         call)
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> Weight(
    ValueView<const Real, DeviceMemorySpace> samples) const;

  [[nodiscard]] const std::shared_ptr<const FunctionSpace>& Space()
    const noexcept
  {
    return space_;
  }
  [[nodiscard]] const Quadrature& GetQuadrature() const noexcept
  {
    return quadrature_;
  }
  /// @return E: the space at the quadrature's points
  [[nodiscard]] const PointEvaluator<Real>& Evaluator() const noexcept
  {
    return *evaluator_;
  }

private:
  std::shared_ptr<const FunctionSpace> space_;
  Quadrature quadrature_;
  std::unique_ptr<PointEvaluator<Real>> evaluator_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> weighted_;
};

/// @brief B = E_test^T W E_trial: the trial at the points of a quadrature,
///        the weights, then the transpose of the test evaluation.
class BilinearForm
{
public:
  /// @param trial a space, or any evaluator source with a geometry
  /// @throws pcms_error as LinearForm, if the trial has no geometry, or if an
  ///         evaluator cannot be created
  BilinearForm(std::shared_ptr<const FunctionSpace> test,
               std::shared_ptr<const PointEvaluatorFactory> trial,
               FormOptions options = {});

  /// @brief The square form of `space`.
  explicit BilinearForm(std::shared_ptr<const FunctionSpace> space,
                        FormOptions options = {});

  /// @return the form of `trial` on the same test points; the map to the
  ///         trial is not kept
  [[nodiscard]] BilinearForm WithTrial(
    std::shared_ptr<const PointEvaluatorFactory> trial) const;

  /// @brief functional = B u, a field of the test space's layout
  /// @throws pcms_error if the test evaluator does not support the transpose
  void Apply(const Field<Real>& u, Field<Real>& functional) const;

  /// @brief functional += B^T v, a functional of the trial's layout
  /// @throws pcms_error unless SupportsTranspose()
  void ApplyTranspose(const Field<Real>& v, Field<Real>& functional) const;

  /// @return whether the test evaluator supports the transpose
  [[nodiscard]] bool SupportsApply() const noexcept
  {
    return test_->Evaluator().SupportsTranspose();
  }

  /// @return whether the trial evaluator supports the transpose
  [[nodiscard]] bool SupportsTranspose() const noexcept
  {
    return TrialEvaluator().SupportsTranspose();
  }

  /// @return v^T B u
  [[nodiscard]] Real Pair(const Field<Real>& v, const Field<Real>& u) const;

  /// @return whether the form is square and its evaluator gives the matrix E
  ///         (the diagonal and the COO assembly need it)
  [[nodiscard]] bool SupportsMatrix() const noexcept;

  /// @brief diagonal = the diagonal of a square form at the owned DOF
  ///        holders
  /// @throws pcms_error unless SupportsMatrix(), or if the layout is
  ///         distributed
  void AssembleDiagonal(Field<Real>& diagonal) const;

  /// @brief lumped = B 1, the row sums of a square form
  /// @throws pcms_error if the form is not square or a basis function's
  ///         integral is not positive
  void AssembleLumped(Field<Real>& lumped) const;

  /// @brief A square form as coordinate entries with global DOF-holder ids;
  ///        entries with equal indices are to be summed.
  struct Coo
  {
    Kokkos::View<GO*, DeviceMemorySpace> rows;
    Kokkos::View<GO*, DeviceMemorySpace> columns;
    Kokkos::View<Real*, DeviceMemorySpace> values;
  };
  /// @throws pcms_error unless SupportsMatrix()
  [[nodiscard]] Coo AssembleCoo() const;

  [[nodiscard]] const std::shared_ptr<const FunctionSpace>& TestSpace()
    const noexcept
  {
    return test_->Space();
  }
  [[nodiscard]] const Quadrature& GetQuadrature() const noexcept
  {
    return test_->GetQuadrature();
  }

private:
  /// @return the basis of the trial's values for a field stored in `stored`
  [[nodiscard]] ValueBasis TrialBasis(const ValueBasis& stored) const
  {
    return TrialEvaluator().OutputBasis(stored);
  }

  /// @return E_test v, in a scratch the form owns
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> EvaluateTest(
    const Field<Real>& v) const;

  /// @return E_trial u, in a scratch the form owns
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> EvaluateTrial(
    const Field<Real>& u) const;

  /// @return the sum over the points of w a . b
  [[nodiscard]] Real Contract(ValueView<const Real, DeviceMemorySpace> a,
                              ValueView<const Real, DeviceMemorySpace> b) const;

  // the evaluations and the contraction, for the projections and integrals
  friend ValueBasis TrialBasis(const BilinearForm& form,
                               const ValueBasis& stored)
  {
    return form.TrialBasis(stored);
  }
  friend ValueView<const Real, DeviceMemorySpace> EvaluateTest(
    const BilinearForm& form, const Field<Real>& v)
  {
    return form.EvaluateTest(v);
  }
  friend ValueView<const Real, DeviceMemorySpace> EvaluateTrial(
    const BilinearForm& form, const Field<Real>& u)
  {
    return form.EvaluateTrial(u);
  }
  friend Real Contract(const BilinearForm& form,
                       ValueView<const Real, DeviceMemorySpace> a,
                       ValueView<const Real, DeviceMemorySpace> b)
  {
    return form.Contract(a, b);
  }

  BilinearForm(std::shared_ptr<const LinearForm> test,
               std::shared_ptr<const PointEvaluatorFactory> trial,
               FormOptions options);

  [[nodiscard]] const PointEvaluator<Real>& TrialEvaluator() const noexcept
  {
    return trial_evaluator_ != nullptr ? *trial_evaluator_ : test_->Evaluator();
  }

  std::shared_ptr<const LinearForm> test_;
  std::shared_ptr<const PointEvaluatorFactory> trial_;
  FormOptions options_;
  // null when the trial is the test space without a map
  std::unique_ptr<PointEvaluator<Real>> trial_evaluator_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> test_values_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> trial_values_;
};

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_FORMS_HPP
