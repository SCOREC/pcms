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

/// @brief The quadrature of a BilinearForm and the map to its trial.
struct FormOptions
{
  /// the quadrature, on the test space's geometry; unset: the test space's
  /// GetQuadrature()
  std::optional<Quadrature> quadrature;
  /// map from the test space's geometry to the trial's, or null
  std::shared_ptr<const CoordinateMap> to_trial;
};

/// @brief B = E_test^T W E_trial: the trial evaluated at the points of a
///        quadrature on the test space, weighted, then the transpose of the
///        test evaluation. The evaluations return views of scratch the form
///        owns, valid until the next evaluation of the same side.
class BilinearForm
{
public:
  /// @param trial a space, or any evaluator source with a geometry
  /// @throws pcms_error if the quadrature is not on the test space's
  ///         geometry, the trial has no geometry or an evaluator cannot be
  ///         created
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
  void Apply(const Field<Real>& u, Field<Real>& functional) const;

  /// @brief functional = B^T v + beta functional, a functional of the
  ///        trial's layout; beta == 0 does not read the functional
  void ApplyTranspose(const Field<Real>& v, Field<Real>& functional,
                      Real beta = 0) const;

  /// @brief functional = E_test^T W samples
  void Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                Field<Real>& functional) const;

  /// @return v^T B u, over the ranks of the quadrature's geometry
  ///         (collective)
  [[nodiscard]] Real Pair(const Field<Real>& v, const Field<Real>& u) const;

  /// @return E_test v
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> EvaluateTest(
    const Field<Real>& v) const;

  /// @return E_trial u
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> EvaluateTrial(
    const Field<Real>& u) const;

  /// @return the sum over the points of w a . b, over the ranks of the
  ///         quadrature's geometry (collective)
  [[nodiscard]] Real Contract(ValueView<const Real, DeviceMemorySpace> a,
                              ValueView<const Real, DeviceMemorySpace> b) const;

  /// @return the basis of the trial's values for a field stored in `stored`
  [[nodiscard]] ValueBasis TrialBasis(const ValueBasis& stored) const
  {
    return TrialEvaluator().OutputBasis(stored);
  }

  /// @return whether the trial is the test space, without a map
  [[nodiscard]] bool IsSquare() const noexcept
  {
    return trial_evaluator_ == nullptr;
  }

  [[nodiscard]] const std::shared_ptr<const FunctionSpace>& TestSpace()
    const noexcept
  {
    return test_->space;
  }
  [[nodiscard]] const Quadrature& GetQuadrature() const noexcept
  {
    return test_->quadrature;
  }
  /// @return E_test: the test space at the quadrature's points
  [[nodiscard]] const PointEvaluator<Real>& TestEvaluator() const noexcept
  {
    return *test_->evaluator;
  }
  /// @return E_trial: the trial at the quadrature's points
  [[nodiscard]] const PointEvaluator<Real>& TrialEvaluator() const noexcept
  {
    return trial_evaluator_ != nullptr ? *trial_evaluator_ : *test_->evaluator;
  }

private:
  /// The test side, shared by the forms WithTrial makes.
  struct Test
  {
    std::shared_ptr<const FunctionSpace> space;
    Quadrature quadrature;
    std::unique_ptr<PointEvaluator<Real>> evaluator;
  };

  BilinearForm(std::shared_ptr<const Test> test,
               std::shared_ptr<const PointEvaluatorFactory> trial,
               FormOptions options);

  /// @return W samples, in scratch the form owns
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> Weight(
    ValueView<const Real, DeviceMemorySpace> samples) const;

  std::shared_ptr<const Test> test_;
  std::shared_ptr<const PointEvaluatorFactory> trial_;
  FormOptions options_;
  // null when the trial is the test space without a map
  std::unique_ptr<PointEvaluator<Real>> trial_evaluator_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> test_values_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> trial_values_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> weighted_;
};

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_FORMS_HPP
