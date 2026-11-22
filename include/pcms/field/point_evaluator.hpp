#ifndef PCMS_FIELD_POINT_EVALUATOR_HPP
#define PCMS_FIELD_POINT_EVALUATOR_HPP

#include "pcms/field/differential_operator.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <string>
#include "pcms/field/detail/point_evaluator.hpp"
#include <memory>
#include <variant>

namespace pcms
{

template <typename T>
class Field;

/// @brief The evaluation at the points as a CSR matrix whose columns are the
///        global ids of the DOF holders.
template <typename T>
struct PointMatrix
{
  Kokkos::View<const LO*, DeviceMemorySpace> row_offsets; // [points + 1]
  Kokkos::View<const GO*, DeviceMemorySpace> columns;     // [entries]
  Kokkos::View<const T*, DeviceMemorySpace> values;       // [entries]
  /// local DOF holders of the columns
  Kokkos::View<const LO*, DeviceMemorySpace> holders; // [entries]
};

/// @brief Evaluates fields at the points of a BoundPoints.
///
/// Created by PointEvaluatorFactory::CreatePointEvaluator. Evaluate writes
/// OutputBasis(stored basis of the field): the basis `OutputBasis` gives for
/// the evaluator's operator. The output has one row per query point and the
/// field's number of components times the coordinate dimension per added
/// index as columns. EvaluateImpl writes the rows of the
/// BoundPoints's FoundPoints(); in OutOfBoundsMode::FILL, Evaluate then sets
/// the rows of MissingPoints() to the fill value.
/// @tparam T field value type
/// @tparam LayoutPolicy layout of the output view
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class PointEvaluator
{
public:
  /// Writes the field's values at the query points into `values`.
  /// @param field a field of a compatible space; today that is the space
  ///        this evaluator was created from (a field of its layout), or any
  ///        field for an evaluator that accepts every layout
  /// @param values output, tagged with OutputBasis of the field's basis
  /// @throws std::invalid_argument if `values` is not shaped
  ///         [NumPoints()][number of output components]
  /// @throws pcms_error if `values` does not claim the basis this evaluator
  ///         writes for the field, or if the evaluator accepts the fields of
  ///         one layout and `field` is not one of them
  void Evaluate(const Field<T>& field,
                ValueView<T, DeviceMemorySpace, LayoutPolicy> values) const
  {
    CheckArguments("PointEvaluator::Evaluate", field, values);
    EvaluateImpl(field, values);
    const auto filled = FilledPoints();
    if (filled.extent(0) > 0) {
      detail::FillRows(values.GetValues(), filled,
                       static_cast<T>(bound_->Policy().fill_value));
    }
  }

  /// @return whether EvaluateTranspose and the Value operator are supported
  [[nodiscard]] virtual bool SupportsTranspose() const noexcept
  {
    return false;
  }

  /// @brief Adds the transpose of the evaluation to a functional:
  ///        functional += E^T values, with E the map from a field's DOFs to
  ///        its values at the points, summed over the ranks (collective).
  /// @param values shaped and tagged as Evaluate's output for a field of the
  ///        functional's basis
  /// @param functional a field of the evaluator's layout, accumulated into
  /// @throws pcms_error if the transpose is not supported, `functional` is
  ///         not of the evaluator's layout or `values` has another basis
  /// @throws std::invalid_argument if `values` is not shaped like
  ///         Evaluate's output
  void EvaluateTranspose(
    ValueView<const T, DeviceMemorySpace, LayoutPolicy> values,
    Field<T>& functional) const
  {
    if (!SupportsTranspose() || op_ != DifferentialOperator::Value) {
      throw pcms_error("PointEvaluator::EvaluateTranspose: not supported");
    }
    CheckArguments("PointEvaluator::EvaluateTranspose", functional, values);
    EvaluateTransposeImpl(values.GetValues(), functional);
  }

  /// @return whether GetMatrix is supported
  [[nodiscard]] virtual bool SupportsMatrix() const noexcept { return false; }

  /// @brief The evaluation as a sparse matrix E: component c of the values is
  ///        E times component c of the DOFs. For the Gradient operator, row
  ///        p * d + j holds partial j of point p (d the output columns).
  /// @throws pcms_error if the matrix is not supported
  [[nodiscard]] PointMatrix<T> GetMatrix() const
  {
    if (!SupportsMatrix()) {
      throw pcms_error("PointEvaluator::GetMatrix: not supported");
    }
    return GetMatrixImpl();
  }

  /// @return number of query points
  [[nodiscard]] LO NumPoints() const noexcept { return bound_->NumPoints(); }

  /// @param stored value basis of the evaluated field
  /// @return value basis this evaluator writes for `stored`
  /// @throws pcms_error if the operator is not supported for `stored`
  [[nodiscard]] virtual ValueBasis OutputBasis(const ValueBasis& stored) const
  {
    switch (op_) {
      case DifferentialOperator::Value: return stored;
      case DifferentialOperator::Gradient:
        // the partials are taken where the points were located, in the
        // geometry's system; a bound map re-expresses them afterwards
        return detail::GradientBasis(
          stored, bound_->GetGeometry()->GetCoordinateSystem());
      default:
        throw pcms_error(std::string("PointEvaluator::OutputBasis: ") +
                         std::string(ToString(op_)) + " is not implemented");
    }
  }

  /// @return the operator this evaluator applies
  [[nodiscard]] DifferentialOperator Operator() const noexcept { return op_; }

  /// @return the layout whose fields Evaluate accepts, or null for any
  [[nodiscard]] const std::shared_ptr<const FieldLayout>& AcceptedLayout()
    const noexcept
  {
    return accepted_layout_;
  }

  /// @return indices of the query points whose rows Evaluate sets to the fill
  ///         value; empty unless the policy is OutOfBoundsMode::FILL
  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> FilledPoints()
    const noexcept
  {
    if (bound_->Policy().mode != OutOfBoundsMode::FILL) {
      return {};
    }
    return bound_->MissingPoints();
  }

  virtual ~PointEvaluator() noexcept = default;

protected:
  /// @param bound BoundPoints of the query points; not null
  /// @param op the operator this evaluator applies
  /// @param accepted_layout the layout whose fields EvaluateImpl can read, or
  ///        null to accept any field (layout-agnostic evaluators)
  explicit PointEvaluator(
    std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op = DifferentialOperator::Value,
    std::shared_ptr<const FieldLayout> accepted_layout = nullptr)
    : bound_(std::move(bound)),
      op_(op),
      accepted_layout_(std::move(accepted_layout))
  {
  }

  /// @return the BoundPoints of the query points
  [[nodiscard]] const BoundPoints& Bound() const noexcept { return *bound_; }

  /// @brief Writes the field's values at the BoundPoints's FoundPoints().
  /// @param field field to evaluate; its output shape and basis are checked
  /// @param values output view; rows of other points may be left unset
  virtual void EvaluateImpl(
    const Field<T>& field,
    ValueView<T, DeviceMemorySpace, LayoutPolicy> values) const = 0;

  /// @brief Adds E^T values to the functional; called only if
  ///        SupportsTranspose(), with checked arguments.
  virtual void EvaluateTransposeImpl(
    Rank2View<const T, DeviceMemorySpace, LayoutPolicy> /*values*/,
    Field<T>& /*functional*/) const
  {
  }

  /// @return the evaluation matrix; called only if SupportsMatrix()
  virtual PointMatrix<T> GetMatrixImpl() const { return {}; }

private:
  /// @brief Checks that `field` is of the accepted layout and that `values`
  ///        has the basis and shape of the output for it.
  /// @throws pcms_error if the layout or the basis differs
  /// @throws std::invalid_argument if the shape differs
  template <typename Values>
  void CheckArguments(const char* caller, const Field<T>& field,
                      const Values& values) const
  {
    if (accepted_layout_ != nullptr &&
        &field.GetLayout() != accepted_layout_.get()) {
      throw pcms_error(std::string(caller) +
                       ": the field is not a field of the function space "
                       "this evaluator was created from");
    }
    const auto num_points = static_cast<size_t>(NumPoints());
    // per block of the field's components (one per component of a scalar
    // field, one for a vector or tensor), the components of the basis
    // OutputBasis writes: the operator adds indices, a map can change the
    // dimension
    const ValueBasis stored = field.GetData().GetValueBasis();
    const ValueBasis written = OutputBasis(stored);
    auto num_components =
      static_cast<size_t>(field.GetLayout().GetNumComponents());
    for (int i = 0; i < stored.Rank(); ++i) {
      num_components /= static_cast<size_t>(stored.system->Dimension());
    }
    for (int i = 0; i < written.Rank(); ++i) {
      num_components *= static_cast<size_t>(written.system->Dimension());
    }
    if (values.extent(0) != num_points || values.extent(1) != num_components) {
      throw std::invalid_argument(
        std::string(caller) + ": the output view is shaped [" +
        std::to_string(values.extent(0)) + "][" +
        std::to_string(values.extent(1)) + "] but must be [" +
        std::to_string(num_points) + "][" + std::to_string(num_components) +
        "] (query points x output components)");
    }
    detail::CheckWrittenValueBasis(caller, values.GetBasis(), written);
  }

  std::shared_ptr<const BoundPoints> bound_;
  DifferentialOperator op_;
  // layout of the fields Evaluate accepts; null accepts any field
  std::shared_ptr<const FieldLayout> accepted_layout_;
};

// Variant types using default layout for device memory space
using PointEvaluatorVariant =
  std::variant<std::unique_ptr<PointEvaluator<int8_t>>,
               std::unique_ptr<PointEvaluator<int32_t>>,
               std::unique_ptr<PointEvaluator<int64_t>>,
               std::unique_ptr<PointEvaluator<float>>,
               std::unique_ptr<PointEvaluator<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_POINT_EVALUATOR_HPP
