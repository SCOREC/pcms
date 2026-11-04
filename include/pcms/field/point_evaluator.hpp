#ifndef PCMS_FIELD_POINT_EVALUATOR_HPP
#define PCMS_FIELD_POINT_EVALUATOR_HPP

#include "pcms/field/differential_operator.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <string>

namespace pcms
{

template <typename T>
class Field;

namespace detail
{

/// @brief Sets every column of the listed rows of `values` to `value`.
/// @param values output view
/// @param rows indices of the rows to set
/// @param value value to write
template <typename T, typename LayoutPolicy>
void FillRows(Rank2View<T, DeviceMemorySpace, LayoutPolicy> values,
              Kokkos::View<const LO*, DeviceMemorySpace> rows, T value)
{
  const auto num_cols = static_cast<int>(values.extent(1));
  const auto fill = KOKKOS_LAMBDA(const LO k)
  {
    for (int c = 0; c < num_cols; ++c) {
      values(rows(k), c) = value;
    }
  };
  detail::ForEach("PointEvaluator::FillRows", static_cast<LO>(rows.extent(0)),
                  fill);
}

} // namespace detail

/// @brief Evaluates fields at the points of a bound source.
///
/// Created by PointEvaluatorFactory::CreatePointEvaluator. Evaluate writes
/// OutputBasis(stored basis of the field): the basis `OutputBasis` gives for
/// the evaluator's operator. The output has one row per query point and the
/// field's number of components times the coordinate dimension per added
/// index as columns. EvaluateImpl writes the rows of the
/// bound source's FoundPoints(); in OutOfBoundsMode::FILL, Evaluate then sets
/// the rows of MissingPoints() to the fill value.
/// @tparam T field value type
/// @tparam LayoutPolicy layout of the output view
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
/// Implementer API: not covered by the public-API stability promise.
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
    if (accepted_layout_ != nullptr &&
        &field.GetLayout() != accepted_layout_.get()) {
      throw pcms_error("PointEvaluator::Evaluate: the field is not a field of "
                       "the function space this evaluator was created from");
    }
    CheckOutputShape(field, values);
    detail::CheckWrittenValueBasis(
      "PointEvaluator::Evaluate", values.GetBasis(),
      OutputBasis(field.GetData().GetValueBasis()));
    EvaluateImpl(field, values);
    const auto filled = FilledPoints();
    if (filled.extent(0) > 0) {
      detail::FillRows(values.GetValues(), filled,
                       static_cast<T>(bound_->Policy().fill_value));
    }
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
  /// @param bound bound source of the query points; not null
  /// @param op the operator this evaluator applies
  /// @param accepted_layout the layout whose fields EvaluateImpl can read, or
  ///        null to accept any field (layout-agnostic evaluators)
  explicit PointEvaluator(
    std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op = DifferentialOperator::Value,
    std::shared_ptr<const FieldLayout> accepted_layout = nullptr)
    : bound_(std::move(bound)),
      op_(op),
      accepted_layout_(std::move(accepted_layout))
  {
  }

  /// @return the bound source of the query points
  [[nodiscard]] const BoundSource& Bound() const noexcept { return *bound_; }

  /// @brief Writes the field's values at the bound source's FoundPoints().
  /// @param field field to evaluate; its output shape and basis are checked
  /// @param values output view; rows of other points may be left unset
  virtual void EvaluateImpl(
    const Field<T>& field,
    ValueView<T, DeviceMemorySpace, LayoutPolicy> values) const = 0;

private:
  /// @brief Checks that `values` has one row per query point and one column
  ///        per component of the basis OutputBasis gives for the field.
  /// @throws std::invalid_argument otherwise
  void CheckOutputShape(
    const Field<T>& field,
    const ValueView<T, DeviceMemorySpace, LayoutPolicy>& values) const
  {
    const auto num_points = static_cast<size_t>(NumPoints());
    // the field's components times the dimension per index OutputBasis adds
    // for the evaluator's operator
    const ValueBasis stored = field.GetData().GetValueBasis();
    const ValueBasis written = OutputBasis(stored);
    auto num_components =
      static_cast<size_t>(field.GetLayout().GetNumComponents());
    for (int i = stored.Rank(); i < written.Rank(); ++i) {
      num_components *= static_cast<size_t>(written.system->Dimension());
    }
    if (values.extent(0) != num_points || values.extent(1) != num_components) {
      throw std::invalid_argument(
        "PointEvaluator::Evaluate: the output view is shaped [" +
        std::to_string(values.extent(0)) + "][" +
        std::to_string(values.extent(1)) + "] but must be [" +
        std::to_string(num_points) + "][" + std::to_string(num_components) +
        "] (query points x output components)");
    }
  }

  std::shared_ptr<const BoundSource> bound_;
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
