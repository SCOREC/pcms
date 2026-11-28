#ifndef PCMS_FIELD_POINT_ROWS_HPP
#define PCMS_FIELD_POINT_ROWS_HPP

// Private: the evaluation of a field at a set of bound points as sparse rows
// of DOF weights, and the evaluator and evaluator source built on them.

#include "pcms/field/detail/field_data_ops.hpp"
#include "pcms/field/detail/point_evaluator_factory.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <functional>
#include <memory>
#include <set>
#include <utility>

namespace pcms::detail
{

/// @brief The evaluation at a set of points as CSR rows: row p * num_out + r
///        holds the (DOF holder, weight) entries of output r at point p, the
///        same for every component; num_out is 1 for Value and the dimension
///        for Gradient. A holder indexes the layout's DOF holders, ghosts
///        included, in the order of FieldLayout::GetOwnedHost.
struct PointRows
{
  Kokkos::View<const LO*, DeviceMemorySpace> offsets; // [rows + 1]
  Kokkos::View<const LO*, DeviceMemorySpace> holders;
  Kokkos::View<const Real*, DeviceMemorySpace> values;
  int num_out = 1;
};

/// @brief The offsets of `num_rows` rows of `stride` entries each.
inline Kokkos::View<LO*, DeviceMemorySpace> StrideOffsets(LO num_rows,
                                                          int stride)
{
  Kokkos::View<LO*, DeviceMemorySpace> offsets(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "point_rows_offsets"),
    num_rows + 1);
  ForEach(
    "pcms::detail::StrideOffsets", num_rows + 1,
    KOKKOS_LAMBDA(const LO r) { offsets(r) = r * stride; });
  return offsets;
}

/// @brief out(p, c * num_out + r) = row (p, r) times component c of the DOFs.
/// @param dofs dofs(holder, c)
template <typename Dofs, typename Out>
void Apply(const PointRows& rows, Dofs dofs, Out out)
{
  const auto offsets = rows.offsets;
  const auto holders = rows.holders;
  const auto values = rows.values;
  const int num_out = rows.num_out;
  const int num_components = static_cast<int>(out.extent(1)) / num_out;
  ForEach(
    "pcms::detail::Apply", static_cast<LO>(out.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int r = 0; r < num_out; ++r) {
        const LO row = p * num_out + r;
        for (int c = 0; c < num_components; ++c) {
          Real sum = 0;
          for (LO j = offsets(row); j < offsets(row + 1); ++j) {
            sum += values(j) * dofs(holders(j), c);
          }
          out(p, c * num_out + r) = sum;
        }
      }
    });
}

/// @brief total(holder, c) += the transpose of the Value rows applied to
///        component c of `in`.
template <typename In, typename Total>
void ApplyTranspose(const PointRows& rows, In in, Total total)
{
  const auto offsets = rows.offsets;
  const auto holders = rows.holders;
  const auto values = rows.values;
  const int num_components = static_cast<int>(total.extent(1));
  ForEach(
    "pcms::detail::ApplyTranspose", static_cast<LO>(in.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (LO j = offsets(p); j < offsets(p + 1); ++j) {
        for (int c = 0; c < num_components; ++c) {
          Kokkos::atomic_add(&total(holders(j), c), values(j) * in(p, c));
        }
      }
    });
}

/// @brief An evaluator whose evaluation is a PointRows.
class RowsSource
{
public:
  [[nodiscard]] virtual const PointRows& Rows() const = 0;

  /// @return the global id of each DOF holder the rows index; a ghost has
  ///         its owner's id
  [[nodiscard]] virtual GlobalIDView<DeviceMemorySpace> HolderGids() const = 0;

protected:
  ~RowsSource() = default;
};

/// @return `evaluator` as a RowsSource, or null if its evaluation has no rows
[[nodiscard]] inline const RowsSource* AsRowsSource(
  const PointEvaluator<Real>& evaluator)
{
  return dynamic_cast<const RowsSource*>(&evaluator);
}

/// @brief Evaluates fields through the PointRows of the bound points, with
///        the transpose (Value).
class RowsPointEvaluator
  : public PointEvaluator<Real>
  , public RowsSource
{
public:
  /// @param layout layout of the evaluated fields, without ghosts
  RowsPointEvaluator(std::shared_ptr<const BoundPoints> bound,
                     std::shared_ptr<const FieldLayout> layout,
                     DifferentialOperator op, PointRows rows)
    : PointEvaluator<Real>(std::move(bound), op, std::move(layout)),
      rows_(std::move(rows))
  {
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }

  // Public because CUDA forbids extended lambdas in non-public members.
  void EvaluateImpl(const Field<Real>& field,
                    ValueView<Real, DeviceMemorySpace> values) const override
  {
    Apply(rows_, field.GetDOFHolderData().GetValues(), values.GetValues());
  }

  void EvaluateTransposeImpl(Rank2View<const Real, DeviceMemorySpace> values,
                             Field<Real>& functional, Real beta) const override
  {
    const auto& layout = functional.GetLayout();
    Kokkos::View<Real**, DeviceMemorySpace> total("point_rows_transpose",
                                                  layout.GetNumDofHolder(),
                                                  layout.GetNumComponents());
    if (beta != 0) {
      CopyDeviceRank2ViewToDeviceView(
        total, functional.GetDOFHolderData().GetValues());
      ScaleValues(total, beta, total);
    }
    ApplyTranspose(rows_, values, total);
    functional.SetDOFHolderDataUnchecked(MakeConstRank2View(total));
    functional.Synchronize();
  }

  [[nodiscard]] const PointRows& Rows() const noexcept override
  {
    return rows_;
  }

  /// @return the layout's owned gids; the layout has no ghosts
  [[nodiscard]] GlobalIDView<DeviceMemorySpace> HolderGids() const override
  {
    return AcceptedLayout()->GetGids();
  }

private:
  PointRows rows_;
};

/// @brief The evaluator source of a space whose evaluators are
///        RowsPointEvaluators of the rows a function makes from the bound
///        points.
/// @tparam Bound BoundPoints kind the geometry's Bind creates
template <typename Bound>
class RowsEvaluatorFactory final
  : public TypedPointEvaluatorFactory<Real, Bound>
{
public:
  using MakeRows = std::function<PointRows(const Bound&, DifferentialOperator)>;

  /// @param supported operators besides Value
  /// @param nodal whether a DOF is the value of the function at its DOF holder
  RowsEvaluatorFactory(std::shared_ptr<const Geometry> geometry,
                       std::shared_ptr<const FieldLayout> layout,
                       std::set<DifferentialOperator> supported,
                       MakeRows make_rows, bool nodal)
    : TypedPointEvaluatorFactory<Real, Bound>(std::move(geometry),
                                              std::move(supported)),
      layout_(std::move(layout)),
      make_rows_(std::move(make_rows)),
      nodal_(nodal)
  {
  }

  [[nodiscard]] bool HasNodalDOFs() const noexcept override { return nodal_; }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound, DifferentialOperator op) const override
  {
    auto rows = make_rows_(*bound, op);
    return std::make_unique<RowsPointEvaluator>(std::move(bound), layout_, op,
                                                std::move(rows));
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  MakeRows make_rows_;
  bool nodal_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_POINT_ROWS_HPP
