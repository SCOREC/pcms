#ifndef PCMS_FIELD_POINT_ROWS_HPP
#define PCMS_FIELD_POINT_ROWS_HPP

// Private: the evaluation of a field at a set of bound points as sparse rows
// of DOF weights, and the evaluator and evaluator source built on them.

#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/detail/point_evaluator_factory.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <functional>
#include <memory>
#include <set>
#include <utility>

namespace pcms::detail
{

/// @brief The evaluation at a set of points as CSR rows: row p * num_out + r
///        holds the (local DOF holder, weight) entries of output r at point
///        p, the same for every component; num_out is 1 for Value and the
///        dimension for Gradient
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

/// @brief The rows as a PointMatrix.
/// @param gids the global id of each holder the rows use
template <typename Gids>
PointMatrix<Real> ToPointMatrix(const PointRows& rows, Gids gids)
{
  const auto holders = rows.holders;
  Kokkos::View<GO*, DeviceMemorySpace> columns(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "point_rows_columns"),
    holders.extent(0));
  ForEach(
    "pcms::detail::ToPointMatrix", static_cast<LO>(holders.extent(0)),
    KOKKOS_LAMBDA(const LO k) { columns(k) = gids[holders(k)]; });
  return {rows.offsets, columns, rows.values, holders};
}

/// @brief The global ids of a layout's DOF holders, on the device.
inline Kokkos::View<GO*, DeviceMemorySpace> DeviceGids(
  const FieldLayout& layout)
{
  const auto gids = layout.GetGidsHost();
  Kokkos::View<GO*, HostMemorySpace> host("point_rows_gids", gids.size());
  for (std::size_t i = 0; i < gids.size(); ++i) {
    host(i) = gids[i];
  }
  return Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, host);
}

/// @brief Evaluates fields through the PointRows of the bound points, with
///        the transpose (Value) and the matrix.
class RowsPointEvaluator : public PointEvaluator<Real>
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
  [[nodiscard]] bool SupportsMatrix() const noexcept override { return true; }

  // Public because CUDA forbids extended lambdas in non-public members.
  void EvaluateImpl(const Field<Real>& field,
                    ValueView<Real, DeviceMemorySpace> values) const override
  {
    Apply(rows_, field.GetDOFHolderData().GetValues(), values.GetValues());
  }

  void EvaluateTransposeImpl(Rank2View<const Real, DeviceMemorySpace> values,
                             Field<Real>& functional) const override
  {
    const auto& layout = functional.GetLayout();
    Kokkos::View<Real**, DeviceMemorySpace> total("point_rows_transpose",
                                                  layout.GetNumDofHolder(),
                                                  layout.GetNumComponents());
    ApplyTranspose(rows_, values, total);
    AddToField(functional, total);
  }

  [[nodiscard]] PointMatrix<Real> GetMatrixImpl() const override
  {
    return ToPointMatrix(rows_, DeviceGids(*AcceptedLayout()));
  }

protected:
  [[nodiscard]] const PointRows& Rows() const noexcept { return rows_; }

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
  RowsEvaluatorFactory(std::shared_ptr<const Geometry> geometry,
                       std::shared_ptr<const FieldLayout> layout,
                       std::set<DifferentialOperator> supported,
                       MakeRows make_rows)
    : TypedPointEvaluatorFactory<Real, Bound>(std::move(geometry),
                                              std::move(supported)),
      layout_(std::move(layout)),
      make_rows_(std::move(make_rows))
  {
  }

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
};

} // namespace pcms::detail

#endif // PCMS_FIELD_POINT_ROWS_HPP
