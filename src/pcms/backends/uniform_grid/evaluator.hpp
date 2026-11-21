#ifndef PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP

#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"

#include <memory>
#include <utility>

namespace pcms::uniform_grid
{

namespace detail
{

/// @brief The evaluation of uniform-grid Lagrange fields at points as sparse
///        rows: row p * num_out + r holds `stride` (DOF holder, weight)
///        entries, the same for every component; num_out is 1 for Value and
///        the dimension for Gradient (row r: the partial along axis r).
struct UniformGridRows
{
  Kokkos::View<LO*, DeviceMemorySpace> holders;
  Kokkos::View<Real*, DeviceMemorySpace> weights;
  int stride = 1;
  int num_out = 1;
};

/// @brief Builds the rows of order `order` at the points.
/// @param cells the cell of each point
/// @param t the position of each point in its cell
/// @param missing points whose rows are zero
/// @return the rows
template <unsigned Dim>
UniformGridRows BuildGridRows(
  UniformGrid<Dim> grid, int order, bool gradient,
  Kokkos::View<const LO*, DeviceMemorySpace> cells,
  Kokkos::View<const Real**, DeviceMemorySpace> t,
  Kokkos::View<const LO*, DeviceMemorySpace> missing)
{
  UniformGridRows rows;
  rows.stride = order == 0 ? 1 : 1 << Dim;
  rows.num_out = gradient ? static_cast<int>(Dim) : 1;
  const int stride = rows.stride;
  const int num_out = rows.num_out;
  const auto num_points = static_cast<LO>(cells.extent(0));
  const auto size = static_cast<std::size_t>(num_points) * num_out * stride;
  rows.holders = Kokkos::View<LO*, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniform_grid_holders"),
    size);
  rows.weights = Kokkos::View<Real*, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniform_grid_weights"),
    size);
  const auto holders = rows.holders;
  const auto weights = rows.weights;
  pcms::detail::ForEach(
    "uniform_grid_rows", num_points, KOKKOS_LAMBDA(const LO i) {
      if (order == 0) {
        for (int r = 0; r < num_out; ++r) {
          holders(i * num_out + r) = cells(i);
          weights(i * num_out + r) = gradient ? 0 : 1;
        }
        return;
      }
      // cell indices in axis order (GetDimensionedIndex gives them reversed),
      // and the vertex numbering with the first axis fastest
      const auto reversed = grid.GetDimensionedIndex(cells(i));
      LO first_vertex = 0;
      LO vertex_stride[Dim];
      Real h[Dim];
      LO next = 1;
      for (unsigned d = 0; d < Dim; ++d) {
        vertex_stride[d] = next;
        next *= grid.divisions[d] + 1;
        h[d] = grid.edge_length[d] / grid.divisions[d];
        first_vertex += reversed[Dim - 1 - d] * vertex_stride[d];
      }
      for (int corner = 0; corner < stride; ++corner) {
        LO vertex = first_vertex;
        for (unsigned d = 0; d < Dim; ++d) {
          vertex += ((corner >> d) & 1) * vertex_stride[d];
        }
        for (int r = 0; r < num_out; ++r) {
          Real weight = 1;
          for (unsigned d = 0; d < Dim; ++d) {
            const bool upper = (corner >> d) & 1;
            weight *= gradient && static_cast<int>(d) == r
                        ? (upper ? 1 : -1) / h[d]
                        : (upper ? t(i, d) : 1 - t(i, d));
          }
          const auto k =
            (static_cast<std::size_t>(i) * num_out + r) * stride + corner;
          holders(k) = vertex;
          weights(k) = weight;
        }
      }
    });
  const auto block = static_cast<std::size_t>(num_out) * stride;
  pcms::detail::ForEach(
    "uniform_grid_rows_missing", static_cast<LO>(missing.extent(0)),
    KOKKOS_LAMBDA(const LO m) {
      const auto first = static_cast<std::size_t>(missing(m)) * block;
      for (std::size_t k = first; k < first + block; ++k) {
        weights(k) = 0;
      }
    });
  return rows;
}

/// @brief Applies the rows to the DOFs.
/// @param values output: values(p, c * num_out + r) = row (p, r) times
///        component c of `dofs`
template <typename Dofs, typename Values>
void ApplyGridRows(const UniformGridRows& rows, Dofs dofs, Values values)
{
  const auto holders = rows.holders;
  const auto weights = rows.weights;
  const int stride = rows.stride;
  const int num_out = rows.num_out;
  const int num_components = dofs.extent(1);
  pcms::detail::ForEach(
    "uniform_grid_apply_rows", static_cast<LO>(values.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      for (int r = 0; r < num_out; ++r) {
        const auto first = (static_cast<std::size_t>(i) * num_out + r) * stride;
        for (int c = 0; c < num_components; ++c) {
          Real sum = 0;
          for (int k = 0; k < stride; ++k) {
            sum += weights(first + k) * dofs(holders(first + k), c);
          }
          values(i, c * num_out + r) = sum;
        }
      }
    });
}

/// @brief Adds the transpose of the Value rows applied to `values` to
///        `total`.
template <typename Values>
void ApplyGridRowsTranspose(const UniformGridRows& rows, Values values,
                            Kokkos::View<Real**, DeviceMemorySpace> total)
{
  const auto holders = rows.holders;
  const auto weights = rows.weights;
  const int stride = rows.stride;
  const int num_components = total.extent(1);
  pcms::detail::ForEach(
    "uniform_grid_apply_rows_transpose", static_cast<LO>(values.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      const auto first = static_cast<std::size_t>(i) * stride;
      for (int k = 0; k < stride; ++k) {
        for (int c = 0; c < num_components; ++c) {
          Kokkos::atomic_add(&total(holders(first + k), c),
                             weights(first + k) * values(i, c));
        }
      }
    });
}

/// @brief The rows as a PointMatrix.
/// @param gids the global id of each DOF holder
inline PointMatrix<Real> GridRowsMatrix(
  const UniformGridRows& rows, Kokkos::View<const GO*, DeviceMemorySpace> gids)
{
  const auto num_rows = static_cast<LO>(rows.holders.extent(0) / rows.stride);
  const int stride = rows.stride;
  Kokkos::View<LO*, DeviceMemorySpace> offsets("uniform_grid_offsets",
                                               num_rows + 1);
  Kokkos::View<GO*, DeviceMemorySpace> columns("uniform_grid_columns",
                                               rows.holders.extent(0));
  const auto holders = rows.holders;
  pcms::detail::ForEach(
    "uniform_grid_offsets", num_rows + 1,
    KOKKOS_LAMBDA(const LO r) { offsets(r) = r * stride; });
  pcms::detail::ForEach(
    "uniform_grid_columns", static_cast<LO>(columns.extent(0)),
    KOKKOS_LAMBDA(const LO k) { columns(k) = gids(holders(k)); });
  return {offsets, columns, rows.weights};
}

} // namespace detail

/// @brief Evaluates order-0 (cell value) or order-1 (multilinear) fields on a
///        uniform grid and their gradients, with the transpose (Value) and the
///        matrix.
/// @tparam Dim grid dimension
/// @tparam LayoutPolicy layout of the output view
template <unsigned Dim = 2,
          typename LayoutPolicy =
            pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class UniformGridPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  /// @param layout uniform-grid layout of the evaluated fields
  /// @param binding cells of the query points
  UniformGridPointEvaluator(
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout,
    std::shared_ptr<const pcms::detail::UniformGridBinding<Dim>> binding,
    DifferentialOperator op = DifferentialOperator::Value)
    : PointEvaluator<Real, LayoutPolicy>(binding, op, layout),
      layout_(std::move(layout)),
      binding_(std::move(binding))
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  // TODO: add passkey or other mechanism to protect this
  void EvaluateImpl(
    const Field<Real>& field,
    ValueView<Real, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto dof_data = field.GetDOFHolderData().GetValues();
    PCMS_ALWAYS_ASSERT(dof_data.size() ==
                       static_cast<size_t>(layout_->GetNumDofHolder() *
                                           layout_->GetNumComponents()));
    // Every point has a (closest) cell, so all rows are computed; Evaluate
    // overwrites the rows of points outside the grid in FILL mode.
    detail::ApplyGridRows(Rows(), dof_data, out.GetValues());
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }

  // Public because CUDA forbids extended lambdas in non-public members.
  void EvaluateTransposeImpl(
    Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> values,
    Field<Real>& functional) const override
  {
    PCMS_FUNCTION_TIMER;
    Kokkos::View<Real**, DeviceMemorySpace> total("uniform_grid_transpose",
                                                  layout_->GetNumDofHolder(),
                                                  layout_->GetNumComponents());
    detail::ApplyGridRowsTranspose(Rows(), values, total);
    pcms::detail::AddToField(functional, total);
  }

  [[nodiscard]] bool SupportsMatrix() const noexcept override { return true; }

  PointMatrix<Real> GetMatrixImpl() const override
  {
    const auto gids = layout_->GetGidsHost();
    Kokkos::View<GO*, HostMemorySpace> host("uniform_grid_gids", gids.size());
    for (std::size_t i = 0; i < gids.size(); ++i) {
      host(i) = gids[i];
    }
    return detail::GridRowsMatrix(
      Rows(), Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, host));
  }

private:
  // built per call: a Rebind moves the points
  [[nodiscard]] detail::UniformGridRows Rows() const
  {
    return detail::BuildGridRows<Dim>(
      layout_->GetGrid(), layout_->GetOrder(),
      this->Operator() == DifferentialOperator::Gradient,
      binding_->CellIndices(), binding_->Parametric(), this->FilledPoints());
  }

  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
  std::shared_ptr<const pcms::detail::UniformGridBinding<Dim>> binding_;
};

/// @brief Point evaluation of fields on a structured uniform grid.
/// @tparam Dim grid dimension
template <unsigned Dim = 2>
class UniformGridEvaluatorFactory
  : public pcms::detail::TypedPointEvaluatorFactory<
      Real, pcms::detail::UniformGridBinding<Dim>>
{
  using Bound = pcms::detail::UniformGridBinding<Dim>;

public:
  /// @param geometry geometry of the layout's grid
  /// @param layout uniform-grid layout of the evaluated fields
  UniformGridEvaluatorFactory(
    std::shared_ptr<const UniformGridGeometry<Dim>> geometry,
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout)
    : pcms::detail::TypedPointEvaluatorFactory<Real, Bound>(
        std::move(geometry), {DifferentialOperator::Gradient}),
      layout_(std::move(layout))
  {
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound, DifferentialOperator op) const override
  {
    return std::make_unique<UniformGridPointEvaluator<Dim>>(
      layout_, std::move(bound), op);
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
};

using UniformGridEvaluatorFactory2D = UniformGridEvaluatorFactory<2>;

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP
