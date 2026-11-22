#ifndef PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP

#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/field/point_rows.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>

#include <memory>
#include <set>

namespace pcms::uniform_grid
{

namespace detail
{

/// @brief The rows of the order 0 or 1 Lagrange fields of a uniform grid at
///        points.
/// @param wrapped bit d set if the last node of axis d is its first
/// @param cells the cell of each point
/// @param t the position of each point in its cell
/// @param missing points whose rows are zero
/// @return the rows
template <unsigned Dim>
pcms::detail::PointRows BuildGridRows(
  UniformGrid<Dim> grid, unsigned wrapped, int order, bool gradient,
  Kokkos::View<const LO*, DeviceMemorySpace> cells,
  Kokkos::View<const Real**, DeviceMemorySpace> t,
  Kokkos::View<const LO*, DeviceMemorySpace> missing)
{
  const int stride = order == 0 ? 1 : 1 << Dim;
  const int num_out = gradient ? static_cast<int>(Dim) : 1;
  const auto num_points = static_cast<LO>(cells.extent(0));
  const auto size = static_cast<std::size_t>(num_points) * num_out * stride;
  Kokkos::View<LO*, DeviceMemorySpace> holders(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniform_grid_holders"),
    size);
  Kokkos::View<Real*, DeviceMemorySpace> weights(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniform_grid_weights"),
    size);
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
      LO nodes[Dim];
      LO vertex_stride[Dim];
      Real h[Dim];
      LO next = 1;
      for (unsigned d = 0; d < Dim; ++d) {
        nodes[d] = grid.divisions[d] + (((wrapped >> d) & 1U) ? 0 : 1);
        vertex_stride[d] = next;
        next *= nodes[d];
        h[d] = grid.edge_length[d] / grid.divisions[d];
      }
      for (int corner = 0; corner < stride; ++corner) {
        LO vertex = 0;
        for (unsigned d = 0; d < Dim; ++d) {
          const LO node = reversed[Dim - 1 - d] + ((corner >> d) & 1);
          vertex += (node % nodes[d]) * vertex_stride[d];
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
  return {pcms::detail::StrideOffsets(num_points * num_out, stride), holders,
          weights, num_out};
}

} // namespace detail

/// @brief Point evaluation of order 0 (cell value) or order 1 (multilinear)
///        fields on a uniform grid and of their gradients.
template <unsigned Dim>
std::shared_ptr<const PointEvaluatorFactory> CreateLagrangeEvaluatorFactory(
  std::shared_ptr<const UniformGridGeometry<Dim>> geometry,
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout)
{
  using Bound = pcms::detail::UniformGridBinding<Dim>;
  const unsigned wrapped = geometry->WrappedAxes();
  return std::make_shared<pcms::detail::RowsEvaluatorFactory<Bound>>(
    std::move(geometry), layout, std::set{DifferentialOperator::Gradient},
    [grid = layout->GetGrid(), wrapped,
     order = layout->GetOrder()](const Bound& bound, DifferentialOperator op) {
      return detail::BuildGridRows<Dim>(
        grid, wrapped, order, op == DifferentialOperator::Gradient,
        bound.CellIndices(), bound.Parametric(),
        bound.Policy().mode == OutOfBoundsMode::FILL
          ? bound.MissingPoints()
          : Kokkos::View<const LO*, DeviceMemorySpace>{});
    });
}

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_EVALUATOR_HPP
