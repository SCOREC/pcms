#ifndef PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H
#define PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H

#include "pcms/field/evaluator/uniform_grid_geometry.hpp"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include "pcms/transfer/linear_interpolant.hpp"
#include "pcms/transfer/multidimarray.hpp"

#include <memory>
#include <utility>

namespace pcms
{

/// @brief Evaluates order-0 (cell value) or order-1 (multilinear) fields on a
///        uniform grid.
/// @tparam Dim grid dimension
/// @tparam LayoutPolicy layout of the output view
template <unsigned Dim = 2,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class UniformGridPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  /// @param layout uniform-grid layout of the evaluated fields
  /// @param binding cells of the query points
  UniformGridPointEvaluator(
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout,
    std::shared_ptr<const UniformGridBinding<Dim>> binding)
    : PointEvaluator<Real, LayoutPolicy>(binding, DifferentialOperator::Value,
                                         layout),
      layout_(std::move(layout)),
      grid_(layout_->GetGrid()),
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
    const auto values = out.GetValues();
    const auto dof_data = field.GetDOFHolderData().GetValues();
    LO num_points = this->NumPoints();
    int n_comp = layout_->GetNumComponents();

    PCMS_ALWAYS_ASSERT(dof_data.size() ==
                       static_cast<size_t>(layout_->GetNumDofHolder() *
                                           layout_->GetNumComponents()));

    // Every point has a (closest) cell, so all rows are computed; Evaluate
    // overwrites the rows of points outside the grid in FILL mode.
    auto cell_indices = binding_->CellIndices();

    if (layout_->GetOrder() == 0) {
      Kokkos::parallel_for(
        "evaluate_order0_device",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_points, n_comp}),
        KOKKOS_LAMBDA(const LO i, const int c) {
          values(i, c) = dof_data(cell_indices(i), c);
        });
      return;
    }

    // Order-1: multilinear interpolation
    // Create local copy to avoid capturing reference in device lambdas in some
    // compilers
    auto grid_copy = grid_;
    Kokkos::View<LO**, DeviceMemorySpace> cell_dim_indices("cell_dim_indices",
                                                           num_points, Dim);
    Kokkos::parallel_for(
      "compute_cell_dim_indices",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_CLASS_LAMBDA(const LO i) {
        auto dim_idx = grid_copy.GetDimensionedIndex(cell_indices(i));
        for (unsigned d = 0; d < Dim; ++d) {
          cell_dim_indices(i, d) = dim_idx[d];
        }
      });

    auto cell_divisions = grid_.divisions;
    IntVecView dimensions_view("dimensions", Dim);
    Kokkos::parallel_for(
      "set_dimensions_view",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, Dim),
      KOKKOS_LAMBDA(const unsigned d) {
        // cell_dim_indices uses GetDimensionedIndex(), which returns indices in
        // the opposite order of the coordinates, so we reverse the order here.
        dimensions_view(d) = cell_divisions[Dim - 1 - d] + 1;
      });

    // the binding's coordinates are in axis order; cell_dim_indices is filled
    // from GetDimensionedIndex(), which returns indices in the opposite order
    RealMatView parametric_coords("parametric_coords", num_points, Dim);
    const auto parametric = binding_->Parametric();
    Kokkos::parallel_for(
      "compute_parametric_coords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_LAMBDA(const LO i) {
        for (unsigned d = 0; d < Dim; ++d) {
          parametric_coords(i, Dim - 1 - d) = parametric(i, d);
        }
      });

    IntMatView cell_indices_interp("cell_indices_interp", num_points, Dim);
    Kokkos::parallel_for(
      "copy_cell_dim_indices_to_interp",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_LAMBDA(const LO i) {
        for (unsigned d = 0; d < Dim; ++d)
          cell_indices_interp(i, d) = cell_dim_indices(i, d);
      });

    RealMatView values_interp("values_interp",
                              static_cast<size_t>(layout_->GetNumDofHolder()),
                              static_cast<size_t>(n_comp));
    Kokkos::parallel_for(
      "copy_dof_data_to_values_interp",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
        {0, 0}, {static_cast<LO>(layout_->GetNumDofHolder()), n_comp}),
      KOKKOS_LAMBDA(const LO i, const int c) {
        values_interp(i, c) = dof_data(i, c);
      });

    auto interpolator = RegularGridInterpolator(
      parametric_coords, values_interp, cell_indices_interp, dimensions_view);
    auto interpolated = interpolator.linear_interpolation();

    Kokkos::parallel_for(
      "evaluate_order1_device",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_points, n_comp}),
      KOKKOS_LAMBDA(const LO i, const int c) {
        values(i, c) = interpolated(i, c);
      });
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
  UniformGrid<Dim> grid_;
  std::shared_ptr<const UniformGridBinding<Dim>> binding_;
};

/// @brief Point evaluation of fields on a structured uniform grid.
/// @tparam Dim grid dimension
template <unsigned Dim = 2>
class UniformGridEvaluatorFactory
  : public TypedPointEvaluatorFactory<Real, UniformGridBinding<Dim>>
{
  using Bound = UniformGridBinding<Dim>;

public:
  /// @param geometry geometry of the layout's grid
  /// @param layout uniform-grid layout of the evaluated fields
  UniformGridEvaluatorFactory(
    std::shared_ptr<const UniformGridGeometry<Dim>> geometry,
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout)
    : TypedPointEvaluatorFactory<Real, Bound>(std::move(geometry)),
      layout_(std::move(layout))
  {
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<UniformGridPointEvaluator<Dim>>(layout_,
                                                            std::move(bound));
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
};

using UniformGridEvaluatorFactory2D = UniformGridEvaluatorFactory<2>;

} // namespace pcms

#endif // PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H
