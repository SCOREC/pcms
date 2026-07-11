#ifndef PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H
#define PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H

#include "pcms/field/layout/uniform_grid.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include "pcms/transfer/linear_interpolant.hpp"
#include "pcms/transfer/multidimarray.hpp"

#include <memory>

namespace pcms
{

// ---------------------------------------------------------------------------
// Localization hint — computed once per query point set, reused across Evaluate
// ---------------------------------------------------------------------------
template <unsigned Dim>
struct UniformGridFieldLocalizationHint
{
  UniformGridFieldLocalizationHint(
    Kokkos::View<LO*, DeviceMemorySpace> cell_indices,
    Kokkos::View<Real**, DeviceMemorySpace> coordinates, OutOfBoundsMode mode,
    Kokkos::View<bool*, DeviceMemorySpace> is_out_of_bounds,
    size_t num_out_of_bounds)
    : cell_indices_(cell_indices),
      coordinates_(coordinates),
      mode_(mode),
      is_out_of_bounds_(is_out_of_bounds),
      num_out_of_bounds_(num_out_of_bounds)
  {
  }

  Kokkos::View<LO*, DeviceMemorySpace> cell_indices_;
  Kokkos::View<Real**, DeviceMemorySpace> coordinates_;
  OutOfBoundsMode mode_;
  Kokkos::View<bool*, DeviceMemorySpace> is_out_of_bounds_;
  size_t num_out_of_bounds_;
};

// UniformGridPointEvaluator<Dim> implements PointEvaluator<Real> for structured
// uniform grids. Localization results (cell indices, parametric coordinates)
// are computed once at construction and cached for repeated Evaluate calls.
//
// Evaluation logic is taken directly from UniformGridField::Evaluate.
// Output shape: [num_query_points][num_components].
// Preconditions:
//   - field is compatible with layout_.
//   - values has extents [num_query_points][num_components].
template <unsigned Dim = 2,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class UniformGridPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  UniformGridPointEvaluator(
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout,
    UniformGridFieldLocalizationHint<Dim> hint, Real fill_value)
    : layout_(std::move(layout)),
      grid_(layout_->GetGrid()),
      hint_(std::move(hint)),
      fill_value_(fill_value)
  {
  }

  void Evaluate(
    const Field<Real>& field,
    Rank2View<Real, DeviceMemorySpace, LayoutPolicy> values) const override
  {
    PCMS_FUNCTION_TIMER;
    auto dof_data = field.GetDOFHolderData();
    LO num_points = static_cast<LO>(hint_.coordinates_.extent(0));
    int n_comp = layout_->GetNumComponents();

    PCMS_ALWAYS_ASSERT(values.extent(0) == static_cast<size_t>(num_points));
    PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(n_comp));
    PCMS_ALWAYS_ASSERT(dof_data.size() ==
                       static_cast<size_t>(layout_->GetNumOwnedDofHolder() *
                                           layout_->GetNumComponents()));

    auto cell_indices = hint_.cell_indices_;
    auto coordinates = hint_.coordinates_;
    auto is_out_of_bounds = hint_.is_out_of_bounds_;

    if (layout_->GetOrder() == 0) {
      OutOfBoundsMode mode = hint_.mode_;
      Real fv = fill_value_;
      Kokkos::parallel_for(
        "evaluate_order0_device",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_points, n_comp}),
        KOKKOS_LAMBDA(const LO i, const int c) {
          if (is_out_of_bounds(i) && mode == OutOfBoundsMode::FILL) {
            values(i, c) = fv;
          } else {
            LO dof_idx = cell_indices(i);
            values(i, c) = dof_data(dof_idx, c);
          }
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
        dimensions_view(d) = cell_divisions[d] + 1;
      });

    RealMatView parametric_coords("parametric_coords", num_points, Dim);
    Kokkos::parallel_for(
      "compute_parametric_coords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_CLASS_LAMBDA(const LO i) {
        auto cell_bbox = grid_copy.GetCellBBOX(cell_indices(i));
        for (unsigned d = 0; d < Dim; ++d) {
          Real coord = coordinates(i, d);
          Real cell_min = cell_bbox.center[d] - cell_bbox.half_width[d];
          Real cell_max = cell_bbox.center[d] + cell_bbox.half_width[d];
          parametric_coords(i, d) = (coord - cell_min) / (cell_max - cell_min);
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

    // n_comp == 1 for now (multi-component path would need per-component calls)
    PCMS_ALWAYS_ASSERT(
      n_comp == 1 &&
      "UniformGridPointEvaluator: multi-component order-1 not yet supported");

    RealVecView values_interp("values_interp",
                              static_cast<size_t>(dof_data.size()));
    Kokkos::parallel_for(
      "copy_dof_data_to_values_interp",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
        0, static_cast<LO>(dof_data.size())),
      KOKKOS_LAMBDA(const LO i) { values_interp(i) = dof_data(i, 0); });

    auto interpolator = RegularGridInterpolator(
      parametric_coords, values_interp, cell_indices_interp, dimensions_view);
    auto interpolated = interpolator.linear_interpolation();

    OutOfBoundsMode mode = hint_.mode_;
    Real fv = fill_value_;
    Kokkos::parallel_for(
      "evaluate_order1_device",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_LAMBDA(const LO i) {
        if (is_out_of_bounds(i) && mode == OutOfBoundsMode::FILL) {
          values(i, 0) = fv;
        } else {
          values(i, 0) = interpolated(i);
        }
      });
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
  UniformGrid<Dim> grid_;
  UniformGridFieldLocalizationHint<Dim> hint_;
  Real fill_value_;
};

// UniformGridEvaluatorFactory<Dim> implements FieldEvaluatorFactory<Real> for
// structured uniform grids. It owns the layout and creates
// UniformGridPointEvaluator instances on demand.
//
// Localization logic is taken directly from
// UniformGridField::GetLocalizationHint.
template <unsigned Dim = 2>
class UniformGridEvaluatorFactory : public FieldEvaluatorFactory<Real>
{
public:
  explicit UniformGridEvaluatorFactory(
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout)
    : layout_(std::move(layout)), grid_(layout_->GetGrid())
  {
  }

  const FieldLayout& GetLayout() const override { return *layout_; }

  CoordinateSystem GetCoordinateSystem() const override
  {
    return layout_->GetDOFHolderCoordinates().GetCoordinateSystem();
  }

  bool HasDOFHolderCoordinates() const override { return true; }

  bool SupportsNearestBoundary() const override { return true; }

  std::unique_ptr<PointEvaluator<Real>> CreatePointEvaluator(
    const EvaluationRequest& request) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto coords = request.coords;
    const auto policy = request.policy;
    PCMS_FUNCTION_TIMER;
    if (coords.GetCoordinateSystem() != GetCoordinateSystem()) {
      throw pcms_error(
        "UniformGridEvaluatorFactory: coordinate system mismatch");
    }
    if (policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY &&
        !SupportsNearestBoundary()) {
      throw pcms_error(
        "UniformGridEvaluatorFactory: nearest-boundary evaluation is not "
        "supported");
    }
    if (layout_->GetOrder() == 1 && layout_->GetNumComponents() != 1) {
      throw pcms_error(
        "UniformGridEvaluatorFactory: order-1 multi-component evaluation is "
        "not implemented");
    }

    auto coordinates = coords.GetValues();
    LO num_points = static_cast<LO>(coordinates.extent(0));

    Kokkos::View<LO*, DeviceMemorySpace> cell_indices_device("cell_indices",
                                                             num_points);
    Kokkos::View<Real**, DeviceMemorySpace> coordinates_device(
      "coordinates_device", num_points, Dim);
    Kokkos::parallel_for(
      "copy_coordinates_to_device", num_points, KOKKOS_LAMBDA(const LO i) {
        for (unsigned d = 0; d < Dim; ++d) {
          coordinates_device(i, d) = coordinates(i, d);
        }
      });

    Kokkos::View<bool*, DeviceMemorySpace> is_out_of_bounds_device(
      "is_out_of_bounds", num_points);

    // Parallel reduction to compute cell indices and count out-of-bounds points
    // Create local copy to avoid capturing reference in device lambdas in some
    // compilers
    auto grid_copy = grid_;
    size_t num_out_of_bounds = 0;
    Kokkos::parallel_reduce(
      "localize_points_on_device",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
      KOKKOS_CLASS_LAMBDA(const LO i, size_t& local_count) {
        Omega_h::Vector<Dim> point;
        for (unsigned d = 0; d < Dim; ++d) {
          point[d] = coordinates_device(i, d);
        }

        bool out_of_bounds = !grid_copy.IsPointInBounds(point);
        is_out_of_bounds_device(i) = out_of_bounds;
        if (out_of_bounds) {
          local_count += 1;
        }

        cell_indices_device(i) = grid_copy.ClosestCellID(point);
      },
      num_out_of_bounds);

    // Check for errors after kernel execution
    if (num_out_of_bounds > 0 && policy.mode == OutOfBoundsMode::ERROR) {
      throw pcms_error(
        "UniformGridEvaluatorFactory: " + std::to_string(num_out_of_bounds) +
        " point(s) found outside uniform grid domain");
    }

    UniformGridFieldLocalizationHint<Dim> hint(
      cell_indices_device, coordinates_device, policy.mode,
      is_out_of_bounds_device, num_out_of_bounds);

    return std::make_unique<UniformGridPointEvaluator<Dim>>(
      layout_, std::move(hint), policy.fill_value);
  }

  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return layout_->GetDOFHolderCoordinates();
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
  UniformGrid<Dim> grid_;
};

using UniformGridEvaluatorFactory2D = UniformGridEvaluatorFactory<2>;

} // namespace pcms

#endif // PCMS_UNIFORM_GRID_EVALUATOR_FACTORY_H
