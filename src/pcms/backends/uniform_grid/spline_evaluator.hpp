#ifndef PCMS_BACKENDS_UNIFORM_GRID_SPLINE_EVALUATOR_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_SPLINE_EVALUATOR_HPP

#include "pcms/backends/uniform_grid/spline_interpolator.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/types.hpp"

#include <array>
#include <memory>

namespace pcms::uniform_grid
{

namespace detail
{

struct InitCoordsFunctor
{
  Kokkos::View<Real*, DeviceMemorySpace> coords;
  Real bot_left;
  Real delta;

  KOKKOS_FUNCTION void operator()(LO i) const
  {
    coords(i) = bot_left + delta * i;
  }
};

/// @brief Copies one coordinate of the listed points.
/// @param points coordinates, one row per point
/// @param rows indices of the points to copy
/// @param axis the coordinate to copy
/// @param out the coordinates of the listed points
inline void GatherPoints(Kokkos::View<const Real**, DeviceMemorySpace> points,
                         Kokkos::View<const LO*, DeviceMemorySpace> rows,
                         int axis, Kokkos::View<Real*, DeviceMemorySpace> out)
{
  Kokkos::parallel_for(
    "gather_spline_points",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(rows.extent(0))),
    KOKKOS_LAMBDA(const LO k) { out(k) = points(rows(k), axis); });
}

} // namespace detail

/// @brief Evaluates the not-a-knot cubic (1D) or bicubic (2D) spline of a
///        field's grid values, or its gradient for
///        DifferentialOperator::Gradient.
/// @tparam Dim grid dimension, 1 or 2
template <unsigned Dim,
          typename LayoutPolicy =
            pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class UniformGridSplinePointEvaluator
  : public PointEvaluator<Real, LayoutPolicy>
{
  using Vec = Kokkos::View<Real*, DeviceMemorySpace>;
  static constexpr int NumDerivatives = (Dim + 1) * (Dim + 2) / 2;

public:
  /// @param layout order-1 single-component uniform-grid layout
  /// @param binding cells of the query points
  /// @param op Value or Gradient
  UniformGridSplinePointEvaluator(
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout,
    std::shared_ptr<const pcms::detail::UniformGridBinding<Dim>> binding,
    DifferentialOperator op = DifferentialOperator::Value)
    : PointEvaluator<Real, LayoutPolicy>(binding, op, layout),
      grid_(layout->GetGrid()),
      found_(binding->FoundPoints()),
      spline_output_("uniform_grid_spline_output", found_.extent(0),
                     op == DifferentialOperator::Value ? 1 : Dim),
      select_("uniform_grid_spline_select", NumDerivatives)
  {
    for (unsigned d = 0; d < Dim; ++d) {
      x_coords_[d] = Vec("uniform_grid_spline_x", grid_.divisions[d] + 1);
      eval_x_[d] = Vec("uniform_grid_spline_eval_x", found_.extent(0));
      Kokkos::parallel_for(
        "init_x_coords",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
          0, grid_.divisions[d] + 1),
        detail::InitCoordsFunctor{x_coords_[d], grid_.bot_left[d],
                                  grid_.edge_length[d] / grid_.divisions[d]});
      detail::GatherPoints(binding->LocatedPointValues(), found_, d,
                           eval_x_[d]);
    }
    // The selector picks f, then df/dx_d, then the second derivatives.
    auto select = Kokkos::create_mirror_view(select_);
    for (int k = 0; k < NumDerivatives; ++k) {
      select(k) = op == DifferentialOperator::Value
                    ? (k == 0)
                    : (k >= 1 && k <= static_cast<int>(Dim));
    }
    Kokkos::deep_copy(select_, select);
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  // TODO: add passkey or other mechanism to protect this
  void EvaluateImpl(
    const Field<Real>& field,
    ValueView<Real, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    const auto values = out.GetValues();
    const auto dof_data = field.GetDOFHolderData().GetValues();
    const auto found = found_;
    const auto num_in_bounds = static_cast<LO>(found.extent(0));
    if (num_in_bounds == 0) {
      return;
    }
    // x-fastest dof order becomes the interpolators' order
    const LO nx = grid_.divisions[0] + 1;
    const LO ny = Dim == 2 ? grid_.divisions[1] + 1 : 1;
    PCMS_ALWAYS_ASSERT(static_cast<LO>(dof_data.size()) == nx * ny);
    Vec spline_values("uniform_grid_spline_values",
                      static_cast<size_t>(nx * ny));
    Kokkos::parallel_for(
      "copy_spline_values_device",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ny, nx}),
      KOKKOS_LAMBDA(const LO iy, const LO ix) {
        spline_values(ix * ny + iy) = dof_data(iy * nx + ix, 0);
      });

    const auto view = [](const Vec& v) {
      return Rank1View<Real, DeviceMemorySpace>(v.data(), v.extent(0));
    };
    // one column per selected spline output: f, or the gradient
    const auto spline_output = spline_output_;
    const auto num_outputs = static_cast<int>(spline_output.extent(1));
    const Rank2View<Real, DeviceMemorySpace> output(
      spline_output.data(), spline_output.extent(0), num_outputs);
    if constexpr (Dim == 1) {
      ExplicitCubicSplineInterpolator<Real, DeviceMemorySpace> interpolator(
        view(x_coords_[0]), view(spline_values), BoundaryCondition::NOT_A_KNOT,
        0, BoundaryCondition::NOT_A_KNOT, 0);
      interpolator.evaluate(select_, view(eval_x_[0]), output);
    } else {
      const Rank1View<Real, DeviceMemorySpace> empty_bc(nullptr, 0);
      CompactBiCubicSplineInterpolator<Real, DeviceMemorySpace> interpolator(
        view(x_coords_[0]), view(x_coords_[1]), view(spline_values),
        BoundaryCondition::NOT_A_KNOT, empty_bc, BoundaryCondition::NOT_A_KNOT,
        empty_bc, BoundaryCondition::NOT_A_KNOT, empty_bc,
        BoundaryCondition::NOT_A_KNOT, empty_bc);
      interpolator.evaluate(select_, view(eval_x_[0]), view(eval_x_[1]),
                            output);
    }

    Kokkos::parallel_for(
      "copy_in_bounds_values",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_in_bounds),
      KOKKOS_LAMBDA(LO k) {
        for (int c = 0; c < num_outputs; ++c) {
          values(found(k), c) = spline_output(k, c);
        }
      });
  }

private:
  UniformGrid<Dim> grid_;
  Kokkos::View<const LO*, DeviceMemorySpace> found_;
  std::array<Vec, Dim> x_coords_;
  // per bound point: found points' coordinates and spline output
  std::array<Vec, Dim> eval_x_;
  Kokkos::View<Real**, DeviceMemorySpace> spline_output_;
  Kokkos::View<LO*, DeviceMemorySpace> select_;
};

/// @brief Point evaluation of the cubic (1D) or bicubic (2D) spline of fields
///        on a uniform grid, and of its gradient.
/// @tparam Dim grid dimension, 1 or 2
template <unsigned Dim>
class UniformGridSplineEvaluatorFactory
  : public pcms::detail::TypedPointEvaluatorFactory<
      Real, pcms::detail::UniformGridBinding<Dim>>
{
  using Bound = pcms::detail::UniformGridBinding<Dim>;

public:
  /// @param geometry geometry of the layout's grid
  /// @param layout order-1 single-component uniform-grid layout
  /// @throws pcms_error if the layout is not order 1 with one component
  UniformGridSplineEvaluatorFactory(
    std::shared_ptr<const UniformGridGeometry<Dim>> geometry,
    std::shared_ptr<const UniformGridFieldLayout<Dim>> layout)
    : pcms::detail::TypedPointEvaluatorFactory<Real, Bound>(
        std::move(geometry), {DifferentialOperator::Gradient}),
      layout_(std::move(layout))
  {
    if (layout_->GetOrder() != 1) {
      throw pcms_error(
        "UniformGridSplineEvaluatorFactory: spline evaluation requires an "
        "order-1 uniform-grid layout");
    }
    if (layout_->GetNumComponents() != 1) {
      throw pcms_error(
        "UniformGridSplineEvaluatorFactory: spline evaluation only supports "
        "single-component fields");
    }
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const Bound> bound, DifferentialOperator op) const override
  {
    return std::make_unique<UniformGridSplinePointEvaluator<Dim>>(
      layout_, std::move(bound), op);
  }

private:
  std::shared_ptr<const UniformGridFieldLayout<Dim>> layout_;
};

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_SPLINE_EVALUATOR_HPP
