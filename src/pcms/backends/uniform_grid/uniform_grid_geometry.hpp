#ifndef PCMS_BACKENDS_UNIFORM_GRID_UNIFORM_GRID_GEOMETRY_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_UNIFORM_GRID_GEOMETRY_HPP

// Private: the uniform-grid geometry class and its binding.

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/backends/uniform_grid/geometry.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_search.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>
#include <utility>

namespace pcms::detail
{

/// Fills, per query point, whether it is found (inside the grid, or any point
/// when `nearest`), its cell and its position in the cell from the search
/// results. (A free function: CUDA does not allow device lambdas in
/// constructors or non-public member functions.)
template <unsigned Dim>
void FillUniformGridBinding(
  Kokkos::View<const typename UniformGridPointSearch<Dim>::Result*> results,
  bool nearest, Kokkos::View<bool*, DeviceMemorySpace> found,
  Kokkos::View<LO*, DeviceMemorySpace> cell_indices,
  Kokkos::View<Real**, DeviceMemorySpace> parametric)
{
  using Result = typename UniformGridPointSearch<Dim>::Result;
  pcms::detail::ForEach(
    "UniformGridBinding", static_cast<LO>(results.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      const Result result = results(i);
      found(i) = nearest || !result.IsOutside();
      cell_indices(i) = result.CellId();
      for (unsigned d = 0; d < Dim; ++d) {
        parametric(i, d) = result.parametric[d];
      }
    });
}

/// Uniform-grid BoundPoints: the cell of each query point and its position in
/// the cell. Points outside the grid are not found, except with
/// OutOfBoundsMode::NEAREST_BOUNDARY, where they use the closest cell and
/// their position relative to it.
template <unsigned Dim>
class UniformGridBinding final : public BoundPoints
{
public:
  UniformGridBinding(std::shared_ptr<const Geometry> geometry,
                     const BindRequest& request, const UniformGrid<Dim>& grid)
    : BoundPoints(std::move(geometry), request), search_(grid)
  {
    Locate();
  }

  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> CellIndices()
    const noexcept
  {
    return Kokkos::subview(cell_indices_,
                           Kokkos::make_pair(LO{0}, NumPoints()));
  }

  /// [point][axis]: the position of each query point in its cell, each
  /// coordinate in [0, 1] for a point inside the grid
  [[nodiscard]] Kokkos::View<const Real**, DeviceMemorySpace> Parametric()
    const noexcept
  {
    return Kokkos::subview(parametric_, Kokkos::make_pair(LO{0}, NumPoints()),
                           Kokkos::ALL);
  }

protected:
  void RebindImpl() override { Locate(); }

  [[nodiscard]] bool UnsharedImpl() const noexcept override
  {
    return cell_indices_.use_count() <= 1 && parametric_.use_count() <= 1;
  }

private:
  /// Locates the query points with the search into cell_indices_ and
  /// parametric_ and sets the found points.
  void Locate()
  {
    const auto results = search_(LocatedPointValues());
    const auto n = static_cast<LO>(results.extent(0));
    const bool nearest = Policy().mode == OutOfBoundsMode::NEAREST_BOUNDARY;
    Kokkos::View<bool*, DeviceMemorySpace> found(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "uniform_grid_found"), n);
    pcms::detail::ReserveRows(cell_indices_, "uniform_grid_cell_indices", n);
    pcms::detail::ReserveRows(parametric_, "uniform_grid_parametric", n, Dim);
    FillUniformGridBinding<Dim>(results, nearest, found, cell_indices_,
                                parametric_);
    SetFound(found);
  }

  UniformGridPointSearch<Dim> search_;
  Kokkos::View<LO*, DeviceMemorySpace> cell_indices_;
  Kokkos::View<Real**, DeviceMemorySpace> parametric_;
};

} // namespace pcms::detail

namespace pcms::uniform_grid
{

/// Geometry of a uniform grid: an affine map of the unit cells, located in
/// closed form (no search to build). Create it with CreateGeometry.
template <unsigned Dim>
class UniformGridGeometry final : public Geometry
{
public:
  /// Own it by std::shared_ptr, as CreateGeometry does: Bind uses
  /// shared_from_this().
  UniformGridGeometry(UniformGrid<Dim> grid,
                      std::shared_ptr<const CoordinateSystem> system)
    : grid_(grid),
      system_(pcms::detail::ResolveCoordinateSystem(std::move(system), Dim))
  {
  }

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return system_;
  }

  [[nodiscard]] bool SupportsNearestBoundary() const noexcept override
  {
    return true;
  }

  [[nodiscard]] const UniformGrid<Dim>& GetGrid() const noexcept
  {
    return grid_;
  }

  /// @brief The tensor Gauss-Legendre rule of the cells, exact for
  ///        polynomials of degree `order` in each coordinate, before binding.
  /// @return the points [point][axis] and their weights in the coordinate
  ///         measure
  /// @throws pcms_error for an order below 1
  [[nodiscard]] std::pair<Kokkos::View<Real**, DeviceMemorySpace>,
                          Kokkos::View<Real*, DeviceMemorySpace>>
  CellRule(int order) const;

protected:
  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

  /// @brief The Gauss rule (CellRule), with the coordinate system's volume
  ///        factor in the weights if the coordinate system has a metric.
  [[nodiscard]] Quadrature CreateQuadrature(
    const pcms::detail::QuadratureRule& rule) const override;

private:
  UniformGrid<Dim> grid_;
  std::shared_ptr<const CoordinateSystem> system_;
};

// instantiated in geometry.cpp
extern template class UniformGridGeometry<1>;
extern template class UniformGridGeometry<2>;
extern template class UniformGridGeometry<3>;

/// @brief Calls `f` with `geometry` as its typed uniform-grid geometry.
/// @return the result of `f`, or an empty one if `geometry` is not a uniform
///         grid
template <typename F>
[[nodiscard]] auto VisitUniformGrid(
  const std::shared_ptr<const Geometry>& geometry,
  F f) -> decltype(f(std::shared_ptr<const UniformGridGeometry<1>>{}))
{
  switch (geometry->GetCoordinateSystem()->Dimension()) {
    case 1:
      if (auto grid =
            std::dynamic_pointer_cast<const UniformGridGeometry<1>>(geometry)) {
        return f(std::move(grid));
      }
      break;
    case 2:
      if (auto grid =
            std::dynamic_pointer_cast<const UniformGridGeometry<2>>(geometry)) {
        return f(std::move(grid));
      }
      break;
    case 3:
      if (auto grid =
            std::dynamic_pointer_cast<const UniformGridGeometry<3>>(geometry)) {
        return f(std::move(grid));
      }
      break;
    default: break;
  }
  return {};
}

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_UNIFORM_GRID_GEOMETRY_HPP
