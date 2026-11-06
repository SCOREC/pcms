#ifndef PCMS_FIELD_EVALUATOR_UNIFORM_GRID_GEOMETRY_HPP
#define PCMS_FIELD_EVALUATOR_UNIFORM_GRID_GEOMETRY_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/localization/uniform_grid_search.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <Kokkos_Core.hpp>
#include <Omega_h_vector.hpp>
#include <memory>
#include <string>
#include <utility>

namespace pcms::backend::uniform_grid
{

namespace detail
{

/// Fills, per query point, whether it is found (inside the grid, or any point
/// when `nearest`), its cell and its position in the cell from the search
/// results. (A free function: CUDA does not allow device lambdas in
/// constructors or non-public member functions.)
template <unsigned Dim>
void FillUniformGridBinding(
  Kokkos::View<
    const typename localization::UniformGridPointSearch<Dim>::Result*>
    results,
  bool nearest, Kokkos::View<bool*, DeviceMemorySpace> found,
  Kokkos::View<LO*, DeviceMemorySpace> cell_indices,
  Kokkos::View<Real**, DeviceMemorySpace> parametric)
{
  using Result = typename localization::UniformGridPointSearch<Dim>::Result;
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

} // namespace detail

/// Uniform-grid bound source: the cell of each query point and its position in
/// the cell. Points outside the grid are not found, except with
/// OutOfBoundsMode::NEAREST_BOUNDARY, where they use the closest cell and
/// their position relative to it.
template <unsigned Dim>
class UniformGridBinding final : public BoundSource
{
public:
  UniformGridBinding(std::shared_ptr<const Geometry> geometry,
                     const BindRequest& request, const UniformGrid<Dim>& grid)
    : BoundSource(std::move(geometry), request), search_(grid)
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
    detail::FillUniformGridBinding<Dim>(results, nearest, found, cell_indices_,
                                        parametric_);
    SetFound(found);
  }

  localization::UniformGridPointSearch<Dim> search_;
  Kokkos::View<LO*, DeviceMemorySpace> cell_indices_;
  Kokkos::View<Real**, DeviceMemorySpace> parametric_;
};

/// Geometry of a uniform grid: an affine map of the unit cells, located in
/// closed form (no search to build). Create it with CreateUniformGridGeometry
/// (pcms/geometry/geometry.hpp).
template <unsigned Dim>
class UniformGridGeometry final : public Geometry
{
public:
  /// Own it by std::shared_ptr, as CreateUniformGridGeometry does: Bind uses
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

protected:
  [[nodiscard]] std::shared_ptr<BoundSource> BindImpl(
    const BindRequest& request) const override
  {
    PCMS_FUNCTION_TIMER;
    return std::make_shared<UniformGridBinding<Dim>>(shared_from_this(),
                                                     request, grid_);
  }

private:
  UniformGrid<Dim> grid_;
  std::shared_ptr<const CoordinateSystem> system_;
};

} // namespace pcms::backend::uniform_grid

#endif // PCMS_FIELD_EVALUATOR_UNIFORM_GRID_GEOMETRY_HPP
