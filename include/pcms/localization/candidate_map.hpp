#ifndef PCMS_LOCALIZATION_CANDIDATE_MAP_HPP
#define PCMS_LOCALIZATION_CANDIDATE_MAP_HPP

#include <Kokkos_Core.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

#include "pcms/localization/detail/candidate_map.hpp"
#include "pcms/localization/overlap.hpp"
#include "pcms/localization/shapes.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

namespace pcms
{

/** An axis-aligned bounding box type. */
template <typename T>
concept IsAABBox = detail::IsAABBox<std::remove_cvref_t<T>>;

/**
 * A shape a candidate map can bin: copyable to the device, with a bounding box
 * (`BoundingBox`, found by argument-dependent lookup). Binning also tests it
 * against the cells with `Overlaps`, which reports a missing test at compile
 * time.
 */
template <typename S>
concept Shape = std::is_trivially_copyable_v<S> && requires(const S& s) {
  { BoundingBox(s) } -> IsAABBox;
};

/**
 * The shapes that may overlap each cell of a uniform grid, the bins: the
 * candidates a query tests instead of every shape.
 */
template <int Dim>
class CandidateMap
{
public:
  /** Per row, shape indices, ascending. */
  using Bins = detail::CandidateBins;

  CandidateMap() = default;

  /**
   * \param grid the grid
   * \param bins per cell of `grid`, the indices of the shapes that may overlap
   *        it, ascending
   */
  CandidateMap(const UniformGrid<Dim>& grid, Bins bins)
    : grid_(grid), bins_(std::move(bins))
  {
  }

  [[nodiscard]] const UniformGrid<Dim>& GetGrid() const noexcept
  {
    return grid_;
  }

  /** Per grid cell, the indices of the shapes that may overlap it. */
  [[nodiscard]] const Bins& GetBins() const noexcept { return bins_; }

  /**
   * Per point, the shapes that may contain it: the bin of its cell, or of the
   * nearest cell for a point outside the grid.
   */
  [[nodiscard]] Bins GetCandidates(Kokkos::View<const Real* [Dim]> points) const
  {
    Bins candidates{};
    Kokkos::count_and_fill_crs(
      candidates, static_cast<LO>(points.extent(0)),
      detail::PointCandidatesFunctor<Dim>{points, grid_, bins_});
    return candidates;
  }

  /**
   * Per box, the shapes that may overlap it: those in the bins of the cells
   * it overlaps, each once, ascending.
   */
  [[nodiscard]] Bins GetCandidates(Kokkos::View<const AABBox<Dim>*> boxes) const
  {
    Bins candidates{};
    Kokkos::count_and_fill_crs(
      candidates, static_cast<LO>(boxes.extent(0)),
      detail::BoxCandidatesFunctor<Dim>{boxes, grid_, bins_});
    return detail::UniqueRows(std::move(candidates));
  }

private:
  UniformGrid<Dim> grid_{};
  Bins bins_{};
};

/**
 * Bins `shapes` on `grid`: a shape is in the bin of every cell it overlaps.
 * \param inflation distance by which the shapes and the cells are widened, so
 *        a shape within it of a cell is in the cell's bin (e.g. a search
 *        tolerance)
 */
template <Shape S>
[[nodiscard]] CandidateMap<ShapeDimension<S>> CreateCandidateMap(
  Kokkos::View<const S*> shapes, const UniformGrid<ShapeDimension<S>>& grid,
  Real inflation = 0)
{
  return CandidateMap<ShapeDimension<S>>(
    grid, detail::BuildCandidateBins<S>(shapes, grid, inflation));
}

} // namespace pcms

#endif // PCMS_LOCALIZATION_CANDIDATE_MAP_HPP
