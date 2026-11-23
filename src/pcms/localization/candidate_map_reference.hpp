#ifndef PCMS_LOCALIZATION_CANDIDATE_MAP_REFERENCE_HPP
#define PCMS_LOCALIZATION_CANDIDATE_MAP_REFERENCE_HPP

#include "pcms/localization/candidate_map.hpp"

namespace pcms::detail
{

/**
 * Cell-major brute-force build: for every grid cell, test every shape.
 *
 * This is the implementation the shape-major build replaced. It is
 * O(num_cells * nshapes) and must not be used in production, but it is the
 * ground truth for what the candidate map should contain, so it is kept as an
 * oracle for the equivalence test.
 */
template <typename S>
struct ReferenceCellMajorFunctor
{
  static constexpr int Dim = ShapeDimension<S>;

  ReferenceCellMajorFunctor(Kokkos::View<const S*> shapes,
                            const UniformGrid<Dim>& grid)
    : shapes_(shapes), grid_(grid)
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO row, LO* fill) const
  {
    const auto cell_bbox = grid_.GetCellBBOX(row);
    LO num_intersections = 0;
    for (LO i = 0; i < static_cast<LO>(shapes_.extent(0)); ++i) {
      if (Overlaps(shapes_(i), cell_bbox)) {
        if (fill) {
          fill[num_intersections] = i;
        }
        ++num_intersections;
      }
    }
    return num_intersections;
  }

private:
  Kokkos::View<const S*> shapes_;
  UniformGrid<Dim> grid_;
};

// The candidate map CreateCandidateMap builds (without inflation), by brute
// force. A test oracle; never use in production.
template <Shape S>
CandidateMap<ShapeDimension<S>> CreateCandidateMapReference(
  Kokkos::View<const S*> shapes, const UniformGrid<ShapeDimension<S>>& grid)
{
  CandidateBins bins{};
  Kokkos::count_and_fill_crs(bins, static_cast<LO>(grid.GetNumCells()),
                             ReferenceCellMajorFunctor<S>{shapes, grid});
  return CandidateMap<ShapeDimension<S>>(grid, std::move(bins));
}

} // namespace pcms::detail

#endif // PCMS_LOCALIZATION_CANDIDATE_MAP_REFERENCE_HPP
