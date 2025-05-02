#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <unordered_map>
#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include "types.h"
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>
#include "pcms/uniform_grid.h"
#include "pcms/bounding_box.h"

namespace pcms
{
//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail {

/**
 * Check if a point is within a radial cutoff distance of an axis-aligned bounding box (AABB).
 *
 * Used to quickly reject boxes that are too far from the point during neighbor search.
 */
template <unsigned dim>
KOKKOS_INLINE_FUNCTION
bool radial_intersects_bbox(const double pt[dim], const AABBox<dim>& bbox, double cutoff_squared)
{
  double d_min = 0.0; // Accumulates squared distance from point to box (only in directions where point is outside the box)
  for (unsigned d = 0; d < dim; ++d) {
    double dist = fabs(pt[d] - bbox.center[d]); // Distance from point to box center along axis d
    double excess = dist - bbox.half_width[d]; // How far point lies outside the box in this axis
    if (excess > 0.0) {
      d_min += excess * excess;  // Add squared excess if point is outside box in this axis
      if (d_min > cutoff_squared) return false; 
    }
  }
  return true; 
}

template <unsigned dim>
struct GridRadialNeighborFunctor {
    Kokkos::View<const double**> target_points;
    Kokkos::View<const double**> source_points;
    Kokkos::View<const pcms::UniformGrid<dim>[1]> grid;
    Kokkos::View<const LO*> cell_ptrs;
    Kokkos::View<const LO*> cell_indices;
    Kokkos::View<const double[dim]> cell_size;
    double cutoff;
    double cutoff_squared;
    LO num_cells;

    KOKKOS_FUNCTION
    GridRadialNeighborFunctor(
        Kokkos::View<const double**> tgt_pts,
        Kokkos::View<const double**> src_pts,
        Kokkos::View<const pcms::UniformGrid<dim>[1]> grid_in,
        Kokkos::View<const LO*> cell_ptrs_in,
        Kokkos::View<const LO*> cell_indices_in,
        double cutoff_in,
        LO num_cells_in,
        Kokkos::View<const double[dim]> cell_size_in)
        : target_points(tgt_pts),
          source_points(src_pts),
          grid(grid_in),
          cell_ptrs(cell_ptrs_in),
          cell_indices(cell_indices_in),
          cell_size(cell_size_in),
          cutoff(cutoff_in),
          cutoff_squared(cutoff_in * cutoff_in),
          num_cells(num_cells_in) {}

    KOKKOS_INLINE_FUNCTION
    LO operator()(LO target_idx, LO* fill) const {
        double pt[dim];
        for (int d = 0; d < dim; ++d)
            pt[d] = target_points(target_idx, d);

        LO count = 0;
        const auto& grid_obj = grid(0);

        // Compute min/max grid indices around the cutoff
        int min_idx[dim], max_idx[dim];
        for (int d = 0; d < dim; ++d) {
            const double min_coord = pt[d] - cutoff;
            const double max_coord = pt[d] + cutoff;

            min_idx[d] = static_cast<int>((min_coord - grid_obj.bot_left[d]) / cell_size(d));
            min_idx[d] = (min_idx[d] < 0) ? 0 : (min_idx[d] >= grid_obj.divisions[d]) ? grid_obj.divisions[d] - 1 : min_idx[d];

            max_idx[d] = static_cast<int>((max_coord - grid_obj.bot_left[d]) / cell_size(d));
            max_idx[d] = (max_idx[d] < 0) ? 0 : (max_idx[d] >= grid_obj.divisions[d]) ? grid_obj.divisions[d] - 1 : max_idx[d];
        }

        // Iterate over intersecting grid cells
        if constexpr (dim == 3) {
            for (int z = min_idx[2]; z <= max_idx[2]; ++z)
                for (int y = min_idx[1]; y <= max_idx[1]; ++y)
                    for (int x = min_idx[0]; x <= max_idx[0]; ++x)
                        process_cell(x, y, z, pt, count, fill);
        } else if constexpr (dim == 2) {
            for (int y = min_idx[1]; y <= max_idx[1]; ++y)
                for (int x = min_idx[0]; x <= max_idx[0]; ++x)
                    process_cell(x, y, 0, pt, count, fill);
        } else {
            for (int x = min_idx[0]; x <= max_idx[0]; ++x)
                process_cell(x, 0, 0, pt, count, fill);
        }

        return count;
    }

private:
    KOKKOS_INLINE_FUNCTION
    void process_cell(int x, int y, int z, const double pt[dim], LO& count, LO* fill) const {
        const auto& grid_obj = grid(0);
        LO cell_id;

        if constexpr (dim == 3) {
            cell_id = z * (grid_obj.divisions[0] * grid_obj.divisions[1]) + y * grid_obj.divisions[0] + x;
        } else if constexpr (dim == 2) {
            cell_id = y * grid_obj.divisions[0] + x;
        } else {
            cell_id = x;
        }

        if (cell_id >= num_cells || cell_ptrs[cell_id] == cell_ptrs[cell_id + 1]) return;

        const auto bbox = grid_obj.GetCellBBOX(cell_id);
        if (!radial_intersects_bbox<dim>(pt, bbox, cutoff_squared)) return;

        for (LO i = cell_ptrs[cell_id]; i < cell_ptrs[cell_id + 1]; ++i) {
            const LO src_idx = cell_indices[i];
            double r2 = 0.0;
            for (int d = 0; d < dim; ++d) {
                double diff = pt[d] - source_points(src_idx, d);
                r2 += diff * diff;
                if (r2 > cutoff_squared) break;
            }
            if (r2 <= cutoff_squared) {
                if (fill) fill[count] = src_idx;
                ++count;
            }
        }
    }
};

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh, Kokkos::View<Uniform2DGrid[1]> grid, int num_grid_cells);
}
KOKKOS_FUNCTION
Omega_h::Vector<3> barycentric_from_global(
  const Omega_h::Vector<2>& point, const Omega_h::Matrix<2, 3>& vertex_coords);

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

class GridPointSearch
{
  using CandidateMapT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  static constexpr auto dim = 2;
  struct Result {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2
    };

    Dimensionality dimensionality;
    LO tri_id;
    Omega_h::Vector<dim + 1> parametric_coords;
  };

  GridPointSearch(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  Kokkos::View<Result*> operator()(Kokkos::View<Real*[dim] > point) const;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Kokkos::View<Uniform2DGrid[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

} // namespace detail
#endif // PCMS_COUPLING_POINT_SEARCH_H
