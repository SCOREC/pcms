#ifndef PCMS_TRANSFER_MESH_INTERSECTION_HPP
#define PCMS_TRANSFER_MESH_INTERSECTION_HPP

#include <pcms/localization/point_search.h>
#include <pcms/transfer/queue_visited.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_int_scan.hpp>
#include <r3d.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms
{
constexpr static double abs_tol = 1e-18; /// abs tolerance
constexpr static double rel_tol = 1e-12; /// rel tolerance

[[nodiscard]] OMEGA_H_INLINE r3d::Few<r3d::Vector<2>, 3>
get_vert_coords_of_elem(const Omega_h::Reals& coords,
                        const Omega_h::LOs& faces2nodes, const int id)
{
  const auto elm_verts = Omega_h::gather_verts<3>(faces2nodes, id);

  const Omega_h::Matrix<2, 3> elm_vert_coords =
    Omega_h::gather_vectors<3, 2>(coords, elm_verts);

  r3d::Few<r3d::Vector<2>, 3> r3d_vector;
  for (int i = 0; i < 3; ++i) {
    r3d_vector[i][0] = elm_vert_coords[i][0];
    r3d_vector[i][1] = elm_vert_coords[i][1];
  }

  return r3d_vector;
}

/**
 * @brief Stores results of mesh element intersections for conservative
 * transfer.
 *
 * Contains mappings from each target element to the list of source elements
 * that intersect with it. Used to guide integration over overlapping regions.
 *
 * - `tgt2src_offsets[i]` is the offset into `tgt2src_indices` where source
 *    elements for target element `i` begin.
 * - `tgt2src_indices` contains flattened indices of source elements per target.
 */
struct IntersectionResults
{
  Omega_h::LOs tgt2src_offsets;
  Omega_h::LOs tgt2src_indices;
};

class FindIntersections
{
private:
  Omega_h::Mesh& source_mesh_;
  Omega_h::Mesh& target_mesh_;

public:
  FindIntersections(Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh)
    : source_mesh_(source_mesh), target_mesh_(target_mesh)
  {
  }

  /**
   * @brief Performs adjacency-based intersection search between target and
   * source elements.
   *
   * For each target element, starting from the source element that contains its
   * centroid, a queue-based BFS traversal is used over the adjacency graph of
   * source elements. If an element intersects the target triangle (based on
   * area tolerance), it is included.
   *
   * @param tgt2src_offsets Offsets array (only used when writing indices).
   * @param[out] nIntersections Number of intersecting source elements per
   * target element.
   * @param[out] tgt2src_indices Indices of intersecting source elements.
   * @param is_count_only If true, only counts intersections; if false, also
   * fills tgt2src_indices.
   *
   * @note This method assumes 2D linear triangles and uses
   * `r3d::intersect_simplices` for geometric intersection.
   *
   * @see r3d::intersect_simplices, intersectTargets
   */
  void adjBasedIntersectSearch(const Omega_h::LOs& tgt2src_offsets,
                               Omega_h::Write<Omega_h::LO>& nIntersections,
                               Omega_h::Write<Omega_h::LO>& tgt2src_indices,
                               bool is_count_only);
};

/**
 * @brief Computes source-target element intersections for conservative
 * projection.
 *
 * For each target element in the target mesh, this function identifies source
 * elements from the source mesh that geometrically intersect with it using an
 * adjacency-based breadth-first search strategy. The result is returned as a
 * compact mapping.
 *
 * @param source_mesh The source Omega_h mesh.
 * @param target_mesh The target Omega_h mesh.
 * @return An IntersectionResults struct containing target-to-source mapping
 * data.
 *
 * @note The intersection test is done using 2D polygon intersection routines
 * from r3d. Only valid (non-degenerate) polygonal intersections are included.
 *
 * @see FindIntersections::adjBasedIntersectSearch
 */

IntersectionResults intersectTargets(Omega_h::Mesh& source_mesh,
                                     Omega_h::Mesh& target_mesh);
} // namespace pcms
#endif // PCMS_TRANSFER_MESH_INTERSECTION_HPP
