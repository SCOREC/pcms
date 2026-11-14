#ifndef PCMS_TRANSFER_MESH_INTERSECTION_HPP
#define PCMS_TRANSFER_MESH_INTERSECTION_HPP

#include <pcms/configuration.h>
#include <pcms/localization/overlap.h>
#include <pcms/localization/point_search.h>
#include <pcms/localization/queue_visited.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_int_scan.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms
{
// Gather the vertex coordinates of a simplex element (triangle for Dim==2,
// tetrahedron for Dim==3) into an r3d simplex, ready for
// r3d::intersect_simplices.
template <int Dim>
[[nodiscard]] OMEGA_H_INLINE r3d::Few<r3d::Vector<Dim>, Dim + 1>
get_vert_coords_of_elem(const Omega_h::Reals& coords,
                        const Omega_h::LOs& elems2nodes, const int id)
{
  const auto elm_verts = Omega_h::gather_verts<Dim + 1>(elems2nodes, id);

  const Omega_h::Matrix<Dim, Dim + 1> elm_vert_coords =
    Omega_h::gather_vectors<Dim + 1, Dim>(coords, elm_verts);

  r3d::Few<r3d::Vector<Dim>, Dim + 1> r3d_vector;
  for (int i = 0; i < Dim + 1; ++i) {
    for (int d = 0; d < Dim; ++d) {
      r3d_vector[i][d] = elm_vert_coords[i][d];
    }
  }

  return r3d_vector;
}

/**
 * @brief Conservative test for a zero-volume overlap between two simplices.
 *
 * If verticies of one simplex lies in the closed outer half-space of some face
 * plane of the other, their intersection is contained in that plane and so has
 * no volume.
 *
 * @return true when the overlap is provably degenerate and the clip can be
 * skipped
 */
template <int Dim>
[[nodiscard]] OMEGA_H_INLINE bool simplices_have_degenerate_overlap(
  const r3d::Few<r3d::Vector<Dim>, Dim + 1>& a,
  const r3d::Few<r3d::Vector<Dim>, Dim + 1>& b)
{
  return Overlap<Dim, SimplexVsSimplex>::apply(a, b);
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
   * @param start_elements For each target element, the source element
   * containing its centroid (from a point search over the source mesh)
   * @param tgt2src_offsets Offsets array (only used when writing indices).
   * @param[out] nIntersections Number of intersecting source elements per
   * target element.
   * @param[out] tgt2src_indices Indices of intersecting source elements.
   * @param is_count_only If true, only counts intersections; if false, also
   * fills tgt2src_indices.
   * @param use_prefilter if true, skip degenerate clips
   *
   * @note Templated on spatial dimension `Dim`: linear triangles (Dim==2) or
   * linear tetrahedra (Dim==3), using `r3d::intersect_simplices` for geometric
   * intersection.
   *
   * @see r3d::intersect_simplices, intersectTargets
   */
  template <int Dim>
  void adjBasedIntersectSearch(const Kokkos::View<const LO*>& start_elements,
                               const Omega_h::LOs& tgt2src_offsets,
                               Omega_h::Write<Omega_h::LO>& nIntersections,
                               Omega_h::Write<Omega_h::LO>& tgt2src_indices,
                               bool is_count_only, bool use_prefilter = true);
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
                                     Omega_h::Mesh& target_mesh,
                                     bool use_prefilter = true);

/**
 */
IntersectionResults intersectTargets(
  Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
  const GridPointSearchVariant& source_search, bool use_prefilter = true);
} // namespace pcms
#endif // PCMS_TRANSFER_MESH_INTERSECTION_HPP
