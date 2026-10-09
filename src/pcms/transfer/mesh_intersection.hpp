#ifndef PCMS_TRANSFER_MESH_INTERSECTION_HPP
#define PCMS_TRANSFER_MESH_INTERSECTION_HPP

#include <pcms/configuration.h>
#include <pcms/localization/point_search.h>
#include <pcms/discretization/discretization/omega_h.hpp>
#include <pcms/field/function_space.h>
#include <memory>
#include <pcms/localization/queue_visited.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_int_scan.hpp>
#include <r3d.hpp>
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
 * Uses the face planes of both simplices as candidate separating axes: if
 * every vertex of one simplex lies in the closed outer half-space of some face
 * plane of the other, their intersection is contained in that plane and so has
 * no volume. Rejecting such a pair here avoids an r3d clip, which is an order
 * of magnitude more expensive. Face-adjacent elements -- the common case when
 * walking a source mesh's dual graph -- are separated exactly by their shared
 * face's plane, so they are caught after one or two tests.
 *
 * The face planes alone are an incomplete set of separating axes (a complete
 * test also needs the edge-edge cross products), so a genuinely disjoint pair
 * may not be detected. That is safe: an undetected pair falls through to the
 * exact clip. The test never reports separation for a pair with positive
 * overlap volume, which is what makes it usable as a pre-filter.
 *
 * @return true when the overlap is provably degenerate and the clip can be
 * skipped; false when the exact clip is still required.
 */
template <int Dim>
[[nodiscard]] OMEGA_H_INLINE bool simplices_have_degenerate_overlap(
  const r3d::Few<r3d::Vector<Dim>, Dim + 1>& a,
  const r3d::Few<r3d::Vector<Dim>, Dim + 1>& b)
{
  // Each iteration treats one simplex as the half-space owner. Both are
  // needed: a face of `a` may separate the pair while no face of `b` does.
  for (int side = 0; side < 2; ++side) {
    const auto& owner = (side == 0) ? a : b;
    const auto& other = (side == 0) ? b : a;

    // Face `f` is the face opposite vertex `f`.
    for (int f = 0; f < Dim + 1; ++f) {
      // Vertices spanning the face, in index order with `f` skipped.
      int face_verts[Dim];
      int n_face_verts = 0;
      for (int v = 0; v < Dim + 1; ++v) {
        if (v != f) {
          face_verts[n_face_verts++] = v;
        }
      }
      const auto& origin = owner[face_verts[0]];

      // Face normal, oriented to point away from the opposite vertex so the
      // owner's interior lies on the negative side.
      r3d::Vector<Dim> normal;
      if constexpr (Dim == 2) {
        const auto edge = owner[face_verts[1]] - origin;
        normal[0] = edge[1];
        normal[1] = -edge[0];
      } else {
        const auto e1 = owner[face_verts[1]] - origin;
        const auto e2 = owner[face_verts[2]] - origin;
        normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
        normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
        normal[2] = e1[0] * e2[1] - e1[1] * e2[0];
      }

      Omega_h::Real opposite_side = 0.0;
      Omega_h::Real normal_sq = 0.0;
      for (int d = 0; d < Dim; ++d) {
        opposite_side += normal[d] * (owner[f][d] - origin[d]);
        normal_sq += normal[d] * normal[d];
      }
      // A degenerate face gives no usable axis.
      if (normal_sq == 0.0) {
        continue;
      }
      const Omega_h::Real orientation = (opposite_side > 0.0) ? -1.0 : 1.0;

      bool all_outside = true;
      for (int v = 0; v < Dim + 1 && all_outside; ++v) {
        Omega_h::Real signed_distance = 0.0;
        for (int d = 0; d < Dim; ++d) {
          signed_distance += normal[d] * (other[v][d] - origin[d]);
        }
        all_outside = (orientation * signed_distance) >= 0.0;
      }
      if (all_outside) {
        return true;
      }
    }
  }
  return false;
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
   * containing its centroid (from a point search over the source mesh); the
   * BFS starts there. Located once by the caller and shared by both passes.
   * @param tgt2src_offsets Offsets array (only used when writing indices).
   * @param[out] nIntersections Number of intersecting source elements per
   * target element.
   * @param[out] tgt2src_indices Indices of intersecting source elements.
   * @param is_count_only If true, only counts intersections; if false, also
   * fills tgt2src_indices.
   * @param use_prefilter If true, skip the exact clip for pairs that
   * simplices_have_degenerate_overlap rejects. Only exposed so tests can
   * compare against the unfiltered path; production callers want it on.
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
 * @brief As above, locating target centroids with a caller-owned search over
 * the source mesh instead of building one.
 *
 * The source function space already owns such a search for point evaluation,
 * and constructing another costs a full candidate-map build and holds a second
 * copy of it for the operator's lifetime. The resulting map does not depend on
 * which search located the centroids.
 *
 * @param source_search Must have been built over `source_mesh`; that is the
 * caller's responsibility. Its dimension must match the meshes (pcms_error
 * otherwise).
 */
IntersectionResults intersectTargets(
  Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
  const GridPointSearchVariant& source_search, bool use_prefilter = true);
/// Overlap map between a source and a target discretization: for each target
/// element, the source elements it intersects. Built once per mesh pair and
/// shared by every operator on that pair.
class MeshIntersection
{
public:
  virtual ~MeshIntersection() = default;
  virtual std::shared_ptr<const Discretization> GetSourceDiscretization()
    const noexcept = 0;
  virtual std::shared_ptr<const Discretization> GetTargetDiscretization()
    const noexcept = 0;
  /// Target-to-source CRS map; offsets are sized to the target element count.
  virtual const IntersectionResults& GetTargetToSource() const noexcept = 0;
};

/// MeshIntersection of two Omega_h simplex meshes of the same dimension.
class OmegaHMeshIntersection final : public MeshIntersection
{
public:
  /// Reuses the source space's point search when the space owns one.
  OmegaHMeshIntersection(const FunctionSpace& source_space,
                         const FunctionSpace& target_space,
                         bool use_prefilter = true);
  /// @param source_search search over the source mesh to reuse; may be null
  OmegaHMeshIntersection(std::shared_ptr<const OmegaHDiscretization> source,
                         std::shared_ptr<const OmegaHDiscretization> target,
                         const GridPointSearchVariant* source_search = nullptr,
                         bool use_prefilter = true);

  std::shared_ptr<const Discretization> GetSourceDiscretization()
    const noexcept override;
  std::shared_ptr<const Discretization> GetTargetDiscretization()
    const noexcept override;
  const IntersectionResults& GetTargetToSource() const noexcept override;
  Omega_h::Mesh& GetSourceMesh() const noexcept;
  Omega_h::Mesh& GetTargetMesh() const noexcept;

private:
  std::shared_ptr<const OmegaHDiscretization> source_;
  std::shared_ptr<const OmegaHDiscretization> target_;
  IntersectionResults results_;
};

/// The point search owned by an Omega_h Lagrange space; null for other spaces.
const GridPointSearchVariant* SourceSearchFromSpace(const FunctionSpace& space);

/// MeshIntersection for the spaces' discretizations, dispatched on their type.
/// Throws pcms_error when no implementation covers the pair.
std::shared_ptr<MeshIntersection> IntersectMeshes(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  bool use_prefilter = true);

} // namespace pcms
#endif // PCMS_TRANSFER_MESH_INTERSECTION_HPP
