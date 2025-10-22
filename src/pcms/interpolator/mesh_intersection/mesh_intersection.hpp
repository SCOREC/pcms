#ifndef PCMS_INTERPOLATOR_MESH_INTERSECTION_HPP
#define PCMS_INTERPOLATOR_MESH_INTERSECTION_HPP

#include <pcms/point_search.h>
#include <pcms/interpolator/queue_visited.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_int_scan.hpp>
#include <r3d.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

constexpr static double abs_tol = 1e-18; /// abs tolerance
constexpr static double rel_tol = 1e-12; /// rel tolerance

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
    : source_mesh_(source_mesh), target_mesh_(target_mesh) {};

  void adjBasedIntersectSearch(const Omega_h::LOs& tgt2src_offsets,
                               Omega_h::Write<Omega_h::LO>& nIntersections,
                               Omega_h::Write<Omega_h::LO>& tgt2src_indices,
                               bool is_count_only);
};

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

inline Kokkos::View<Omega_h::Real* [2]> compute_centroid(Omega_h::Mesh& mesh)
{

  const auto& mesh_coords = mesh.coords();
  const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& nfaces = mesh.nfaces();

  Kokkos::View<Omega_h::Real* [2]> cell_centroids("centroid coordinates",
                                                  nfaces);
  Omega_h::parallel_for(
    "calculate the centroid in each tri element", nfaces,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      const auto current_el_verts = Omega_h::gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
        Omega_h::gather_vectors<3, 2>(mesh_coords, current_el_verts);
      auto centroid = Omega_h::average(current_el_vert_coords);
      cell_centroids(id, 0) = centroid[0];
      cell_centroids(id, 1) = centroid[1];
    });

  return cell_centroids;
}
/**
 * @brief Performs adjacency-based intersection search between target and source
 * elements.
 *
 * For each target element, starting from the source element that contains its
 * centroid, a queue-based BFS traversal is used over the adjacency graph of
 * source elements. If an element intersects the target triangle (based on area
 * tolerance), it is included.
 *
 * @param tgt2src_offsets Offsets array (only used when writing indices).
 * @param[out] nIntersections Number of intersecting source elements per target
 * element.
 * @param[out] tgt2src_indices Indices of intersecting source elements.
 * @param is_count_only If true, only counts intersections; if false, also fills
 * tgt2src_indices.
 *
 * @note This method assumes 2D linear triangles and uses
 * `r3d::intersect_simplices` for geometric intersection.
 *
 * @see r3d::intersect_simplices, intersectTargets
 */

void FindIntersections::adjBasedIntersectSearch(
  const Omega_h::LOs& tgt2src_offsets,
  Omega_h::Write<Omega_h::LO>& nIntersections,
  Omega_h::Write<Omega_h::LO>& tgt2src_indices, bool is_count_only)
{

  const auto& tgt_coords = target_mesh_.coords();
  const auto& src_coords = source_mesh_.coords();
  const auto& tgt_faces2nodes =
    target_mesh_.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& src_faces2nodes =
    source_mesh_.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& src_elem_areas = measure_elements_real(&source_mesh_);
  const auto& tgt_elem_areas = measure_elements_real(&target_mesh_);
  const auto& t2t =
    source_mesh_.ask_dual(); // gives connected element neighbors
  const auto& t2tt = t2t.a2ab;
  const auto& tt2t = t2t.ab2b;

  auto centroids = compute_centroid(target_mesh_);

  pcms::GridPointSearch search_cell(source_mesh_, 20, 20);
  auto results = search_cell(centroids);

  auto nfaces_target = target_mesh_.nfaces();
  Omega_h::parallel_for(
    nfaces_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;

      auto current_cell_id = results(id).tri_id;
      auto current_tgt_elm_area = tgt_elem_areas[id];

      OMEGA_H_CHECK_PRINTF(current_cell_id >= 0,
                           "ERROR: source cell id not found for given target "
                           "centroid %d (%f, %f)\n",
                           id, centroids(id, 0), centroids(id, 1));

      auto tgt_elm_vert_coords =
        get_vert_coords_of_elem(tgt_coords, tgt_faces2nodes, id);

      Omega_h::LO start_counter;
      if (!is_count_only) {
        start_counter = tgt2src_offsets[id];
      }

      int count = 0;

      count++;
      visited.push_back(current_cell_id);
      queue.push_back(current_cell_id);

      if (!is_count_only) {
        int idx_count = count - 1;
        tgt2src_indices[start_counter + idx_count] = current_cell_id;
      }

      while (!queue.isEmpty()) {
        Omega_h::LO currentElm = queue.front();
        queue.pop_front();
        auto start = t2tt[currentElm];
        auto end = t2tt[currentElm + 1];

        for (int i = start; i < end; ++i) {
          auto neighborElmId = tt2t[i];

          if (visited.notVisited(neighborElmId)) {
            visited.push_back(neighborElmId);
            auto elm_vert_coords = get_vert_coords_of_elem(
              src_coords, src_faces2nodes, neighborElmId);
            r3d::Polytope<2> intersection;
            r3d::intersect_simplices(intersection, tgt_elm_vert_coords,
                                     elm_vert_coords);
            auto intersected_area = r3d::measure(intersection);
            auto current_src_elm_area = src_elem_areas[neighborElmId];
            auto scale =
              Kokkos::fmax(current_tgt_elm_area, current_src_elm_area);
            auto eps = Kokkos::fmax(abs_tol, rel_tol * scale);
            if (intersection.nverts >= 3 && intersected_area >= eps) {
              count++;

              OMEGA_H_CHECK_PRINTF(
                count < 500, "WARNING: count exceeds 500 for target %d", id);

              queue.push_back(neighborElmId);

              if (!is_count_only) {
                Omega_h::LO idx_count = count - 1;
                tgt2src_indices[start_counter + idx_count] = neighborElmId;

              } // end of tgt2src_indices check

            } // end of intersection with bbox check

          } // end of not visited check

        } // end of loop over adj elements to the current element

      } // end of while loop

      nIntersections[id] = count;
    }, // end of lambda
    "count the number of intersections for each target element");
}

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
                                     Omega_h::Mesh& target_mesh)
{
  FindIntersections intersect(source_mesh, target_mesh);

  auto nfaces_target = target_mesh.nfaces();

  Omega_h::Write<Omega_h::LO> nIntersections(
    nfaces_target, 0, "number of intersections in each target vertex");

  Omega_h::Write<Omega_h::LO> tgt2src_indices;

  intersect.adjBasedIntersectSearch(Omega_h::LOs(), nIntersections,
                                    tgt2src_indices, true);

  Kokkos::fence();
  auto tgt2src_offsets = Omega_h::offset_scan(Omega_h::Read(nIntersections),
                                              "offsets for intersections");
  auto ntotal_intersections = tgt2src_offsets.last();

  Kokkos::fence();

  tgt2src_indices = Omega_h::Write<Omega_h::LO>(
    ntotal_intersections, 0,
    "indices of the source elements that intersect the given target element");

  intersect.adjBasedIntersectSearch(tgt2src_offsets, nIntersections,
                                    tgt2src_indices, false);
  return {.tgt2src_offsets = tgt2src_offsets,
          .tgt2src_indices = Omega_h::read(tgt2src_indices)};
}

#endif
