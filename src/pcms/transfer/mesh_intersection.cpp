#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/utility/mesh_geometry.h"

namespace pcms
{
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

  const auto flat_centroids =
    pcms::get_entity_centroids(target_mesh_, Omega_h::FACE);
  auto centroids = Kokkos::View<const Omega_h::Real* [2]>(
    flat_centroids.data(), target_mesh_.nfaces());

  pcms::GridPointSearch2D search_cell(source_mesh_, 20, 20);
  auto results = search_cell(centroids);

  auto nfaces_target = target_mesh_.nfaces();
  Omega_h::parallel_for(
    nfaces_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;

      auto current_cell_id = results(id).element_id;
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
} // namespace pcms
