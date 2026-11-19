#include "pcms/backends/omega_h/mesh_intersection.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/mesh_geometry.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"

namespace pcms
{
using utility::ConvertCoordsTo2D;
using utility::get_entity_centroids;
} // namespace pcms

namespace pcms::omega_h
{
namespace
{
// For each target element, the source element containing or nearest its
// centroid, where the BFS over the source dual graph starts. An outside
// centroid does not imply that the target element has no overlap.
template <int Dim>
Kokkos::View<const LO*, DeviceMemorySpace> LocateTargetCentroids(
  Omega_h::Mesh& target_mesh, const OmegaHGeometry& source_geometry)
{
  const auto flat_centroids = pcms::get_entity_centroids(target_mesh, Dim);
  // Convert layout_right 1D Omega_h array to 2D Kokkos view with correct layout
  auto centroids = ConvertCoordsTo2D(flat_centroids, target_mesh.nelems(), Dim);
  return source_geometry.FindElementSeeds(centroids);
}
} // namespace

template <int Dim>
void FindIntersections::adjBasedIntersectSearch(
  const Kokkos::View<const LO*>& start_elements,
  const Omega_h::LOs& tgt2src_offsets,
  Omega_h::Write<Omega_h::LO>& nIntersections,
  Omega_h::Write<Omega_h::LO>& tgt2src_indices, bool is_count_only,
  bool use_prefilter)
{
  // Element entity dimension equals the spatial dimension (FACE for 2D, REGION
  // for 3D); measures are triangle areas (2D) or tet volumes (3D).
  const auto& tgt_coords = target_mesh_.coords();
  const auto& src_coords = source_mesh_.coords();
  const auto& tgt_elems2nodes = target_mesh_.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto& src_elems2nodes = source_mesh_.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto& src_elem_measures = measure_elements_real(&source_mesh_);
  const auto& tgt_elem_measures = measure_elements_real(&target_mesh_);
  const auto& t2t =
    source_mesh_.ask_dual(); // gives connected element neighbors
  const auto& t2tt = t2t.a2ab;
  const auto& tt2t = t2t.ab2b;

  auto nelems_target = target_mesh_.nelems();
  PCMS_ALWAYS_ASSERT(static_cast<Omega_h::LO>(start_elements.extent(0)) ==
                     nelems_target);
  Omega_h::parallel_for(
    nelems_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;

      auto current_cell_id = start_elements(id);
      auto current_tgt_elm_measure = tgt_elem_measures[id];

      // The point search supplied no seed for this centroid's grid cell.
      if (current_cell_id < 0) {
        nIntersections[id] = 0;
        return;
      }

      auto tgt_elm_vert_coords =
        get_vert_coords_of_elem<Dim>(tgt_coords, tgt_elems2nodes, id);

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
            // If the visited buffer is full, skip this neighbor so the BFS
            // terminates. Without this, an unrecorded neighbor stays "not
            // visited" and is re-queued forever (infinite loop). Mirrors the
            // guard in adj_search.cpp. Raise PCMS_INTERSECTION_TRACK_SIZE if
            // this fires.
            if (!visited.push_back(neighborElmId)) {
              printf("ERROR: visited buffer full "
                     "(PCMS_INTERSECTION_TRACK_SIZE=%d) for "
                     "target %d; some intersections may be missed\n",
                     PCMS_INTERSECTION_TRACK_SIZE, id);
              continue;
            }
            auto elm_vert_coords = get_vert_coords_of_elem<Dim>(
              src_coords, src_elems2nodes, neighborElmId);

            // use plane to check if element is already outside clipped region
            // this avoids construction of unnecessary clip
            if (use_prefilter && simplices_have_degenerate_overlap<Dim>(
                                   tgt_elm_vert_coords, elm_vert_coords)) {
              continue;
            }

            r3d::Polytope<Dim> intersection;
            r3d::intersect_simplices(intersection, tgt_elm_vert_coords,
                                     elm_vert_coords);
            // Take the magnitude: r3d::measure is signed by the orientation of
            // the target simplex used to initialize the polytope, which can be
            // negative for tetrahedra. This mirrors the fabs applied in the
            // sub-simplex decomposition and mass assembly; without it a
            // negatively-oriented target element would reject all of its real
            // overlaps and break conservation.
            auto intersected_measure = Kokkos::fabs(r3d::measure(intersection));
            auto current_src_elm_measure = src_elem_measures[neighborElmId];
            auto scale =
              Kokkos::fmax(current_tgt_elm_measure, current_src_elm_measure);
            auto eps = Kokkos::fmax(PCMS_INTERSECTION_ABS_TOL,
                                    PCMS_INTERSECTION_REL_TOL * scale);
            // A valid intersection is a non-degenerate simplex-simplex overlap:
            // at least Dim+1 vertices (a polygon in 2D, a polyhedron in 3D).
            if (intersection.nverts >= Dim + 1 && intersected_measure >= eps) {
              count++;

              OMEGA_H_CHECK_PRINTF(
                count < PCMS_INTERSECTION_QUEUE_SIZE,
                "intersection count for target %d reached the cap %d; raise "
                "PCMS_INTERSECTION_QUEUE_SIZE",
                id, PCMS_INTERSECTION_QUEUE_SIZE);

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

// Explicit instantiations for the supported spatial dimensions.
template void FindIntersections::adjBasedIntersectSearch<2>(
  const Kokkos::View<const LO*>&, const Omega_h::LOs&,
  Omega_h::Write<Omega_h::LO>&, Omega_h::Write<Omega_h::LO>&, bool, bool);
template void FindIntersections::adjBasedIntersectSearch<3>(
  const Kokkos::View<const LO*>&, const Omega_h::LOs&,
  Omega_h::Write<Omega_h::LO>&, Omega_h::Write<Omega_h::LO>&, bool, bool);

namespace
{
template <int Dim>
IntersectionResults intersectTargetsImpl(const OmegaHGeometry& source,
                                         Omega_h::Mesh& target_mesh,
                                         bool use_prefilter)
{
  FindIntersections intersect(source.GetMesh(), target_mesh);

  auto nelems_target = target_mesh.nelems();

  const auto start_elements = LocateTargetCentroids<Dim>(target_mesh, source);

  Omega_h::Write<Omega_h::LO> nIntersections(
    nelems_target, 0, "number of intersections in each target element");

  Omega_h::Write<Omega_h::LO> tgt2src_indices;

  intersect.adjBasedIntersectSearch<Dim>(start_elements, Omega_h::LOs(),
                                         nIntersections, tgt2src_indices, true,
                                         use_prefilter);

  Kokkos::fence();
  auto tgt2src_offsets = Omega_h::offset_scan(Omega_h::Read(nIntersections),
                                              "offsets for intersections");
  auto ntotal_intersections = tgt2src_offsets.last();

  Kokkos::fence();

  tgt2src_indices = Omega_h::Write<Omega_h::LO>(
    ntotal_intersections, 0,
    "indices of the source elements that intersect the given target element");

  intersect.adjBasedIntersectSearch<Dim>(start_elements, tgt2src_offsets,
                                         nIntersections, tgt2src_indices, false,
                                         use_prefilter);
  return {.tgt2src_offsets = tgt2src_offsets,
          .tgt2src_indices = Omega_h::read(tgt2src_indices)};
}

} // namespace

IntersectionResults intersectTargets(const OmegaHGeometry& source,
                                     Omega_h::Mesh& target_mesh,
                                     bool use_prefilter)
{
  if (source.GetMesh().dim() != target_mesh.dim()) {
    throw pcms_error("intersectTargets: the source geometry's mesh and the "
                     "target mesh differ in dimension");
  }
  if (target_mesh.dim() == 3) {
    return intersectTargetsImpl<3>(source, target_mesh, use_prefilter);
  }
  return intersectTargetsImpl<2>(source, target_mesh, use_prefilter);
}
} // namespace pcms::omega_h
