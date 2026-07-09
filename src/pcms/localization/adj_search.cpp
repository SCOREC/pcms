#include "pcms/localization/adj_search.hpp"
#include <Omega_h_array_ops.hpp>

namespace pcms
{

// Debug helper retained intentionally: useful for diagnosing point-localization
// failures while developing support-search logic.
[[maybe_unused]] static void checkTargetPoints(
  const Kokkos::View<pcms::GridPointSearch2D::Result*>& results)
{
  Kokkos::fence();
  pcms::printInfo("INFO: Checking target points...\n");
  auto check_target_points = OMEGA_H_LAMBDA(Omega_h::LO i)
  {
    if (results(i).element_id < 0) {
      OMEGA_H_CHECK_PRINTF(results(i).element_id >= 0,
                           "ERROR: Source cell id not found for target %d\n",
                           i);
      printf("%d, ", i);
    }
  };
  Omega_h::parallel_for(results.size(), check_target_points,
                        "check_target_points");
  Kokkos::fence();
  pcms::printInfo("\n");
}

// Debug helper retained intentionally: useful for tracing CSR support contents
// for a specific target id during local debugging.
[[maybe_unused]] static void printSupportsForTarget(
  const Omega_h::LO target_id, const Omega_h::Write<Omega_h::LO>& supports_ptr,
  const Omega_h::Write<Omega_h::LO>& nSupports,
  const Omega_h::Write<Omega_h::LO>& support_idx)
{
  Omega_h::parallel_for(
    nSupports.size(), OMEGA_H_LAMBDA(const Omega_h::LO id) {
      if (id == target_id) {
        Omega_h::LO start_ptr = supports_ptr[id];
        Omega_h::LO end_ptr = supports_ptr[id + 1];
        printf("Target vertex: %d\n with %d num supports: nSupports[id]=%d", id,
               end_ptr - start_ptr, nSupports[id]);
        for (Omega_h::LO i = start_ptr; i < end_ptr; ++i) {
          Omega_h::LO cell_id = support_idx[i];
          printf(", %d", cell_id);
        }
        printf("\n");
      }
    });
}

static Omega_h::Write<Omega_h::LO> build_support_offsets(
  const Omega_h::LO nvertices_target,
  const Omega_h::Write<Omega_h::LO>& nSupports, Omega_h::LO& total_supports)
{
  auto supports_ptr = Omega_h::Write<Omega_h::LO>(
    nvertices_target + 1, 0, "number of support source vertices in CSR format");

  total_supports = 0;
  Kokkos::parallel_scan(
    nvertices_target,
    OMEGA_H_LAMBDA(int j, int& update, bool final) {
      update += nSupports[j];
      if (final) {
        supports_ptr[j + 1] = update;
      }
    },
    total_supports);

  return supports_ptr;
}

static Omega_h::Write<Omega_h::LO> locate_target_cells(
  Omega_h::Mesh& source_mesh, const Omega_h::Reals& target_coords,
  Omega_h::LO nvertices_target)
{
  const auto dim = source_mesh.dim();
  auto source_cell_ids =
    Omega_h::Write<Omega_h::LO>(nvertices_target, -1, "source cell ids");

  if (dim == 2) {
    Kokkos::View<Omega_h::Real* [2]> target_points("target_points",
                                                   nvertices_target);
    Omega_h::parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
        target_points(i, 0) = target_coords[i * dim];
        target_points(i, 1) = target_coords[i * dim + 1];
      });
    Kokkos::fence();

    pcms::GridPointSearch2D search_cell(source_mesh, 10, 10);
    auto results = search_cell(target_points);
    Omega_h::parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
        auto source_cell_id = results(i).element_id;
        if (source_cell_id < 0)
          source_cell_id = Kokkos::abs(source_cell_id);
        OMEGA_H_CHECK_PRINTF(
          source_cell_id >= 0,
          "ERROR: Source cell id not found for target %d (%f,%f)\n", i,
          target_points(i, 0), target_points(i, 1));
        source_cell_ids[i] = source_cell_id;
      });
  } else if (dim == 3) {
    Kokkos::View<Omega_h::Real* [3]> target_points("target_points",
                                                   nvertices_target);
    Omega_h::parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
        target_points(i, 0) = target_coords[i * dim];
        target_points(i, 1) = target_coords[i * dim + 1];
        target_points(i, 2) = target_coords[i * dim + 2];
      });
    Kokkos::fence();

    pcms::GridPointSearch3D search_cell(source_mesh, 10, 10, 10);
    auto results = search_cell(target_points);
    Omega_h::parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
        auto source_cell_id = results(i).element_id;
        if (source_cell_id < 0)
          source_cell_id = Kokkos::abs(source_cell_id);
        OMEGA_H_CHECK_PRINTF(source_cell_id >= 0,
                             "ERROR: Source cell id not found for target %d\n",
                             i);
        source_cell_ids[i] = source_cell_id;
      });
  } else {
    throw pcms_error("Unsupported dimension in locate_target_cells");
  }

  return source_cell_ids;
}

OMEGA_H_INLINE void copy_from_rank1_view(const Omega_h::Reals& view,
                                         const Omega_h::LO id,
                                         const Omega_h::LO dim,
                                         Omega_h::Real* out)
{
  for (Omega_h::LO k = 0; k < dim; ++k) {
    out[k] = view[id * dim + k];
  }
}

OMEGA_H_INLINE void copy_from_rank2_view(
  const Kokkos::View<Omega_h::Real* [2]>& view, const Omega_h::LO id,
  const Omega_h::LO dim, Omega_h::Real* out)
{
  for (Omega_h::LO k = 0; k < dim; ++k) {
    out[k] = view(id, k);
  }
}

OMEGA_H_INLINE void add_support_if_unvisited_and_within_cutoff(
  const Omega_h::LO candidate_id, const Omega_h::Real cutoff_distance,
  const Omega_h::Real* target_coords, const Omega_h::Reals& support_coords_view,
  const Omega_h::LO dim, Track& visited, Queue& queue, int& count,
  const bool is_build_csr_call, const Omega_h::LO start_counter,
  Omega_h::Write<Omega_h::LO> support_idx)
{
  if (!visited.notVisited(candidate_id))
    return;

  // If the visited buffer is full, stop here so the BFS terminates. Without
  // this, unrecorded vertices are treated as unvisited and re-queued forever.
  if (!visited.push_back(candidate_id))
    return;

  Omega_h::Real candidate_coords[max_dim];
  for (Omega_h::LO k = 0; k < dim; ++k) {
    candidate_coords[k] = support_coords_view[candidate_id * dim + k];
  }

  const Omega_h::Real dist =
    pcms::distance_squared(target_coords, candidate_coords, dim);
  if (dist <= cutoff_distance) {
    count++;
    queue.push_back(candidate_id);
    if (!is_build_csr_call) {
      const Omega_h::LO idx_count = count - 1;
      support_idx[start_counter + idx_count] = candidate_id;
    }
  }
}

static void adjBasedSearchFromPoints(
  Omega_h::Mesh& source_mesh, const Omega_h::Reals& target_coords,
  const Omega_h::Write<Omega_h::LO>& source_cell_ids,
  Omega_h::Write<Omega_h::LO>& supports_ptr,
  Omega_h::Write<Omega_h::LO>& nSupports,
  Omega_h::Write<Omega_h::LO>& support_idx,
  Omega_h::Write<Omega_h::Real>& radii2, bool is_build_csr_call)
{
  const auto& sourcePoints_coords = source_mesh.coords();
  const auto dim = source_mesh.dim();
  const auto nvertices_target = nSupports.size();
  OMEGA_H_CHECK(radii2.size() == nvertices_target);

  const auto& vert2vert = source_mesh.ask_star(Omega_h::VERT);
  const auto& v2v_ptr = vert2vert.a2ab;
  const auto& v2v_data = vert2vert.ab2b;
  const auto& cells2verts = source_mesh.ask_verts_of(dim);

  Omega_h::parallel_for(
    nvertices_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;
      Omega_h::Real cutoffDistance = radii2[id];
      Omega_h::LO source_cell_id = source_cell_ids[id];

      OMEGA_H_CHECK_PRINTF(source_cell_id >= 0,
                           "ERROR: Source cell id not found for target %d\n",
                           id);

      const Omega_h::LO num_verts_in_dim = dim + 1;
      Omega_h::LO start_ptr = source_cell_id * num_verts_in_dim;
      Omega_h::LO end_ptr = start_ptr + num_verts_in_dim;
      Omega_h::Real target_coords_i[max_dim];
      copy_from_rank1_view(target_coords, id, dim, target_coords_i);

      Omega_h::LO start_counter = 0;
      if (!is_build_csr_call) {
        start_counter = supports_ptr[id];
      }

      int count = 0;
      for (Omega_h::LO i = start_ptr; i < end_ptr; ++i) {
        Omega_h::LO vert_id = cells2verts[i];
        const int count_before = count;
        add_support_if_unvisited_and_within_cutoff(
          vert_id, cutoffDistance, target_coords_i, sourcePoints_coords, dim,
          visited, queue, count, is_build_csr_call, start_counter, support_idx);
        if (count > count_before && count >= 500) {
          printf("Warning: count exceeds 500 for target %d with %d supports\n",
                 id, end_ptr - start_ptr);
        }
      }

      while (!queue.isEmpty()) {
        Omega_h::LO currentVertex = queue.front();
        queue.pop_front();
        Omega_h::LO start = v2v_ptr[currentVertex];
        Omega_h::LO end = v2v_ptr[currentVertex + 1];

        for (Omega_h::LO i = start; i < end; ++i) {
          auto neighborIndex = v2v_data[i];

          const int count_before = count;
          add_support_if_unvisited_and_within_cutoff(
            neighborIndex, cutoffDistance, target_coords_i, sourcePoints_coords,
            dim, visited, queue, count, is_build_csr_call, start_counter,
            support_idx);
          if (count > count_before && count >= 500) {
            printf("Warning: count exceeds 500 for target %d with start %d "
                   "and end %d radius2 %f adding neighbor %d\n",
                   id, start, end, cutoffDistance, neighborIndex);
          }
        }
      }

      nSupports[id] = count;
    },
    "count the number of supports in each target point");
}

void FindSupports::adjBasedSearch(Omega_h::Write<Omega_h::LO>& supports_ptr,
                                  Omega_h::Write<Omega_h::LO>& nSupports,
                                  Omega_h::Write<Omega_h::LO>& support_idx,
                                  Omega_h::Write<Omega_h::Real>& radii2,
                                  bool is_build_csr_call)
{
  const auto& targetPoints_coords = target_mesh.coords();
  const auto nvertices_target = target_mesh.nverts();
  auto source_cell_ids =
    locate_target_cells(source_mesh, targetPoints_coords, nvertices_target);
  adjBasedSearchFromPoints(source_mesh, targetPoints_coords, source_cell_ids,
                           supports_ptr, nSupports, support_idx, radii2,
                           is_build_csr_call);
}

void FindSupports::adjBasedSearchCentroidNodes(
  Omega_h::Write<Omega_h::LO>& supports_ptr,
  Omega_h::Write<Omega_h::LO>& nSupports,
  Omega_h::Write<Omega_h::LO>& support_idx,
  Omega_h::Write<Omega_h::Real>& radii2, bool is_build_csr_call)
{
  const auto& mesh_coords = source_mesh.coords();
  const auto& nvertices = source_mesh.nverts();
  const auto& dim = source_mesh.dim();

  const auto& nodes2faces = source_mesh.ask_up(Omega_h::VERT, Omega_h::FACE);
  const auto& n2f_ptr = nodes2faces.a2ab;
  const auto& n2f_data = nodes2faces.ab2b;
  const auto& faces2nodes =
    source_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  OMEGA_H_CHECK_PRINTF(source_mesh.dim() == 2,
                       "Only 2D meshes are supported but found %d\n",
                       source_mesh.dim());
  auto cell_centroids = pcms::get_entity_centroids(source_mesh, Omega_h::FACE);

  Omega_h::parallel_for(
    nvertices,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;
      const Omega_h::LO num_verts_in_dim = dim + 1;
      Omega_h::Real target_coords[max_dim];
      Omega_h::Real cutoffDistance = radii2[id];

      copy_from_rank1_view(mesh_coords, id, dim, target_coords);

      Omega_h::LO start_counter = 0;
      if (!is_build_csr_call) {
        start_counter = supports_ptr[id];
      }
      Omega_h::LO start_ptr = n2f_ptr[id];
      Omega_h::LO end_ptr = n2f_ptr[id + 1];

      int count = 0;
      for (Omega_h::LO i = start_ptr; i < end_ptr; ++i) {
        Omega_h::LO cell_id = n2f_data[i];
        add_support_if_unvisited_and_within_cutoff(
          cell_id, cutoffDistance, target_coords, cell_centroids, dim, visited,
          queue, count, is_build_csr_call, start_counter, support_idx);
      }

      while (!queue.isEmpty()) {
        Omega_h::LO currentCell = queue.front();
        queue.pop_front();
        Omega_h::LO start = currentCell * num_verts_in_dim;
        Omega_h::LO end = start + num_verts_in_dim;

        for (Omega_h::LO i = start; i < end; ++i) {
          Omega_h::LO current_vert_id = faces2nodes[i];
          Omega_h::LO start_ptr_current_vert = n2f_ptr[current_vert_id];
          Omega_h::LO end_ptr_vert_current_vert = n2f_ptr[current_vert_id + 1];
          for (Omega_h::LO j = start_ptr_current_vert;
               j < end_ptr_vert_current_vert; ++j) {
            auto neighbor_cell_index = n2f_data[j];

            // check if neighbor index is already in the queue to be checked
            // TODO refactor this into a function

            add_support_if_unvisited_and_within_cutoff(
              neighbor_cell_index, cutoffDistance, target_coords,
              cell_centroids, dim, visited, queue, count, is_build_csr_call,
              start_counter, support_idx);
          }
        }
      }

      nSupports[id] = count;
    },
    "count the number of supports in each target point");
}

SupportResults searchNeighbors(Omega_h::Mesh& source_mesh,
                               const Omega_h::Reals& target_coords,
                               Omega_h::LO nvertices_target,
                               Omega_h::Real& cutoffDistance,
                               Omega_h::LO min_req_support,
                               Omega_h::LO max_allowed_support,
                               bool adapt_radius)
{
  Omega_h::Write<Omega_h::LO> nSupports(
    nvertices_target, 0, "number of supports in each target vertex");
  pcms::printInfo("INFO: Cut off distance: %f\n", cutoffDistance);
  Omega_h::Write<Omega_h::Real> radii2 = Omega_h::Write<Omega_h::Real>(
    nvertices_target, cutoffDistance, "squared radii of the supports");

  Omega_h::Write<Omega_h::LO> supports_ptr;
  Omega_h::Write<Omega_h::LO> supports_idx;
  auto source_cell_ids =
    locate_target_cells(source_mesh, target_coords, nvertices_target);

  if (!adapt_radius) {
    pcms::printInfo("INFO: Fixed radius search *(disregarding required minimum "
                    "support)*... \n");
    adjBasedSearchFromPoints(source_mesh, target_coords, source_cell_ids,
                             supports_ptr, nSupports, supports_idx, radii2,
                             true);
  } else {
    pcms::printInfo("INFO: Adaptive radius search... \n");
    int r_adjust_loop = 0;
    while (true) {
      nSupports = Omega_h::Write<Omega_h::LO>(
        nvertices_target, 0, "number of supports in each target vertex");

      const auto max_radius = Omega_h::get_max(Omega_h::read(radii2));
      pcms::printInfo("INFO: Loop %d: max_radius: %f\n", r_adjust_loop,
                      max_radius);

      adjBasedSearchFromPoints(source_mesh, target_coords, source_cell_ids,
                               supports_ptr, nSupports, supports_idx, radii2,
                               true);
      Kokkos::fence();

      unsigned min_supports_found = 0;
      unsigned max_supports_found = 0;
      minmax(nSupports, min_supports_found, max_supports_found);
      r_adjust_loop++;
      pcms::printInfo(
        "Iter: %d min_nSupports: %d max_nSupports: %d, max_radius %f\n",
        r_adjust_loop, min_supports_found, max_supports_found, max_radius);

      if (within_number_of_support_range(min_supports_found, max_supports_found,
                                         min_req_support,
                                         3 * min_req_support)) {
        break;
      }

      adapt_radii(min_req_support, max_allowed_support, nvertices_target,
                  radii2, nSupports);
    }

    pcms::printInfo("INFO: Took %d loops to adjust the radius\n",
                    r_adjust_loop);
  }

  Omega_h::LO total_supports = 0;
  supports_ptr =
    build_support_offsets(nvertices_target, nSupports, total_supports);

  Kokkos::fence();

  supports_idx = Omega_h::Write<Omega_h::LO>(
    total_supports, 0, "index of source supports of each target node");

  adjBasedSearchFromPoints(source_mesh, target_coords, source_cell_ids,
                           supports_ptr, nSupports, supports_idx, radii2,
                           false);

  return SupportResults{read(supports_ptr), read(supports_idx), radii2};
}

SupportResults searchNeighbors(Omega_h::Mesh& source_mesh,
                               Omega_h::Mesh& target_mesh,
                               Omega_h::Real& cutoffDistance,
                               Omega_h::LO min_req_support,
                               Omega_h::LO max_allowed_support,
                               bool adapt_radius)
{
  auto target_coords = target_mesh.coords();
  auto result = searchNeighbors(
    source_mesh, target_coords, target_mesh.nverts(), cutoffDistance,
    min_req_support, max_allowed_support, adapt_radius);
  target_mesh.add_tag<Omega_h::Real>(Omega_h::VERT, "radii2", 1, result.radii2);
  return result;
}

SupportResults searchNeighbors(Omega_h::Mesh& mesh,
                               Omega_h::Real cutoffDistance,
                               Omega_h::LO min_support, bool adapt_radius)
{
  Omega_h::Write<Omega_h::LO> supports_ptr;
  Omega_h::Write<Omega_h::LO> supports_idx;
  FindSupports search(mesh);
  Omega_h::LO nvertices_target = mesh.nverts();
  Omega_h::Write<Omega_h::LO> nSupports(
    nvertices_target, 0, "number of supports in each target vertex");

  pcms::printInfo("INFO: Inside searchNeighbors 1\n");
  Omega_h::Write<Omega_h::Real> radii2 = Omega_h::Write<Omega_h::Real>(
    nvertices_target, cutoffDistance, "squared radii of the supports");
  pcms::printInfo("INFO: Cutoff distance: %f\n", cutoffDistance);

  if (!adapt_radius) {
    pcms::printInfo("INFO: Fixed radius search *(disregarding required minimum "
                    "support)* ... \n");
    search.adjBasedSearch(supports_ptr, nSupports, supports_idx, radii2, true);
  } else {
    pcms::printInfo("INFO: Adaptive radius search... \n");
    int r_adjust_loop = 0;
    while (true) { // until the number of minimum support is met
      const auto max_radius = Omega_h::get_max(Omega_h::read(radii2));
      pcms::printInfo("INFO: Loop %d: max_radius: %f\n", r_adjust_loop,
                      max_radius);

      nSupports = Omega_h::Write<Omega_h::LO>(
        nvertices_target, 0, "number of supports in each target vertex");
      search.adjBasedSearchCentroidNodes(supports_ptr, nSupports, supports_idx,
                                         radii2, true);

      Kokkos::fence();
      unsigned min_nSupports = 0;
      unsigned max_nSupports = 0;
      minmax(nSupports, min_nSupports, max_nSupports);

      r_adjust_loop++;
      pcms::printInfo(
        "Iter: %d min_nSupports: %d max_nSupports: %d at loop %d, "
        "max_radius %f\n",
        r_adjust_loop, min_nSupports, max_nSupports, r_adjust_loop, max_radius);

      if (within_number_of_support_range(min_nSupports, max_nSupports,
                                         min_support, 3 * min_support)) {
        break;
      }

      adapt_radii(min_support, 3 * min_support, radii2.size(), radii2,
                  nSupports);
    }
    pcms::printInfo("INFO: Took %d loops to adjust the radius\n",
                    r_adjust_loop);
  }

  // offset array for the supports of each target vertex
  Omega_h::LO total_supports = 0;
  supports_ptr =
    build_support_offsets(nvertices_target, nSupports, total_supports);

  pcms::printInfo("INFO: Inside searchNeighbors 3\n");
  Kokkos::fence();

  supports_idx = Omega_h::Write<Omega_h::LO>(
    total_supports, 0, "index of source supports of each target node");
  pcms::printInfo("INFO: Total_supports: %d\n", total_supports);

  search.adjBasedSearchCentroidNodes(supports_ptr, nSupports, supports_idx,
                                     radii2, false);

  mesh.add_tag<Omega_h::Real>(Omega_h::VERT, "support_radius", 1, radii2);
  return SupportResults{read(supports_ptr), read(supports_idx), radii2};
}

} // namespace pcms
