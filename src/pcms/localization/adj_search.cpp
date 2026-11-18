#include "pcms/localization/adj_search.hpp"
#include "pcms/localization/queue_visited.hpp"
#include <Omega_h_array_ops.hpp>

namespace pcms
{

static constexpr int max_dim = 3;

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

OMEGA_H_INLINE void copy_from_rank1_view(const Omega_h::Reals& view,
                                         const Omega_h::LO id,
                                         const Omega_h::LO dim,
                                         Omega_h::Real* out)
{
  for (Omega_h::LO k = 0; k < dim; ++k) {
    out[k] = view[id * dim + k];
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

// One pass of the breadth-first search of every query point: counts the
// supports within each point's squared radius or, if !is_build_csr_call, also
// writes them at the point's CSR offset.
static void adjBasedSearch(Omega_h::Mesh& mesh, int site_dim,
                           const Omega_h::Reals& target_coords, int seed_dim,
                           Kokkos::View<const LO*, DeviceMemorySpace> seeds,
                           Omega_h::Write<Omega_h::LO>& supports_ptr,
                           Omega_h::Write<Omega_h::LO>& nSupports,
                           Omega_h::Write<Omega_h::LO>& support_idx,
                           Omega_h::Write<Omega_h::Real>& radii2,
                           bool is_build_csr_call)
{
  const int dim = mesh.dim();
  const int verts_per_elem = dim + 1;
  const bool vertex_sites = site_dim == 0;
  const bool element_seeds = seed_dim == dim;
  const auto site_coords =
    vertex_sites ? mesh.coords() : pcms::get_entity_centroids(mesh, dim);
  const auto elems2verts = mesh.ask_elem_verts();
  const auto verts2elems = mesh.ask_up(Omega_h::VERT, dim);
  const auto v2e_ptr = verts2elems.a2ab;
  const auto v2e = verts2elems.ab2b;
  const auto star =
    vertex_sites ? mesh.ask_star(Omega_h::VERT) : Omega_h::Adj{};
  const auto v2v_ptr = star.a2ab;
  const auto v2v = star.ab2b;
  OMEGA_H_CHECK(radii2.size() == nSupports.size());

  Omega_h::parallel_for(
    nSupports.size(),
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      const Omega_h::LO seed = seeds(id);
      if (seed < 0) {
        nSupports[id] = 0;
        return;
      }
      Queue queue;
      Track visited;
      const Omega_h::Real cutoff = radii2[id];
      Omega_h::Real target[max_dim];
      copy_from_rank1_view(target_coords, id, dim, target);
      const Omega_h::LO start_counter =
        is_build_csr_call ? 0 : supports_ptr[id];
      int count = 0;
      const auto add = [&](const Omega_h::LO site) {
        add_support_if_unvisited_and_within_cutoff(
          site, cutoff, target, site_coords, dim, visited, queue, count,
          is_build_csr_call, start_counter, support_idx);
      };

      // the sites of the seed's elements
      const Omega_h::LO first = element_seeds ? 0 : v2e_ptr[seed];
      const Omega_h::LO last = element_seeds ? 1 : v2e_ptr[seed + 1];
      for (Omega_h::LO k = first; k < last; ++k) {
        const Omega_h::LO elem = element_seeds ? seed : v2e[k];
        if (vertex_sites) {
          for (int v = 0; v < verts_per_elem; ++v) {
            add(elems2verts[elem * verts_per_elem + v]);
          }
        } else {
          add(elem);
        }
      }
      // then their neighbours: vertices sharing an edge, or elements sharing
      // a vertex
      while (!queue.isEmpty()) {
        const Omega_h::LO site = queue.front();
        queue.pop_front();
        if (vertex_sites) {
          for (Omega_h::LO k = v2v_ptr[site]; k < v2v_ptr[site + 1]; ++k) {
            add(v2v[k]);
          }
        } else {
          for (int v = 0; v < verts_per_elem; ++v) {
            const Omega_h::LO vert = elems2verts[site * verts_per_elem + v];
            for (Omega_h::LO k = v2e_ptr[vert]; k < v2e_ptr[vert + 1]; ++k) {
              add(v2e[k]);
            }
          }
        }
      }
      nSupports[id] = count;
    },
    "count the number of supports in each target point");
}

SupportResults searchNeighbors(
  Omega_h::Mesh& mesh, int site_dim, const Omega_h::Reals& target_coords,
  int seed_dim, Kokkos::View<const LO*, DeviceMemorySpace> seeds,
  Omega_h::Real cutoffDistance, Omega_h::LO min_req_support,
  Omega_h::LO max_allowed_support, bool adapt_radius)
{
  if ((site_dim != 0 && site_dim != mesh.dim()) ||
      (seed_dim != 0 && seed_dim != mesh.dim())) {
    throw pcms_error("searchNeighbors: the sites and seeds must be the "
                     "vertices or the elements of the mesh");
  }
  const auto nvertices_target = static_cast<Omega_h::LO>(seeds.extent(0));
  Omega_h::Write<Omega_h::LO> nSupports(
    nvertices_target, 0, "number of supports in each target vertex");
  pcms::printInfo("INFO: Cut off distance: %f\n", cutoffDistance);
  Omega_h::Write<Omega_h::Real> radii2 = Omega_h::Write<Omega_h::Real>(
    nvertices_target, cutoffDistance, "squared radii of the supports");

  Omega_h::Write<Omega_h::LO> supports_ptr;
  Omega_h::Write<Omega_h::LO> supports_idx;

  if (!adapt_radius) {
    pcms::printInfo("INFO: Fixed radius search *(disregarding required minimum "
                    "support)*... \n");
    adjBasedSearch(mesh, site_dim, target_coords, seed_dim, seeds, supports_ptr,
                   nSupports, supports_idx, radii2, true);
  } else {
    pcms::printInfo("INFO: Adaptive radius search... \n");
    // Cap the adaptation so targets that can never reach the required support
    // count do not spin the loop forever.
    const int max_radius_adjust_loops = 50;
    int r_adjust_loop = 0;
    while (true) {
      nSupports = Omega_h::Write<Omega_h::LO>(
        nvertices_target, 0, "number of supports in each target vertex");

      const auto max_radius = Omega_h::get_max(Omega_h::read(radii2));
      pcms::printInfo("INFO: Loop %d: max_radius: %f\n", r_adjust_loop,
                      max_radius);

      adjBasedSearch(mesh, site_dim, target_coords, seed_dim, seeds,
                     supports_ptr, nSupports, supports_idx, radii2, true);
      Kokkos::fence();

      // points without a seed have no supports at any radius
      Omega_h::LO min_supports_found = 0;
      Omega_h::LO max_supports_found = 0;
      Kokkos::parallel_reduce(
        "support_range", nvertices_target,
        OMEGA_H_LAMBDA(const Omega_h::LO i, Omega_h::LO& lo, Omega_h::LO& hi) {
          if (seeds(i) >= 0) {
            lo = nSupports[i] < lo ? nSupports[i] : lo;
            hi = nSupports[i] > hi ? nSupports[i] : hi;
          }
        },
        Kokkos::Min<Omega_h::LO>(min_supports_found),
        Kokkos::Max<Omega_h::LO>(max_supports_found));
      r_adjust_loop++;
      pcms::printInfo(
        "Iter: %d min_nSupports: %d max_nSupports: %d, max_radius %f\n",
        r_adjust_loop, min_supports_found, max_supports_found, max_radius);

      if (within_number_of_support_range(min_supports_found, max_supports_found,
                                         min_req_support,
                                         3 * min_req_support)) {
        break;
      }

      if (r_adjust_loop >= max_radius_adjust_loops) {
        pcms::printInfo(
          "WARNING: radius adaptation hit the %d-loop cap; some targets remain "
          "outside the [%d, %d] support range\n",
          max_radius_adjust_loops, min_req_support, 3 * min_req_support);
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

  adjBasedSearch(mesh, site_dim, target_coords, seed_dim, seeds, supports_ptr,
                 nSupports, supports_idx, radii2, false);

  return SupportResults{read(supports_ptr), read(supports_idx), radii2};
}

} // namespace pcms
