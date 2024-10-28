#ifndef ADJ_SEARCH_HPP
#define ADJ_SEARCH_HPP

#include <pcms/point_search.h>

#include "queue_visited.hpp"

using namespace Omega_h;

static constexpr int max_dim = 3;

// TODO change this into span/mdspan
OMEGA_H_INLINE
Real calculateDistance(const Real* p1, const Real* p2, const int dim) {
  Real dx, dy, dz;
  dx = p1[0] - p2[0];
  dy = p1[1] - p2[1];
  if (dim != 3) {
    dz = 0.0;
  } else {
    dz = p1[2] - p2[2];
  }

  return dx * dx + dy * dy + dz * dz;
}

class FindSupports {
 private:
  Mesh& source_mesh;
  Mesh& target_mesh;

 public:
  FindSupports(Mesh& source_mesh_, Mesh& target_mesh_)
      : source_mesh(source_mesh_), target_mesh(target_mesh_) {};

  void adjBasedSearch(const Real& cutoffDistance, const Write<LO>& supports_ptr,
                      Write<LO>& nSupports, Write<LO>& support_idx);
};

void FindSupports::adjBasedSearch(const Real& cutoffDistance,
                                  const Write<LO>& supports_ptr,
                                  Write<LO>& nSupports,
                                  Write<LO>& support_idx) {
  //  Source Mesh Info

  const auto& sourcePoints_coords = source_mesh.coords();
  const auto nvertices_source = source_mesh.nverts();
  const auto dim = source_mesh.dim();

  // Target Mesh Info
  const auto& targetPoints_coords = target_mesh.coords();
  const auto nvertices_target = target_mesh.nverts();

  // CSR data structure of adjacent vertex information of each source vertex
  const auto& vert2vert = source_mesh.ask_star(VERT);
  const auto& v2v_ptr = vert2vert.a2ab;
  const auto& v2v_data = vert2vert.ab2b;

  // CSR data structure of vertex information of each cell in source mesh

  // CSR data structure of adjacent vertex information of each source vertex
  // dim == 2; ask vertices of tri (3 vertices for each tri) & if dim ==3; ask
  // vertices of tetrahedron (4 vertices for each tet)
  const auto& cells2verts = source_mesh.ask_verts_of(dim);

  Kokkos::View<pcms::Real* [2]> target_points("test_points", nvertices_target);

  parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const LO i) {
        target_points(i, 0) = targetPoints_coords[i * dim];
        target_points(i, 1) = targetPoints_coords[i * dim + 1];
      });
  Kokkos::fence();
  pcms::GridPointSearch search_cell(source_mesh, 10, 10);

  // get the cell id for each target point
  auto results = search_cell(target_points);

  parallel_for(
      nvertices_target,
      OMEGA_H_LAMBDA(const LO id) {
        queue queue;
        track visited;

        LO source_cell_id = results(id).tri_id;

        const LO num_verts_in_dim = dim + 1;

        LO start_ptr = source_cell_id * num_verts_in_dim;

        LO end_ptr = start_ptr + num_verts_in_dim;

        Real target_coords[max_dim];

        Real support_coords[max_dim];

        for (LO k = 0; k < dim; ++k) {
          target_coords[k] = target_points(id, k);
        }

        LO start_counter;

        if (support_idx.size() > 0) {
          start_counter = supports_ptr[id];
        }

        int count = 0;
        // Initialize queue by pushing the vertices in the neighborhood of the
        // given target point

        for (LO i = start_ptr; i < end_ptr; ++i) {
          LO vert_id = cells2verts[i];
          visited.push_back(vert_id);

          for (LO k = 0; k < dim; ++k) {
            support_coords[k] = sourcePoints_coords[vert_id * dim + k];
          }

          Real dist = calculateDistance(target_coords, support_coords, dim);
          if (dist <= cutoffDistance) {
            count++;
            queue.push_back(vert_id);
            if (support_idx.size() > 0) {
              LO idx_count = count - 1;
              support_idx[start_counter + idx_count] = vert_id;
            }
          }
        }

        while (!queue.isEmpty()) {
          LO currentVertex = queue.front();
          queue.pop_front();
          LO start = v2v_ptr[currentVertex];
          LO end = v2v_ptr[currentVertex + 1];

          for (LO i = start; i < end; ++i) {
            auto neighborIndex = v2v_data[i];

            // check if neighbor index is already in the queue to be checked
            // TODO refactor this into a function

            if (visited.notVisited(neighborIndex)) {
              visited.push_back(neighborIndex);
              for (int k = 0; k < dim; ++k) {
                support_coords[k] =
                    sourcePoints_coords[neighborIndex * dim + k];
              }

              Real dist = calculateDistance(target_coords, support_coords, dim);

              if (dist <= cutoffDistance) {
                count++;
                queue.push_back(neighborIndex);
                if (support_idx.size() > 0) {
                  LO idx_count = count - 1;
                  support_idx[start_counter + idx_count] = neighborIndex;
                }
              }
            }
          }
        }  // end of while loop

        nSupports[id] = count;
      },
      "count the number of supports in each target point");
}

struct SupportResults {
  Write<LO> supports_ptr;
  Write<LO> supports_idx;
};

SupportResults searchNeighbors(Mesh& source_mesh, Mesh& target_mesh,
                               Real& cutoffDistance) {
  SupportResults support;

  FindSupports search(source_mesh, target_mesh);

  LO nvertices_source = source_mesh.nverts();

  LO nvertices_target = target_mesh.nverts();

  Write<LO> nSupports(nvertices_target, 0,
                      "number of supports in each target vertex");

  search.adjBasedSearch(cutoffDistance, support.supports_ptr, nSupports,
                        support.supports_idx);

  Kokkos::fence();

  support.supports_ptr =
      Write<LO>(nvertices_target + 1, 0,
                "number of support source vertices in CSR format");

  LO total_supports = 0;

  Kokkos::parallel_scan(
      nvertices_target,
      OMEGA_H_LAMBDA(int j, int& update, bool final) {
        update += nSupports[j];
        if (final) {
          support.supports_ptr[j + 1] = update;
        }
      },
      total_supports);

  Kokkos::fence();

  support.supports_idx = Write<LO>(
      total_supports, 0, "index of source supports of each target node");

  search.adjBasedSearch(cutoffDistance, support.supports_ptr, nSupports,
                        support.supports_idx);

  return support;
}

#endif
