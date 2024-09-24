#ifndef ADJ_SEARCH_HPP
#define ADJ_SEARCH_HPP

#include <Omega_h_macros.h>
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
  Mesh& mesh;

 public:
  FindSupports(Mesh& mesh_) : mesh(mesh_) {};

  void adjBasedSearch(Write<LO>& supports_ptr, Write<LO>& nSupports,
                      Write<LO>& support_idx, Write<Real>& radii2,
                      bool is_build_csr_call);
};

void FindSupports::adjBasedSearch(Write<LO>& supports_ptr, Write<LO>& nSupports,
                                  Write<LO>& support_idx, Write<Real>& radii2,
                                  bool is_build_csr_call) {
  // Mesh Info
  LO min_num_supports = 20;  // TODO: make this an input parameter
  const auto& mesh_coords = mesh.coords();
  const auto& nvertices = mesh.nverts();
  const auto& dim = mesh.dim();
  const auto& nfaces = mesh.nfaces();

  // CSR data structure of adjacent cell information of each vertex in a mesh
  const auto& nodes2faces = mesh.ask_up(VERT, FACE);
  const auto& n2f_ptr = nodes2faces.a2ab;
  const auto& n2f_data = nodes2faces.ab2b;
  const auto& faces2nodes = mesh.ask_down(FACE, VERT).ab2b;

  Write<Real> cell_centroids(
      dim * nfaces, 0,
      "stores coordinates of cell centroid of each tri element");

  parallel_for(
      "calculate the centroid in each tri element", nfaces,
      OMEGA_H_LAMBDA(const LO id) {
        const auto current_el_verts = gather_verts<3>(faces2nodes, id);
        const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
            gather_vectors<3, 2>(mesh_coords, current_el_verts);
        auto centroid = average(current_el_vert_coords);
        int index = dim * id;
        cell_centroids[index] = centroid[0];
        cell_centroids[index + 1] = centroid[1];
      });
  // * Got the adj data and cell centroids

  parallel_for(
      nvertices,  // for each target vertex which is a node for this case
      OMEGA_H_LAMBDA(const LO id) {
        queue queue;
        track visited;
        const LO num_verts_in_dim = dim + 1;
        Real target_coords[max_dim];
        Real support_coords[max_dim];
        Real cutoffDistance = radii2[id];  // squared radii of the supports

        //? copying the target vertex coordinates
        for (LO k = 0; k < dim; ++k) {
          target_coords[k] = mesh_coords[id * dim + k];
        }

        LO start_counter;  // start of support idx for the current target vertex

        // not being done in the first call
        // where the support_idx start for the current target vertex
        if (!is_build_csr_call) {
          start_counter = supports_ptr[id];  // start support location for the
                                             // current target vertex
        }
        LO start_ptr =
            n2f_ptr[id];  // start loc of the adj cells of the target node
        LO end_ptr =
            n2f_ptr[id + 1];  // end loc of the adj cells of the target node

        int count = 0;  // number of supports for the current target vertex
        // Initialize queue by pushing the cells in the neighborhood of the
        // given target point

        for (LO i = start_ptr; i < end_ptr;
             ++i) {  // loop over adj cells to the target vertex
          LO cell_id = n2f_data[i];
          visited.push_back(cell_id);  // cell added to the visited list

          for (LO k = 0; k < dim; ++k) {  // support vertex coordinates are the
                                          // centroid of the cell
            support_coords[k] = cell_centroids[cell_id * dim + k];
          }

          Real dist = calculateDistance(target_coords, support_coords, dim);
          if (dist <= cutoffDistance) {
            count++;
            queue.push_back(cell_id);
            if (!is_build_csr_call) {  // not being done in the first call
              LO idx_count = count - 1;
              support_idx[start_counter + idx_count] =
                  cell_id;  // add the support cell to the support_idx
            }  // end of support_idx check
          }  // end of distance check
        }  // end of loop over adj cells to the target vertex

        // loops over the queued cells from the neighborhood of the target
        // vertex
        while (!queue.isEmpty()) {  // ? can queue be empty?
          LO currentCell = queue.front();
          queue.pop_front();
          LO start = currentCell *
                     num_verts_in_dim;  // start vert id of the current cell
          LO end = start + num_verts_in_dim;  // end vert id of the current cell

          for (LO i = start; i < end;
               ++i) {  // loop over the vertices of the current cell
            LO current_vert_id = faces2nodes[i];
            LO start_ptr_current_vert =
                n2f_ptr[current_vert_id];  // start loc of adj cells to current
                                           // vertex
            LO end_ptr_vert_current_vert =
                n2f_ptr[current_vert_id +
                        1];  // end loc of adj cells to current vertex
            for (LO j = start_ptr_current_vert; j < end_ptr_vert_current_vert;
                 ++j) {  // loop over adj cells to the current vertex
              auto neighbor_cell_index = n2f_data[j];  // current cell

              // check if neighbor index is already in the queue to be checked
              // TODO refactor this into a function

              if (visited.notVisited(neighbor_cell_index)) {
                visited.push_back(neighbor_cell_index);
                for (int k = 0; k < dim; ++k) {
                  support_coords[k] =
                      cell_centroids[neighbor_cell_index * dim + k];
                }

                Real dist =
                    calculateDistance(target_coords, support_coords, dim);

                if (dist <= cutoffDistance) {
                  count++;
                  queue.push_back(neighbor_cell_index);
                  if (!is_build_csr_call) {
                    LO idx_count = count - 1;
                    support_idx[start_counter + idx_count] =
                        neighbor_cell_index;
                  }  // end of support_idx check
                }  // end of distance check
              }  // end of not visited check
            }  // end of loop over adj cells to the current vertex
          }  // end of loop over nodes

        }  // end of while loop

        nSupports[id] = count;
      },  // end of lambda
      "count the number of supports in each target point");
}

struct SupportResults {
  Write<LO> supports_ptr;
  Write<LO> supports_idx;
  Write<Real> radii2;  // squared radii of the supports
};

SupportResults searchNeighbors(Mesh& mesh, Real cutoffDistance) {
  LO min_support = 12;
  SupportResults support;

  FindSupports search(mesh);

  LO nvertices_target = mesh.nverts();

  Write<LO> nSupports(nvertices_target, 0,
                      "number of supports in each target vertex");

  printf("Inside searchNeighbors 1\n");
  support.radii2 = Write<Real>(nvertices_target, cutoffDistance,
                               "squared radii of the supports");
  printf("cut off distance: %f\n", cutoffDistance);
  // this call gets the number of supports for each target vertex: nSupports and
  // the squared radii of the supports (with might be increased in the search to
  // have enough supports)
  int r_adjust_loop = 0;
  while (true) { // until the number of minimum support is met
    // find maximum radius
  Real max_radius = 0.0;
  Kokkos::parallel_reduce("find max radius",
      nvertices_target,
      OMEGA_H_LAMBDA(const LO i, Real& local_max) {
        local_max = (support.radii2[i] > local_max) ? support.radii2[i] : local_max;
      },
      Kokkos::Max<Real>(max_radius));
  printf("Loop %d: max_radius: %f\n", r_adjust_loop, max_radius);
  
  search.adjBasedSearch(support.supports_ptr, nSupports, support.supports_idx,
                        support.radii2, true);

  //printf("Inside searchNeighbors 2\n");
  Kokkos::fence();
  // * find the minimum number of supports
  LO min_nSupports = 0;
  Kokkos::parallel_reduce( "find min number of supports",
      nvertices_target,
      OMEGA_H_LAMBDA(const LO i, LO& local_min) {
        local_min = (nSupports[i] < local_min) ? nSupports[i] : local_min;
      },
      Kokkos::Min<LO>(min_nSupports));

  printf("min_nSupports: %d at loop %d, max_radius %f\n", min_nSupports, r_adjust_loop, max_radius);
  r_adjust_loop++;

  if (min_nSupports >= min_support) {
    break;
  }

  Kokkos::fence();

  // * update radius if nSupport is less that min_support
  parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const LO i) {
        if (nSupports[i] < min_support) {
          Real factor = Real(min_support) / nSupports[i];
          factor = (nSupports[i] == 0 || factor > 3.0) ? 3.0 : factor;
          support.radii2[i] *= factor;
          //nSupports[i] = 0;  // ? might not be needed
        }
      });

  }

  printf("INFO: Took %d loops to adjust the radius\n", r_adjust_loop);

  search.adjBasedSearch(support.supports_ptr, nSupports, support.supports_idx,
                        support.radii2, true);

  // offset array for the supports of each target vertex
  support.supports_ptr =
      Write<LO>(nvertices_target + 1, 0,
                "number of support source vertices in CSR format");

  LO total_supports = 0;

  // get the total number of supports and fill the offset array
  Kokkos::parallel_scan(
      nvertices_target,
      OMEGA_H_LAMBDA(int j, int& update, bool final) {  // what final does?
        // OMEGA_H_CHECK(nSupports[j] >= 15);
        update += nSupports[j];
        if (final) {
          support.supports_ptr[j + 1] = update;
        }
      },
      total_supports);

  printf("Inside searchNeighbors 3\n");
  Kokkos::fence();

  support.supports_idx = Write<LO>(
      total_supports, 0, "index of source supports of each target node");

  printf("total_supports: %d\n", total_supports);
  // second call to get the actual support indices
  // now sizes of support.supports_ptr and support.supports_idx are known and >
  // 0
  search.adjBasedSearch(support.supports_ptr, nSupports, support.supports_idx,
                        support.radii2, false);

  return support;
}

#endif
