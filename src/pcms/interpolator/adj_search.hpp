#ifndef ADJ_SEARCH_HPP
#define ADJ_SEARCH_HPP

#include <pcms/point_search.h>

#include "queue_visited.hpp"

using namespace Omega_h;

static constexpr int max_dim = 3;

// TODO change this into span/mdspan
OMEGA_H_INLINE
Real calculateDistance(const Real* p1, const Real* p2, const int dim)
{
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

void checkTargetPoints(
  const Kokkos::View<pcms::GridPointSearch::Result*>& results)
{
  Kokkos::fence();
  printf("INFO: Checking target points...\n");
  auto check_target_points = OMEGA_H_LAMBDA(LO i)
  {
    if (results(i).tri_id < 0) {
      OMEGA_H_CHECK_PRINTF(results(i).tri_id >= 0,
                           "ERROR: Source cell id not found for target %d\n",
                           i);
      printf("%d, ", i);
    }
  };
  parallel_for(results.size(), check_target_points, "check_target_points");
  Kokkos::fence();
  printf("\n");
}

void printSupportsForTarget(const LO target_id, const Write<LO>& supports_ptr,
                            const Write<LO>& nSupports,
                            const Write<LO>& support_idx)
{
  parallel_for(
    nSupports.size(), OMEGA_H_LAMBDA(const LO id) {
      if (id == target_id) {
        LO start_ptr = supports_ptr[id]; // start support location for the
                                         // current target vertex
        LO end_ptr =
          supports_ptr[id +
                       1]; // end support location for the current target vertex
        printf("Target vertex: %d\n with %d num supports: nSupports[id]=%d", id,
               end_ptr - start_ptr, nSupports[id]);
        for (LO i = start_ptr; i < end_ptr;
             ++i) { // loop over adj cells to the target vertex
          LO cell_id = support_idx[i];
          printf(", %d", cell_id);
        }
        printf("\n");
      }
    });
}

class FindSupports
{
private:
  Mesh& source_mesh;
  Mesh& target_mesh; // TODO it's null when one mesh is used

public:
  FindSupports(Mesh& source_mesh_, Mesh& target_mesh_)
    : source_mesh(source_mesh_), target_mesh(target_mesh_) {};

  FindSupports(Mesh& mesh_) : source_mesh(mesh_), target_mesh(mesh_) {};

  void adjBasedSearch(Write<LO>& supports_ptr, Write<LO>& nSupports,
                      Write<LO>& support_idx, Write<Real>& radii2,
                      bool is_build_csr_call);

  void adjBasedSearchCentroidNodes(Write<LO>& supports_ptr,
                                   Write<LO>& nSupports, Write<LO>& support_idx,
                                   Write<Real>& radii2, bool is_build_csr_call);
};

void FindSupports::adjBasedSearch(Write<LO>& supports_ptr, Write<LO>& nSupports,
                                  Write<LO>& support_idx, Write<Real>& radii2,
                                  bool is_build_csr_call)
{

  //  Source Mesh Info

  const auto& sourcePoints_coords = source_mesh.coords();
  const auto nvertices_source = source_mesh.nverts();
  const auto dim = source_mesh.dim();

  // Target Mesh Info
  const auto& targetPoints_coords = target_mesh.coords();
  const auto nvertices_target = target_mesh.nverts();
  OMEGA_H_CHECK(radii2.size() == nvertices_target);

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
  pcms::GridPointSearch search_cell(source_mesh, 1000, 1000);

  // get the cell id for each target point
  auto results = search_cell(target_points);

  checkTargetPoints(results);

  parallel_for(
    nvertices_target,
    OMEGA_H_LAMBDA(const LO id) {
      queue queue;
      track visited;
      Real cutoffDistance = radii2[id];

      LO source_cell_id = results(id).tri_id;
      OMEGA_H_CHECK_PRINTF(
        source_cell_id >= 0,
        "ERROR: Source cell id not found for target %d (%f,%f)\n", id,
        target_points(id, 0), target_points(id, 1));

      const LO num_verts_in_dim = dim + 1;

      LO start_ptr = source_cell_id * num_verts_in_dim;

      LO end_ptr = start_ptr + num_verts_in_dim;

      Real target_coords[max_dim];

      Real support_coords[max_dim];

      for (LO k = 0; k < dim; ++k) {
        target_coords[k] = target_points(id, k);
      }

      LO start_counter;

      if (!is_build_csr_call) {
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
          if (count >= 500) {
            printf(
              "Warning: count exceeds 500 for target %d with %d supports\n", id,
              end_ptr - start_ptr);
            printf("Warning: Target %d: coors: (%f, %f) and support %d: "
                   "coords: (%f, %f)\n",
                   id, target_coords[0], target_coords[1], vert_id,
                   support_coords[0], support_coords[1]);
          }
          queue.push_back(vert_id);
          if (!is_build_csr_call) {
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
              support_coords[k] = sourcePoints_coords[neighborIndex * dim + k];
            }

            Real dist = calculateDistance(target_coords, support_coords, dim);

            if (dist <= cutoffDistance) {
              count++;
              if (count >= 500) {
                // printf("Warning qne: count exceeds 500 for target %d with
                // start %d and end %d radius2 %f adding neighbor %d\n", id,
                // start, end, cutoffDistance, neighborIndex);
              }
              queue.push_back(neighborIndex);
              if (!is_build_csr_call) {
                LO idx_count = count - 1;
                support_idx[start_counter + idx_count] = neighborIndex;
              }
            }
          }
        }
      } // end of while loop

      nSupports[id] = count;
    },
    "count the number of supports in each target point");
  if (is_build_csr_call == false) {
    // printSupportsForTarget(2057,  supports_ptr, nSupports, support_idx);
  }
}

void FindSupports::adjBasedSearchCentroidNodes(Write<LO>& supports_ptr,
                                               Write<LO>& nSupports,
                                               Write<LO>& support_idx,
                                               Write<Real>& radii2,
                                               bool is_build_csr_call)
{
  // Mesh Info
  const auto& mesh_coords = source_mesh.coords();
  const auto& nvertices = source_mesh.nverts();
  const auto& dim = source_mesh.dim();
  const auto& nfaces = source_mesh.nfaces();

  // CSR data structure of adjacent cell information of each vertex in a
  // source_mesh
  const auto& nodes2faces = source_mesh.ask_up(VERT, FACE);
  const auto& n2f_ptr = nodes2faces.a2ab;
  const auto& n2f_data = nodes2faces.ab2b;
  const auto& faces2nodes = source_mesh.ask_down(FACE, VERT).ab2b;

  Write<Real> cell_centroids(
    dim * nfaces, 0, "stores coordinates of cell centroid of each tri element");

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
    nvertices, // for each target vertex which is a node for this case
    OMEGA_H_LAMBDA(const LO id) {
      queue queue;
      track visited;
      const LO num_verts_in_dim = dim + 1;
      Real target_coords[max_dim];
      Real support_coords[max_dim];
      Real cutoffDistance = radii2[id]; // squared radii of the supports

      //? copying the target vertex coordinates
      for (LO k = 0; k < dim; ++k) {
        target_coords[k] = mesh_coords[id * dim + k];
      }

      LO start_counter; // start of support idx for the current target vertex

      // not being done in the first call
      // where the support_idx start for the current target vertex
      if (!is_build_csr_call) {
        start_counter = supports_ptr[id]; // start support location for the
                                          // current target vertex
      }
      LO start_ptr =
        n2f_ptr[id]; // start loc of the adj cells of the target node
      LO end_ptr =
        n2f_ptr[id + 1]; // end loc of the adj cells of the target node

      int count = 0; // number of supports for the current target vertex
      // Initialize queue by pushing the cells in the neighborhood of the
      // given target point

      for (LO i = start_ptr; i < end_ptr;
           ++i) { // loop over adj cells to the target vertex
        LO cell_id = n2f_data[i];
        visited.push_back(cell_id); // cell added to the visited list

        for (LO k = 0; k < dim; ++k) { // support vertex coordinates are the
                                       // centroid of the cell
          support_coords[k] = cell_centroids[cell_id * dim + k];
        }

        Real dist = calculateDistance(target_coords, support_coords, dim);
        if (dist <= cutoffDistance) {
          count++;
          queue.push_back(cell_id);
          if (!is_build_csr_call) { // not being done in the first call
            LO idx_count = count - 1;
            support_idx[start_counter + idx_count] =
              cell_id; // add the support cell to the support_idx
          } // end of support_idx check
        } // end of distance check
      } // end of loop over adj cells to the target vertex

      // loops over the queued cells from the neighborhood of the target
      // vertex
      while (!queue.isEmpty()) { // ? can queue be empty?
        LO currentCell = queue.front();
        queue.pop_front();
        LO start =
          currentCell * num_verts_in_dim;  // start vert id of the current cell
        LO end = start + num_verts_in_dim; // end vert id of the current cell

        for (LO i = start; i < end;
             ++i) { // loop over the vertices of the current cell
          LO current_vert_id = faces2nodes[i];
          LO start_ptr_current_vert =
            n2f_ptr[current_vert_id]; // start loc of adj cells to current
                                      // vertex
          LO end_ptr_vert_current_vert =
            n2f_ptr[current_vert_id +
                    1]; // end loc of adj cells to current vertex
          for (LO j = start_ptr_current_vert; j < end_ptr_vert_current_vert;
               ++j) { // loop over adj cells to the current vertex
            auto neighbor_cell_index = n2f_data[j]; // current cell

            // check if neighbor index is already in the queue to be checked
            // TODO refactor this into a function

            if (visited.notVisited(neighbor_cell_index)) {
              visited.push_back(neighbor_cell_index);
              for (int k = 0; k < dim; ++k) {
                support_coords[k] =
                  cell_centroids[neighbor_cell_index * dim + k];
              }

              Real dist = calculateDistance(target_coords, support_coords, dim);

              if (dist <= cutoffDistance) {
                count++;
                queue.push_back(neighbor_cell_index);
                if (!is_build_csr_call) {
                  LO idx_count = count - 1;
                  support_idx[start_counter + idx_count] = neighbor_cell_index;
                } // end of support_idx check
              } // end of distance check
            } // end of not visited check
          } // end of loop over adj cells to the current vertex
        } // end of loop over nodes

      } // end of while loop

      nSupports[id] = count;
    }, // end of lambda
    "count the number of supports in each target point");

  if (is_build_csr_call == false) {
    // printSupportsForTarget(2057,  supports_ptr, nSupports, support_idx);
  }
}

struct SupportResults
{
  Write<LO> supports_ptr;
  Write<LO> supports_idx;
  Write<Real> radii2; // squared radii of the supports
};

SupportResults searchNeighbors(Mesh& source_mesh, Mesh& target_mesh,
                               Real& cutoffDistance, LO min_req_support = 12,
                               bool adapt_radius = true)
{
  SupportResults support;

  FindSupports search(source_mesh, target_mesh);

  LO nvertices_source = source_mesh.nverts();

  LO nvertices_target = target_mesh.nverts();

  Write<LO> nSupports(nvertices_target, 0,
                      "number of supports in each target vertex");
  printf("INFO: Cut off distance: %f\n", cutoffDistance);
  Write<Real> radii2 = Write<Real>(nvertices_target, cutoffDistance,
                                   "squared radii of the supports");

  if (!adapt_radius) {
    printf("INFO: Fixed radius search... \n");
    search.adjBasedSearch(support.supports_ptr, nSupports, support.supports_idx,
                          radii2, true);
  } else {
    printf("INFO: Adaptive radius search... \n");
    int r_adjust_loop = 0;
    while (true) {
      nSupports = Write<LO>(nvertices_target, 0,
                            "number of supports in each target vertex");
      // find maximum radius
      Real max_radius = 0.0;
      Kokkos::parallel_reduce(
        "find max radius", nvertices_target,
        OMEGA_H_LAMBDA(const LO i, Real& local_max) {
          local_max = (radii2[i] > local_max) ? radii2[i] : local_max;
        },
        Kokkos::Max<Real>(max_radius));
      printf("INFO: Loop %d: max_radius: %f\n", r_adjust_loop, max_radius);

      SupportResults support; // create support every time to avoid complexity

      Kokkos::fence();
      search.adjBasedSearch(support.supports_ptr, nSupports,
                            support.supports_idx, radii2, true);

      Kokkos::fence();
      //* find the minimum number of supports
      LO min_supports_found = 0;
      Kokkos::Min<LO> min_reducer(min_supports_found);
      Kokkos::parallel_reduce(
        "find min number of supports", nvertices_target,
        OMEGA_H_LAMBDA(const LO i, LO& local_min) {
          min_reducer.join(local_min, nSupports[i]);
        },
        min_reducer);

      printf("INFO: min_supports_found: %d at loop %d, max_radius %f\n",
             min_supports_found, r_adjust_loop, max_radius);
      r_adjust_loop++;

      Kokkos::fence();
      if (min_supports_found >= min_req_support) {
        break;
      }

      Kokkos::fence();

      // * update radius if nSupport is less that min_req_support
      parallel_for(
        nvertices_target, OMEGA_H_LAMBDA(const LO i) {
          if (nSupports[i] < min_req_support) {
            Real factor = Real(min_req_support) / Real(nSupports[i]);
            factor = (factor > 1.1 || nSupports[i] == 0) ? 1.1 : factor;
            radii2[i] = radii2[i] * factor;
          }
          // nSupports[i] = 0; // reset the number of supports ? not sure if
          // needed
        });
      Kokkos::fence();
    }

    printf("INFO: Took %d loops to adjust the radius\n", r_adjust_loop);
  }

  support.supports_ptr = Write<LO>(
    nvertices_target + 1, 0, "number of support source vertices in CSR format");

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

  search.adjBasedSearch(support.supports_ptr, nSupports, support.supports_idx,
                        radii2, false);

  support.radii2 = radii2;

  target_mesh.add_tag<Real>(VERT, "radii2", 1, support.radii2);

  return support;
}

SupportResults searchNeighbors(Mesh& mesh, Real cutoffDistance,
                               LO min_support = 12, bool adapt_radius = true)
{
  SupportResults support;

  FindSupports search(mesh);

  LO nvertices_target = mesh.nverts();

  Write<LO> nSupports(nvertices_target, 0,
                      "number of supports in each target vertex");

  printf("INFO: Inside searchNeighbors 1\n");
  Write<Real> radii2 = Write<Real>(nvertices_target, cutoffDistance,
                                   "squared radii of the supports");
  printf("INFO: Cutoff distance: %f\n", cutoffDistance);
  // this call gets the number of supports for each target vertex: nSupports and
  // the squared radii of the supports (with might be increased in the search to
  // have enough supports)
  int r_adjust_loop = 0;
  while (true) { // until the number of minimum support is met
    // find maximum radius
    Real max_radius = 0.0;
    Kokkos::parallel_reduce(
      "find max radius", nvertices_target,
      OMEGA_H_LAMBDA(const LO i, Real& local_max) {
        local_max = (radii2[i] > local_max) ? radii2[i] : local_max;
      },
      Kokkos::Max<Real>(max_radius));
    printf("INFO: Loop %d: max_radius: %f\n", r_adjust_loop, max_radius);

    nSupports = Write<LO>(nvertices_target, 0,
                          "number of supports in each target vertex");
    SupportResults support; // create support every time to avoid complexity

    search.adjBasedSearchCentroidNodes(support.supports_ptr, nSupports,
                                       support.supports_idx, radii2, true);

    Kokkos::fence();
    //* find the minimum number of supports
    LO min_nSupports = 0;
    Kokkos::parallel_reduce(
      "find min number of supports", nvertices_target,
      OMEGA_H_LAMBDA(const LO i, LO& local_min) {
        local_min = (nSupports[i] < local_min) ? nSupports[i] : local_min;
      },
      Kokkos::Min<LO>(min_nSupports));

    printf("min_nSupports: %d at loop %d, max_radius %f\n", min_nSupports,
           r_adjust_loop, max_radius);
    r_adjust_loop++;

    if (min_nSupports >= min_support) {
      break;
    }

    Kokkos::fence();

    // * update radius if nSupport is less that min_support
    parallel_for(
      nvertices_target, OMEGA_H_LAMBDA(const LO i) {
        if (nSupports[i] < min_support) {
          Real factor = Real(min_support) / Real(nSupports[i]);
          factor = (nSupports[i] == 0 || factor > 1.5) ? 1.5 : factor;
          radii2[i] *= factor;
        }
        nSupports[i] = 0; // ? might not be needed
      });
  } // while loop

  printf("INFO: Took %d loops to adjust the radius\n", r_adjust_loop);

  // search.adjBasedSearch(support.supports_ptr, nSupports,
  // support.supports_idx,
  //                       radii2, true);

  // offset array for the supports of each target vertex
  support.supports_ptr = Write<LO>(
    nvertices_target + 1, 0, "number of support source vertices in CSR format");

  LO total_supports = 0;

  // get the total number of supports and fill the offset array
  Kokkos::parallel_scan(
    nvertices_target,
    OMEGA_H_LAMBDA(int j, int& update, bool final) { // what final does?
      // OMEGA_H_CHECK(nSupports[j] >= 15);
      update += nSupports[j];
      if (final) {
        support.supports_ptr[j + 1] = update;
      }
    },
    total_supports);

  printf("INFO: Inside searchNeighbors 3\n");
  Kokkos::fence();

  support.supports_idx = Write<LO>(
    total_supports, 0, "index of source supports of each target node");

  printf("INFO: Total_supports: %d\n", total_supports);
  // second call to get the actual support indices
  // now sizes of support.supports_ptr and support.supports_idx are known and >
  // 0
  search.adjBasedSearchCentroidNodes(support.supports_ptr, nSupports,
                                     support.supports_idx, radii2, false);
  support.radii2 = radii2;

  mesh.add_tag<Omega_h::Real>(VERT, "support_radius", 1, radii2);
  return support;
}

#endif
