#ifndef ADJ_SEARCH_HPP
#define ADJ_SEARCH_HPP

#include <pcms/point_search.h>

#include "queue_visited.hpp"

static constexpr int max_dim = 3;

// TODO change this into span/mdspan
OMEGA_H_INLINE
Omega_h::Real calculateDistance(const Omega_h::Real* p1,
                                const Omega_h::Real* p2, const int dim)
{
  Omega_h::Real dx, dy, dz;
  dx = p1[0] - p2[0];
  dy = p1[1] - p2[1];
  if (dim != 3) {
    dz = 0.0;
  } else {
    dz = p1[2] - p2[2];
  }
  return dx * dx + dy * dy + dz * dz;
}

inline void checkTargetPoints(
  const Kokkos::View<pcms::GridPointSearch::Result*>& results)
{
  Kokkos::fence();
  printf("INFO: Checking target points...\n");
  auto check_target_points = OMEGA_H_LAMBDA(Omega_h::LO i)
  {
    if (results(i).tri_id < 0) {
      OMEGA_H_CHECK_PRINTF(results(i).tri_id >= 0,
                           "ERROR: Source cell id not found for target %d\n",
                           i);
      printf("%d, ", i);
    }
  };
  Omega_h::parallel_for(results.size(), check_target_points,
                        "check_target_points");
  Kokkos::fence();
  printf("\n");
}

inline void printSupportsForTarget(
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

class FindSupports
{
private:
  Omega_h::Mesh& source_mesh;
  Omega_h::Mesh& target_mesh; // TODO it's null when one mesh is used

public:
  FindSupports(Omega_h::Mesh& source_mesh_, Omega_h::Mesh& target_mesh_)
    : source_mesh(source_mesh_), target_mesh(target_mesh_) {};

  FindSupports(Omega_h::Mesh& mesh_)
    : source_mesh(mesh_), target_mesh(mesh_) {};

  void adjBasedSearch(Omega_h::Write<Omega_h::LO>& supports_ptr,
                      Omega_h::Write<Omega_h::LO>& nSupports,
                      Omega_h::Write<Omega_h::LO>& support_idx,
                      Omega_h::Write<Omega_h::Real>& radii2,
                      bool is_build_csr_call);

  void adjBasedSearchCentroidNodes(Omega_h::Write<Omega_h::LO>& supports_ptr,
                                   Omega_h::Write<Omega_h::LO>& nSupports,
                                   Omega_h::Write<Omega_h::LO>& support_idx,
                                   Omega_h::Write<Omega_h::Real>& radii2,
                                   bool is_build_csr_call);
};

inline void FindSupports::adjBasedSearch(
  Omega_h::Write<Omega_h::LO>& supports_ptr,
  Omega_h::Write<Omega_h::LO>& nSupports,
  Omega_h::Write<Omega_h::LO>& support_idx,
  Omega_h::Write<Omega_h::Real>& radii2, bool is_build_csr_call)
{

  const auto& sourcePoints_coords = source_mesh.coords();
  const auto dim = source_mesh.dim();

  const auto& targetPoints_coords = target_mesh.coords();
  const auto nvertices_target = target_mesh.nverts();
  OMEGA_H_CHECK(radii2.size() == nvertices_target);

  const auto& vert2vert = source_mesh.ask_star(Omega_h::VERT);
  const auto& v2v_ptr = vert2vert.a2ab;
  const auto& v2v_data = vert2vert.ab2b;
  const auto& cells2verts = source_mesh.ask_verts_of(dim);

  Kokkos::View<Omega_h::Real* [2]> target_points("test_points",
                                                 nvertices_target);
  Omega_h::parallel_for(
    nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
      target_points(i, 0) = targetPoints_coords[i * dim];
      target_points(i, 1) = targetPoints_coords[i * dim + 1];
    });
  Kokkos::fence();

  pcms::GridPointSearch search_cell(source_mesh, 10, 10);
  auto results = search_cell(target_points);
  checkTargetPoints(results);

  Omega_h::parallel_for(
    nvertices_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;
      Omega_h::Real cutoffDistance = radii2[id];

      Omega_h::LO source_cell_id = results(id).tri_id;
      OMEGA_H_CHECK_PRINTF(
        source_cell_id >= 0,
        "ERROR: Source cell id not found for target %d (%f,%f)\n", id,
        target_points(id, 0), target_points(id, 1));

      const Omega_h::LO num_verts_in_dim = dim + 1;
      Omega_h::LO start_ptr = source_cell_id * num_verts_in_dim;
      Omega_h::LO end_ptr = start_ptr + num_verts_in_dim;
      Omega_h::Real target_coords[max_dim];
      Omega_h::Real support_coords[max_dim];

      for (Omega_h::LO k = 0; k < dim; ++k) {
        target_coords[k] = target_points(id, k);
      }

      Omega_h::LO start_counter;
      if (!is_build_csr_call) {
        start_counter = supports_ptr[id];
      }

      // * Method:
      // 1. Get the vertices of the source cell (source cell is the cell in the
      // source mesh in which the target point lies): done above
      // 2. Using those 3 vertices, get the adjacent vertices of those 3
      // vertices and go on until the queue is empty
      // 3. Already visited vertices are stored in visited and the vertices to
      // be checked (dist < cutoff) are stored in the queue
      // 4. If not CSR building call, store the support vertices in support_idx
      // * Method

      int count = 0;
      for (Omega_h::LO i = start_ptr; i < end_ptr; ++i) {
        Omega_h::LO vert_id = cells2verts[i];
        visited.push_back(vert_id);

        for (Omega_h::LO k = 0; k < dim; ++k) {
          support_coords[k] = sourcePoints_coords[vert_id * dim + k];
        }

        Omega_h::Real dist =
          calculateDistance(target_coords, support_coords, dim);
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
            Omega_h::LO idx_count = count - 1;
            support_idx[start_counter + idx_count] = vert_id;
          }
        }
      }

      while (!queue.isEmpty()) {
        Omega_h::LO currentVertex = queue.front();
        queue.pop_front();
        Omega_h::LO start = v2v_ptr[currentVertex];
        Omega_h::LO end = v2v_ptr[currentVertex + 1];

        for (Omega_h::LO i = start; i < end; ++i) {
          auto neighborIndex = v2v_data[i];

          // check if neighbor index is already in the queue to be checked
          // TODO refactor this into a function

          if (visited.notVisited(neighborIndex)) {
            visited.push_back(neighborIndex);
            for (int k = 0; k < dim; ++k) {
              support_coords[k] = sourcePoints_coords[neighborIndex * dim + k];
            }

            Omega_h::Real dist =
              calculateDistance(target_coords, support_coords, dim);

            if (dist <= cutoffDistance) {
              count++;
              if (count >= 500) {
                printf("Warning: count exceeds 500 for target %d with start %d "
                       "and end %d radius2 %f adding neighbor %d\n",
                       id, start, end, cutoffDistance, neighborIndex);
              }
              queue.push_back(neighborIndex);
              if (!is_build_csr_call) {
                Omega_h::LO idx_count = count - 1;
                support_idx[start_counter + idx_count] = neighborIndex;
              }
            }
          }
        }
      } // end of while loop

      nSupports[id] = count;
    }, // lambda
    "count the number of supports in each target point");
  if (is_build_csr_call == false) {
    // printSupportsForTarget(2057,  supports_ptr, nSupports, support_idx);
  }
}

inline void FindSupports::adjBasedSearchCentroidNodes(
  Omega_h::Write<Omega_h::LO>& supports_ptr,
  Omega_h::Write<Omega_h::LO>& nSupports,
  Omega_h::Write<Omega_h::LO>& support_idx,
  Omega_h::Write<Omega_h::Real>& radii2, bool is_build_csr_call)
{
  // Mesh Info
  const auto& mesh_coords = source_mesh.coords();
  const auto& nvertices = source_mesh.nverts();
  const auto& dim = source_mesh.dim();
  const auto& nfaces = source_mesh.nfaces();

  const auto& nodes2faces = source_mesh.ask_up(Omega_h::VERT, Omega_h::FACE);
  const auto& n2f_ptr = nodes2faces.a2ab;
  const auto& n2f_data = nodes2faces.ab2b;
  const auto& faces2nodes =
    source_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  Omega_h::Write<Omega_h::Real> cell_centroids(
    dim * nfaces, 0, "stores coordinates of cell centroid of each tri element");

  Omega_h::parallel_for(
    "calculate the centroid in each tri element", nfaces,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      const auto current_el_verts = Omega_h::gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
        Omega_h::gather_vectors<3, 2>(mesh_coords, current_el_verts);
      auto centroid = Omega_h::average(current_el_vert_coords);
      int index = dim * id;
      cell_centroids[index] = centroid[0];
      cell_centroids[index + 1] = centroid[1];
    });
  // * Got the adj data and cell centroids

  Omega_h::parallel_for(
    nvertices,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;
      const Omega_h::LO num_verts_in_dim = dim + 1;
      Omega_h::Real target_coords[max_dim];
      Omega_h::Real support_coords[max_dim];
      Omega_h::Real cutoffDistance = radii2[id];

      //? copying the target vertex coordinates
      for (Omega_h::LO k = 0; k < dim; ++k) {
        target_coords[k] = mesh_coords[id * dim + k];
      }

      Omega_h::LO start_counter;
      if (!is_build_csr_call) {
        start_counter = supports_ptr[id];
      }
      Omega_h::LO start_ptr = n2f_ptr[id];
      Omega_h::LO end_ptr = n2f_ptr[id + 1];

      int count = 0;
      for (Omega_h::LO i = start_ptr; i < end_ptr; ++i) {
        Omega_h::LO cell_id = n2f_data[i];
        visited.push_back(cell_id);

        for (Omega_h::LO k = 0; k < dim; ++k) {
          support_coords[k] = cell_centroids[cell_id * dim + k];
        }

        Omega_h::Real dist =
          calculateDistance(target_coords, support_coords, dim);
        if (dist <= cutoffDistance) {
          count++;
          queue.push_back(cell_id);
          if (!is_build_csr_call) {
            Omega_h::LO idx_count = count - 1;
            support_idx[start_counter + idx_count] = cell_id;
          }
        }
      }

      while (!queue.isEmpty()) { // ? can queue be empty?
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

            if (visited.notVisited(neighbor_cell_index)) {
              visited.push_back(neighbor_cell_index);
              for (int k = 0; k < dim; ++k) {
                support_coords[k] =
                  cell_centroids[neighbor_cell_index * dim + k];
              }

              Omega_h::Real dist =
                calculateDistance(target_coords, support_coords, dim);

              if (dist <= cutoffDistance) {
                count++;
                queue.push_back(neighbor_cell_index);
                if (!is_build_csr_call) {
                  Omega_h::LO idx_count = count - 1;
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
  Omega_h::LOs supports_ptr;
  Omega_h::LOs supports_idx;
  Omega_h::Write<Omega_h::Real> radii2;
};

inline SupportResults searchNeighbors(Omega_h::Mesh& source_mesh,
                                      Omega_h::Mesh& target_mesh,
                                      Omega_h::Real& cutoffDistance,
                                      Omega_h::LO min_req_support = 12,
                                      bool adapt_radius = true)
{
  FindSupports search(source_mesh, target_mesh);
  Omega_h::LO nvertices_target = target_mesh.nverts();

  Omega_h::Write<Omega_h::LO> nSupports(
    nvertices_target, 0, "number of supports in each target vertex");
  printf("INFO: Cut off distance: %f\n", cutoffDistance);
  Omega_h::Write<Omega_h::Real> radii2 = Omega_h::Write<Omega_h::Real>(
    nvertices_target, cutoffDistance, "squared radii of the supports");

  Omega_h::Write<Omega_h::LO> supports_ptr;
  Omega_h::Write<Omega_h::LO> supports_idx;

  if (!adapt_radius) {
    printf("INFO: Fixed radius search *(disregarding required minimum "
           "support)*... \n");
    search.adjBasedSearch(supports_ptr, nSupports, supports_idx, radii2, true);
  } else {
    printf("INFO: Adaptive radius search... \n");
    int r_adjust_loop = 0;
    while (true) {
      nSupports = Omega_h::Write<Omega_h::LO>(
        nvertices_target, 0, "number of supports in each target vertex");

      Omega_h::Real max_radius = 0.0;
      Kokkos::parallel_reduce(
        "find max radius", nvertices_target,
        OMEGA_H_LAMBDA(const Omega_h::LO i, Omega_h::Real& local_max) {
          local_max = (radii2[i] > local_max) ? radii2[i] : local_max;
        },
        Kokkos::Max<Omega_h::Real>(max_radius));
      printf("INFO: Loop %d: max_radius: %f\n", r_adjust_loop, max_radius);

      // create storage every time to avoid complexity
      // FIXME avoid repeated dynamic allocation
      Omega_h::Write<Omega_h::LO> temp_supports_ptr;
      Omega_h::Write<Omega_h::LO> temp_supports_idx;
      Kokkos::fence();
      search.adjBasedSearch(temp_supports_ptr, nSupports, temp_supports_idx, radii2,
                            true);
      Kokkos::fence();

      Omega_h::LO min_supports_found = 0;
      Kokkos::Min<Omega_h::LO> min_reducer(min_supports_found);
      Kokkos::parallel_reduce(
        "find min number of supports", nvertices_target,
        OMEGA_H_LAMBDA(const Omega_h::LO i, Omega_h::LO& local_min) {
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
      Omega_h::parallel_for(
        nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
          if (nSupports[i] < min_req_support) {
            Omega_h::Real factor =
              Omega_h::Real(min_req_support) / Omega_h::Real(nSupports[i]);
            factor = (factor > 1.1 || nSupports[i] == 0) ? 1.1 : factor;
            radii2[i] = radii2[i] * factor;
          }
        });
      Kokkos::fence();
    }

    printf("INFO: Took %d loops to adjust the radius\n", r_adjust_loop);
  }

  supports_ptr = Omega_h::Write<Omega_h::LO>(
    nvertices_target + 1, 0, "number of support source vertices in CSR format");

  Omega_h::LO total_supports = 0;

  Kokkos::parallel_scan(
    nvertices_target,
    OMEGA_H_LAMBDA(int j, int& update, bool final) {
      update += nSupports[j];
      if (final) {
        supports_ptr[j + 1] = update;
      }
    },
    total_supports);

  Kokkos::fence();

  supports_idx = Omega_h::Write<Omega_h::LO>(
    total_supports, 0, "index of source supports of each target node");

  search.adjBasedSearch(supports_ptr, nSupports, supports_idx, radii2, false);

  target_mesh.add_tag<Omega_h::Real>(Omega_h::VERT, "radii2", 1, radii2);
  return SupportResults{read(supports_ptr), read(supports_idx), radii2};
}

inline SupportResults searchNeighbors(Omega_h::Mesh& mesh,
                                      Omega_h::Real cutoffDistance,
                                      Omega_h::LO min_support = 12,
                                      bool adapt_radius = true)
{
  Omega_h::Write<Omega_h::LO> supports_ptr;
  Omega_h::Write<Omega_h::LO> supports_idx;
  FindSupports search(mesh);
  Omega_h::LO nvertices_target = mesh.nverts();
  Omega_h::Write<Omega_h::LO> nSupports(
    nvertices_target, 0, "number of supports in each target vertex");

  printf("INFO: Inside searchNeighbors 1\n");
  Omega_h::Write<Omega_h::Real> radii2 = Omega_h::Write<Omega_h::Real>(
    nvertices_target, cutoffDistance, "squared radii of the supports");
  printf("INFO: Cutoff distance: %f\n", cutoffDistance);

  if (!adapt_radius) {
    printf("INFO: Fixed radius search *(disregarding required minimum "
           "support)* ... \n");
    search.adjBasedSearch(supports_ptr, nSupports, supports_idx, radii2, true);
  } else {
    printf("INFO: Adaptive radius search... \n");
    int r_adjust_loop = 0;
    while (true) { // until the number of minimum support is met
      Omega_h::Real max_radius = 0.0;
      Kokkos::parallel_reduce(
        "find max radius", nvertices_target,
        OMEGA_H_LAMBDA(const Omega_h::LO i, Omega_h::Real& local_max) {
          local_max = (radii2[i] > local_max) ? radii2[i] : local_max;
        },
        Kokkos::Max<Omega_h::Real>(max_radius));
      printf("INFO: Loop %d: max_radius: %f\n", r_adjust_loop, max_radius);

      nSupports = Omega_h::Write<Omega_h::LO>(
        nvertices_target, 0, "number of supports in each target vertex");
      SupportResults support; // create support every time to avoid complexity
      search.adjBasedSearchCentroidNodes(supports_ptr, nSupports, supports_idx,
                                         radii2, true);

      Kokkos::fence();
      Omega_h::LO min_nSupports = 0;
      Kokkos::parallel_reduce(
        "find min number of supports", nvertices_target,
        OMEGA_H_LAMBDA(const Omega_h::LO i, Omega_h::LO& local_min) {
          local_min = (nSupports[i] < local_min) ? nSupports[i] : local_min;
        },
        Kokkos::Min<Omega_h::LO>(min_nSupports));

      printf("min_nSupports: %d at loop %d, max_radius %f\n", min_nSupports,
             r_adjust_loop, max_radius);
      r_adjust_loop++;

      if (min_nSupports >= min_support) {
        break;
      }

      Kokkos::fence();
      Omega_h::parallel_for(
        nvertices_target, OMEGA_H_LAMBDA(const Omega_h::LO i) {
          if (nSupports[i] < min_support) {
            Omega_h::Real factor =
              Omega_h::Real(min_support) / Omega_h::Real(nSupports[i]);
            factor = (nSupports[i] == 0 || factor > 1.5) ? 1.5 : factor;
            radii2[i] *= factor;
          }
          nSupports[i] = 0; // ? might not be needed
        });
    } // while loop
    printf("INFO: Took %d loops to adjust the radius\n", r_adjust_loop);
  } // adaptive radius search

  // offset array for the supports of each target vertex
  supports_ptr = Omega_h::Write<Omega_h::LO>(
    nvertices_target + 1, 0, "number of support source vertices in CSR format");

  Omega_h::LO total_supports = 0;
  Kokkos::parallel_scan(
    nvertices_target,
    OMEGA_H_LAMBDA(int j, int& update, bool final) {
      update += nSupports[j];
      if (final) {
        supports_ptr[j + 1] = update;
      }
    },
    total_supports);

  printf("INFO: Inside searchNeighbors 3\n");
  Kokkos::fence();

  supports_idx = Omega_h::Write<Omega_h::LO>(
    total_supports, 0, "index of source supports of each target node");
  printf("INFO: Total_supports: %d\n", total_supports);

  search.adjBasedSearchCentroidNodes(supports_ptr, nSupports, supports_idx,
                                     radii2, false);

  mesh.add_tag<Omega_h::Real>(Omega_h::VERT, "support_radius", 1, radii2);
  return SupportResults{read(supports_ptr), read(supports_idx), radii2};
}

#endif
