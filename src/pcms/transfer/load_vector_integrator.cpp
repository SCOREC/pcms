#include "pcms/transfer/load_vector_integrator.hpp"

namespace pcms
{

/**
 * @brief Converts barycentric coordinates to global (physical) coordinates.
 *
 * Given barycentric coordinates within a 2D triangle and the coordinates of
 * the triangle's vertices, this function computes the corresponding global
 * position.
 *
 * @param barycentric_coord The barycentric coordinates \f$(\lambda_1,
 * \lambda_2, \lambda_3)\f$ of the point.
 * @param verts_coord The coordinates of the triangle's three vertices in global
 * space.
 * @return The 2D global coordinates corresponding to the given barycentric
 * position.
 */

[[nodiscard]] OMEGA_H_INLINE Omega_h::Vector<2> global_from_barycentric(
  const MeshField::Vector3& barycentric_coord,
  const Omega_h::Few<Omega_h::Vector<2>, 3>& verts_coord)
{
  Omega_h::Vector<2> real_coords = {0.0, 0.0};

  for (int i = 0; i < 3; ++i) {
    real_coords[0] += barycentric_coord[i] * verts_coord[i][0];
    real_coords[1] += barycentric_coord[i] * verts_coord[i][1];
  }
  return real_coords;
}

/**
 * @brief Computes the barycentric coordinates of a 2D point with respect to a
 * triangle.
 *
 * Given a point in global (x, y) coordinates and the coordinates of the three
 * vertices of a triangle, this function evaluates the barycentric coordinates
 * \f$(\lambda_1, \lambda_2, \lambda_3)\f$ of the point with respect to that
 * triangle.
 *
 * @param point The 2D global coordinates of the point to evaluate (in
 * Omega_h::Vector<2> format).
 * @param verts_coord The vertex coordinates of the triangle (in r3d::Vector<2>
 * format).
 * @return A vector of three barycentric coordinates corresponding to the input
 * point.
 *
 */

[[nodiscard]] OMEGA_H_INLINE Omega_h::Vector<3> evaluate_barycentric(
  const Omega_h::Vector<2>& point,
  const r3d::Few<r3d::Vector<2>, 3>& verts_coord)
{
  Omega_h::Few<Omega_h::Vector<2>, 3> omegah_vector;
  for (int i = 0; i < 3; ++i) {
    omegah_vector[i][0] = verts_coord[i][0];
    omegah_vector[i][1] = verts_coord[i][1];
  }

  auto barycentric_coordinate =
    Omega_h::barycentric_from_global<2, 2>(point, omegah_vector);

  return barycentric_coordinate;
}
/**
 * @brief Evaluates the value of a linear function at a given point using
 * barycentric coordinates.
 *
 * This function computes the interpolated value of a nodal scalar field over a
 * triangle, using barycentric coordinates within the specified element.
 *
 * @param nodal_values The global array of nodal field values.
 * @param faces2nodes The element-to-node connectivity array.
 * @param bary_coords The barycentric coordinates of the evaluation point within
 * the triangle.
 * @param elm_id The ID of the triangle element being evaluated.
 * @return The interpolated function value at the given point.
 *
 * @note This function assumes linear (3-node) triangular element.
 */

[[nodiscard]] OMEGA_H_INLINE double evaluate_function_value(
  const Omega_h::Reals& nodal_values, const Omega_h::LOs& faces2nodes,
  const Omega_h::Vector<3>& bary_coords, const int elm_id)
{

  double value = 0;
  const auto elm_verts = Omega_h::gather_verts<3>(faces2nodes, elm_id);
  for (int i = 0; i < 3; ++i) {
    int nid = elm_verts[i];
    value += nodal_values[nid] * bary_coords[i];
  }

  return value;
}

/**
 * @brief Deduplicate, reorder, and orient a 2D polygon produced by r3d.
 *
 * This function performs all necessary cleanup and reconstruction of a polygon
 * returned by `r3d::intersect_simplices()`, which may contain:
 *  - duplicated vertices,
 *  - unordered vertex lists,
 *  - invalid or inconsistent neighbor links (`pnbrs`),
 *  - negative orientation (CW instead of CCW).
 *
 * The cleanup proceeds with the following stages:
 *
 * **(1) Geometric deduplication:**
 * Vertices whose coordinates are equal within a tolerance `tol` are collapsed
 * into a single unique vertex. A compacted vertex list is built.
 *
 * **(2) CCW vertex reordering:**
 * The remaining unique vertices are sorted by their polar angle around the
 * polygon centroid. This yields a globally consistent counter-clockwise (CCW)
 *
 *
 * **(3) Rebuilding neighbor links:**
 * After sorting, each vertex's two neighbors (`pnbrs[0]` and `pnbrs[1]`) are
 * reassigned to form a closed CCW cycle:
 *
 *      pnbrs[0] = previous vertex in CCW order
 *      pnbrs[1] = next     vertex in CCW order
 *
 *
 *
 * @param poly The polygon to clean, and reorder.
 *
 * @param tol Tolerance for geometric duplicate detection (default: 1e-12).
 *
 * @return The number of unique, CCW-ordered vertices remaining in the polygon.
 *
 * @see r3d::Polytope
 * @see r3d::measure
 * @see r3d::intersect_simplices
 */

[[nodiscard]] OMEGA_H_INLINE int remove_duplicate_vertices_and_fix_links(
  r3d::Polytope<2>& poly, const double tol = 1e-12)
{

  const int old_n = poly.nverts;
  int new_n = 0;

  int old2new[r3d::MaxVerts<2>::value];
  for (int i = 0; i < old_n; ++i) {
    old2new[i] = -1;
  }

  // geometric deupliactes filtering
  for (int i = 0; i < old_n; ++i) {
    const auto& pi = poly.verts[i].pos;
    bool dup = false;
    for (int j = 0; j < new_n; ++j) {
      const auto& pj = poly.verts[j].pos;
      if (Kokkos::fabs(pi[0] - pj[0]) < tol &&
          Kokkos::fabs(pi[1] - pj[1]) < tol) {
        dup = true;
        old2new[i] = j;
        break;
      }
    }
    if (!dup) {
      old2new[i] = new_n;
      poly.verts[new_n] = poly.verts[i];
      ++new_n;
    }
  }
  poly.nverts = new_n;

  //  CCW reorder verts by angle about centroid
  if (new_n >= 3) {
    // centroid
    double cx = 0.0;
    double cy = 0.0;
    for (int i = 0; i < new_n; ++i) {
      cx += poly.verts[i].pos[0];
      cy += poly.verts[i].pos[1];
    }
    cx /= new_n;
    cy /= new_n;

    // angle sort
    int order[r3d::MaxVerts<2>::value];
    for (int i = 0; i < new_n; ++i) {
      order[i] = i;
    }

    for (int i = 0; i < new_n - 1; ++i) {
      for (int j = i + 1; j < new_n; ++j) {
        const double a1 = Kokkos::atan2(poly.verts[order[i]].pos[1] - cy,
                                        poly.verts[order[i]].pos[0] - cx);
        const double a2 = Kokkos::atan2(poly.verts[order[j]].pos[1] - cy,
                                        poly.verts[order[j]].pos[0] - cx);
        if (a1 > a2) {
          int t = order[i];
          order[i] = order[j];
          order[j] = t;
        }
      }
    }

    r3d::Vertex<2> tmp[r3d::MaxVerts<2>::value];
    for (int i = 0; i < new_n; ++i) {
      tmp[i] = poly.verts[order[i]];
    }

    for (int i = 0; i < new_n; ++i) {
      poly.verts[i] = tmp[i];
    }

    //  set circular neighbors consistent with verts[] order
    for (int i = 0; i < new_n; ++i) {
      const int prev = (i - 1 + new_n) % new_n;
      const int next = (i + 1) % new_n;
      poly.verts[i].pnbrs[0] = prev; // CCW prev
      poly.verts[i].pnbrs[1] = next; // CCW next
    }

    //  ensure positive orientation (if measure is still negative, reverse)
    double area = r3d::measure(poly);
    if (area < 0.0) {
      // reverse verts and relink
      for (int i = 0; i < new_n / 2; ++i) {
        r3d::Vertex<2> t = poly.verts[i];
        poly.verts[i] = poly.verts[new_n - 1 - i];
        poly.verts[new_n - 1 - i] = t;
      }
      for (int i = 0; i < new_n; ++i) {
        const int prev = (i - 1 + new_n) % new_n;
        const int next = (i + 1) % new_n;
        poly.verts[i].pnbrs[0] = prev;
        poly.verts[i].pnbrs[1] = next;
      }
    }
  } else {
    // clear links
    for (int i = 0; i < new_n; ++i) {
      poly.verts[i].pnbrs[0] = -1;
      poly.verts[i].pnbrs[1] = -1;
    }
  }

  return new_n;
}

template <typename TriangleOp>
OMEGA_H_INLINE void for_each_intersection_subtriangle(
  const int elm, const IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_faces2nodes, const Omega_h::LOs& src_faces2nodes,
  TriangleOp&& op)
{
  auto tgt_elm_vert_coords =
    get_vert_coords_of_elem(tgt_coords, tgt_faces2nodes, elm);
  const int start = intersection.tgt2src_offsets[elm];
  const int end = intersection.tgt2src_offsets[elm + 1];

  for (int i = start; i < end; ++i) {
    const int current_src_elm = intersection.tgt2src_indices[i];
    auto src_elm_vert_coords =
      get_vert_coords_of_elem(src_coords, src_faces2nodes, current_src_elm);
    r3d::Polytope<2> poly;
    r3d::intersect_simplices(poly, tgt_elm_vert_coords, src_elm_vert_coords);
    auto nverts = remove_duplicate_vertices_and_fix_links(poly, 1e-12);
    ;
    auto poly_area = r3d::measure(poly);

    for (int j = 1; j < nverts - 1; ++j) {
      // build triangle from poly.verts[0], poly.verts[j],
      // poly.verts[j+1]
      auto& p0 = poly.verts[0].pos;
      auto& p1 = poly.verts[j].pos;
      auto& p2 = poly.verts[j + 1].pos;

      Omega_h::Few<Omega_h::Vector<2>, 3> tri_coords;
      tri_coords[0] = {p0[0], p0[1]};
      tri_coords[1] = {p1[0], p1[1]};
      tri_coords[2] = {p2[0], p2[1]};

      Omega_h::Few<Omega_h::Vector<2>, 2> basis;
      basis[0] = tri_coords[1] - tri_coords[0];
      basis[1] = tri_coords[2] - tri_coords[0];

      Omega_h::Real area =
        Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));

      const double EPS_AREA = abs_tol + rel_tol * poly_area;
      if (area <= EPS_AREA)
        continue; // drops duplicates and colinear/degenerates

      op(tri_coords, tgt_elm_vert_coords, src_elm_vert_coords, current_src_elm,
         area);
    }
  }
}

Kokkos::View<MeshField::Real*> buildLoadVector(
  Omega_h::Mesh& target_mesh, Omega_h::Mesh& source_mesh,
  const IntersectionResults& intersection, const Omega_h::Reals& source_values)
{

  const auto& tgt_coords = target_mesh.coords();
  const auto& src_coords = source_mesh.coords();
  const auto& tgt_faces2nodes =
    target_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& src_faces2nodes =
    source_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  IntegrationData<2> integrationPoints;
  int npts = integrationPoints.size();

  // TODO: Make it generalised; hardcoded for liner 2D
  Kokkos::View<MeshField::Real*> elmLoadVector("elmLoadVector",
                                               target_mesh.nelems() * 3);
  Kokkos::parallel_for(
    "calculate load vector", target_mesh.nelems(),
    KOKKOS_LAMBDA(const int& elm) {
      Omega_h::Vector<3> part_integration = {0.0, 0.0, 0.0};
      for_each_intersection_subtriangle(
        elm, intersection, tgt_coords, src_coords, tgt_faces2nodes,
        src_faces2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<2>, 3>& tri_coords,
            const r3d::Few<r3d::Vector<2>, 3>& tgt_elm_vert_coords,
            const r3d::Few<r3d::Vector<2>, 3>& src_elm_vert_coords,
            const int current_src_elm, const Omega_h::Real area) {
          for (int ip = 0; ip < npts; ++ip) {
            auto bary = integrationPoints.bary_coords(ip);
            auto weight = integrationPoints.weights(ip);

            // convert barycentric to real coords in triangle
            auto real_coords = global_from_barycentric(bary, tri_coords);

            // evaluate shape function (barycentric wrt target for linear)
            auto shape_fn =
              evaluate_barycentric(real_coords, tgt_elm_vert_coords);

            // evaluate function at point (barycentric wrt source for linear)
            auto src_bary =
              evaluate_barycentric(real_coords, src_elm_vert_coords);
            auto fval = evaluate_function_value(source_values, src_faces2nodes,
                                                src_bary, current_src_elm);

            // integration
            for (int k = 0; k < 3; ++k) {
              part_integration[k] += shape_fn[k] * fval * weight * 2 * area;
            }
          }
        });

      for (int j = 0; j < 3; ++j) {
        elmLoadVector(elm * 3 + j) = part_integration[j];
      }
    });

  return elmLoadVector;
}
Errors evaluate_proj_and_cons_errors(Omega_h::Mesh& target_mesh,
                                     Omega_h::Mesh& source_mesh,
                                     const IntersectionResults& intersection,
                                     const Omega_h::Reals& target_values,
                                     const Omega_h::Reals& source_values)
{

  const auto& tgt_coords = target_mesh.coords();
  const auto& src_coords = source_mesh.coords();
  const auto& tgt_faces2nodes =
    target_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& src_faces2nodes =
    source_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  IntegrationData<2> integrationPoints;
  int npts = integrationPoints.size();

  constexpr double EPS_DEN = 1e-30;

  Kokkos::View<double*> accum("accum", 4);
  Kokkos::deep_copy(accum, 0.0);

  Kokkos::parallel_for(
    "evaluate relative errors", target_mesh.nelems(),
    KOKKOS_LAMBDA(const int& elm) {
      double N2 = 0.0, D2 = 0.0, C = 0.0, QD = 0.0;
      for_each_intersection_subtriangle(
        elm, intersection, tgt_coords, src_coords, tgt_faces2nodes,
        src_faces2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<2>, 3>& tri_coords,
            const r3d::Few<r3d::Vector<2>, 3>& tgt_elm_vert_coords,
            const r3d::Few<r3d::Vector<2>, 3>& src_elm_vert_coords,
            const int current_src_elm, const Omega_h::Real area) {
          for (int ip = 0; ip < npts; ++ip) {
            auto bary = integrationPoints.bary_coords(ip);
            auto weight = integrationPoints.weights(ip);

            // convert barycentric to real coords in triangle
            auto real_coords = global_from_barycentric(bary, tri_coords);
            auto tgt_bary =
              evaluate_barycentric(real_coords, tgt_elm_vert_coords);

            // evaluate shape function (barycentric wrt target for linear)
            auto tgtVal = evaluate_function_value(
              target_values, tgt_faces2nodes, tgt_bary, elm);

            // evaluate function at point (barycentric wrt source for linear)
            auto src_bary =
              evaluate_barycentric(real_coords, src_elm_vert_coords);
            auto srcVal = evaluate_function_value(
              source_values, src_faces2nodes, src_bary, current_src_elm);

            // integration
            auto diff = srcVal - tgtVal;
            auto w = 2 * weight * area;
            N2 += diff * diff * w;
            D2 += srcVal * srcVal * w;
            C += diff * w;
            QD += srcVal * w;
          }
        });

      Kokkos::atomic_add(&accum(0), N2);
      Kokkos::atomic_add(&accum(1), D2);
      Kokkos::atomic_add(&accum(2), C);
      Kokkos::atomic_add(&accum(3), QD);
    });

  auto h_accum = Kokkos::create_mirror(accum);
  Kokkos::deep_copy(h_accum, accum);
  const double proj_err =
    Kokkos::sqrt(h_accum(0)) / Kokkos::max(Kokkos::sqrt(h_accum(1)), EPS_DEN);
  const double cons_err =
    Kokkos::fabs(h_accum(2)) / Kokkos::max(Kokkos::fabs(h_accum(3)), EPS_DEN);

  return Errors{.proj_err = proj_err, .cons_err = cons_err};
}
} // namespace pcms
