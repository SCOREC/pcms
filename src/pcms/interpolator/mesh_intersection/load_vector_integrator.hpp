/**
 * @file load_vector_integrator.hpp
 * @brief Functions for computing load vectors in conservative field projection.
 *
 * This file implements routines to compute the element-wise load vector
 * (right-hand side) contributions used in Galerkin projection of scalar
 * fields from a source mesh to a target mesh.
 *
 * The integration is performed over polygonal intersections (supermesh) between
 * source and target elements using barycentric quadrature. The resulting values
 * represent unassembled local contributions that can later be combined into a
 * global load vector.
 *
 * @note
 * - Assumes 2D linear triangular meshes.
 * - Intersection data is provided via the `IntersectionResults` structure.
 *
 * @created by Abhiyan Paudel
 * @date August 2025
 */
#ifndef PCMS_INTERPOLATOR_LOAD_VECTOR_INTEGRATOR_HPP
#define PCMS_INTERPOLATOR_LOAD_VECTOR_INTEGRATOR_HPP

#include <Omega_h_shape.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <vector>
#include "mesh_intersection.hpp"
#include <Kokkos_MathematicalFunctions.hpp>

/**
 * @brief Computes the load vector for each target element in the conservative
 * field transfer.
 *
 * This routine is used for constructing the right-hand side (RHS) of the
 * conservative field transfer formulation, projecting field quantities from the
 * source mesh to the target mesh.
 *
 * The underlying algorithm computes contributions to the load vector
 * using geometric intersection data between source and target elements.
 *
 * @note Currently this method works for a two-dimensional linear triangles.
 */

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
  const Omega_h::Vector<3>& bary_coords, const int elm_id);
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
 * @brief Provides barycentric integration points and weights for a triangle
 * element.
 *
 * This templated struct stores the barycentric coordinates
 * and quadrature weights for performing numerical integration over a reference
 * triangle. It is used for integrating functions over elements in the
 * conservative field transfer.
 *
 * @tparam order The quadrature order (number of integration points and
 * polynomial accuracy).
 */
template <int order>
struct IntegrationData
{
  // Barycentric coordinates of integration points
  Kokkos::View<MeshField::Vector3*> bary_coords;

  // Quadrature weights associated with each integration point
  Kokkos::View<Omega_h::Real*> weights;

  /**
   * @brief Constructs the integration data for a given quadrature order
   *
   * Initializes barycentric coordinates and weights using
   * MeshField's predefined triangle quadrature rules.
   */
  IntegrationData()
  {
    auto ip_vec = MeshField::getIntegrationPoints(MeshField::Triangle, order);
    std::size_t num_ip = ip_vec.size();

    bary_coords = Kokkos::View<MeshField::Vector3*>("bary_coords", num_ip);
    weights = Kokkos::View<Omega_h::Real*>("weights", num_ip);

    auto bary_coords_host = Kokkos::create_mirror_view(bary_coords);
    auto weights_host = Kokkos::create_mirror_view(weights);

    for (std::size_t i = 0; i < num_ip; ++i) {
      bary_coords_host(i) = ip_vec[i].param;
      weights_host(i) = ip_vec[i].weight;
    }

    Kokkos::deep_copy(bary_coords, bary_coords_host);
    Kokkos::deep_copy(weights, weights_host);
  }

  /**
   * @brief Returns the number of integration points
   *
   * @return Number of integration points  for the selected order.
   */

  int size() const { return bary_coords.extent(0); }
};

/**
 * @brief Computes the counter-clockwise (CCW) vertex ordering of a 2D polygon.
 *
 * Given a 2D polygon stored in an `r3d::Polytope<2>` struct where each vertex
 * has exactly two neighbors. This function reconstructs the CCW traversal order
 * of the polygon's vertices and stores the result in the provided `order`
 * array.
 *
 * The traversal begins at vertex 0 and proceeds by selecting the neighbor that
 * is not previously visited vertex, thereby completing a full loop around the
 * polygon.
 *
 * @param poly The polygon to process, represented as an `r3d::Polytope<2>`.
 * 			   Each vertex includes a list of two neighbor indices
 * (`pnbrs[2]`) forming the cycle.
 *
 * @param[out] order An output array to store the CCW vertex order. Must be
 * preallocated to hold at least `poly.nverts` size.
 * @return The number of vertices in the polygon (i.e., `poly.nverts`)
 *
 * @see r3d::Polytope
 */

[[nodiscard]] OMEGA_H_INLINE int get_polygon_cycle_ccw(
  const r3d::Polytope<2>& poly, int* order)
{
  const int m = poly.nverts;
  if (m <= 0)
    return 0;
  order[0] = 0;
  int prev = -1;
  int curr = 0;
  for (int s = 1; s < m; ++s) {
    const int a = poly.verts[curr].pnbrs[0];
    const int b = poly.verts[curr].pnbrs[1];
    const int next = (a != prev ? a : b); // pick neighbor that's not 'prev'
    order[s] = next;
    prev = curr;
    curr = next;
  }
  return m;
}

/**
 * @brief Computes the per-element RHS load vectors for conservative field
 * projection from source to target mesh.
 *
 * This function computes local (element-wise) right-hand side (RHS)
 * contributions for the Galerkin projection of a scalar field from the source
 * mesh to the target mesh. It integrates over the polygonal intersection
 * regions between each target element and its intersecting source elements
 * using barycentric quadrature.
 *
 * The output is a flat array containing unassembled load vector contributions
 * at the nodes of each target triangle.
 *
 * @param target_mesh The target mesh object receiving the projected scalar
 * field.
 * @param source_mesh The source mesh object containing the original scalar
 * field values.
 * @param intersection Precomputed intersection data for each target element.
 *                     Includes the number and indices of intersecting source
 * elements.
 * @param source_values Scalar field values defined at the nodes of the source
 * mesh.
 *
 * @return A Kokkos view containing per-element load vectors.
 *         Each triangle contributes 3 values (one per node), so the view has
 * size 3 × (number of target elements).
 *
 * @note
 * - This function assumes 2D linear triangular elements.
 * - Degenerate or near-zero-area intersection polygons are skipped.
 * - Each polygon is triangulated using a fan structure and integrated using
 * barycentric quadrature rules.
 * - The returned vector must be assembled into a global RHS vector in a later
 * step.
 *
 * @see evaluate_barycentric, evaluate_function_value, global_from_barycentric
 * @see IntersectionResults
 */

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
  int count = 0;
  Kokkos::parallel_for(
    "calculate load vector", target_mesh.nelems(),
    KOKKOS_LAMBDA(const int& elm) {
      auto tgt_elm_vert_coords =
        get_vert_coords_of_elem(tgt_coords, tgt_faces2nodes, elm);
      const int start = intersection.tgt2src_offsets[elm];
      const int end = intersection.tgt2src_offsets[elm + 1];
      Omega_h::Vector<3> part_integration = {0.0, 0.0, 0.0};

      for (int i = start; i < end; ++i) {
        const int current_src_elm = intersection.tgt2src_indices[i];
        auto src_elm_vert_coords =
          get_vert_coords_of_elem(src_coords, src_faces2nodes, current_src_elm);
        r3d::Polytope<2> poly;
        r3d::intersect_simplices(poly, tgt_elm_vert_coords,
                                 src_elm_vert_coords);
        auto nverts = poly.nverts;
        int order[r3d::MaxVerts<2>::value];
        auto m = get_polygon_cycle_ccw(poly, order);
        auto poly_area = r3d::measure(poly);
        double sum_area = 0;
        for (int j = 1; j < nverts - 1; ++j) {
          // build triangle from poly.verts[order[0]], poly.verts[order[j]],
          // poly.verts[order[j+1]]
          auto& p0 = poly.verts[order[0]].pos;
          auto& p1 = poly.verts[order[j]].pos;
          auto& p2 = poly.verts[order[j + 1]].pos;

          Omega_h::Few<Omega_h::Vector<2>, 3> tri_coords;

          tri_coords[0] = {p0[0], p0[1]};
          tri_coords[1] = {p1[0], p1[1]};
          tri_coords[2] = {p2[0], p2[1]};

          Omega_h::Few<Omega_h::Vector<2>, 2> basis;

          basis[0] = tri_coords[1] - tri_coords[0];
          basis[1] = tri_coords[2] - tri_coords[0];

          Omega_h::Real area =
            Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
          sum_area += area;

          const double EPS_AREA = abs_tol + rel_tol * poly_area;
          if (area <= EPS_AREA)
            continue; // drops duplicates and colinear/degenerates

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
        }
      }

      for (int j = 0; j < 3; ++j) {
        elmLoadVector(elm * 3 + j) = part_integration[j];
      }
    });

  return elmLoadVector;
}

/// Holds projection and conservation error metrics returned by
/// evaluate_pro_and_cons_errors().
struct Errors
{
  double proj_err; ///< L2 projection error computed on the supermesh.
  double cons_err; ///< Relative conservation error over the supermesh.
};

/**
 * @brief Computes projection and conservation errors over the supermesh for
 * scalar field transfer.
 *
 * This function quantifies the accuracy and conservation properties of
 * conservative field transfer between nonconforming meshes using
 * supermesh-based integration. It returns a struct containing two error
 * metrics:
 *
 * - **Projection Error (`proj_err`)** — Measures the L2 norm of the difference
 * between the projected source field and the target field over the
 * supermesh.This reflects how accurately the field has been projected.
 *
 * - **Conservation Error (`cons_err`)** — Relative difference in the integrated
 * field values between source and target representations. This captures
 * conservation loss across the transfer.
 *
 * ### Mathematical Definitions:
 * \f[
 *   \text{proj\_err} = \frac{ \| q_D - q_T \|_{L_2(\Omega_S)} }
 *                           { \| q_D \|_{L_2(\Omega_S)} }, \quad
 *   \text{cons\_err} = \frac{ \left| \int_{\Omega_S} q_D - \int_{\Omega_S} q_T
 * \right| } { \left| \int_{\Omega_S} q_D \right| }
 * \f]
 *
 * where:
 * - \f$q_D\f$ is the scalar field defined on the source mesh (mesh from where
 * the field is defined),
 * - \f$q_T\f$ is the projected field on the target mesh (mesh to where the
 * field is projected),
 * - \f$\Omega_S\f$ is the supermesh formed by polygonal intersections of source
 * and target elements.
 *
 * Integration is performed by triangulating each intersection region and
 * applying barycentric quadrature. Degenerate or near-zero-area triangles are
 * skipped based on area tolerance.
 *
 *
 * @param target_mesh The target mesh object receiving the projected scalar
 * field.
 * @param source_mesh The source mesh object containing the original scalar
 * field values.
 * @param intersection Precomputed intersection data for each target element.
 *                     Includes the number and indices of intersecting source
 * elements.
 * @param target_values Nodal scalar field values evaluated on the target mesh
 * using galerkin projection.
 * @param source_values Scalar field values defined at the nodes of the source
 * mesh.
 *
 *
 * @return A struct containing:
 *   - `proj_err`: Exact L2 projection error over the supermesh.
 *   - `cons_err`: Relative conservation error over the supermesh.
 *
 * @note
 * - Assumes 2D linear (P1) triangular elements.
 * - Ideal for validating conservative transfer schemes or testing projection
 * fidelity.
 *
 * @see IntersectionResults, buildLoadVector
 */

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
      auto tgt_elm_vert_coords =
        get_vert_coords_of_elem(tgt_coords, tgt_faces2nodes, elm);
      const int start = intersection.tgt2src_offsets[elm];
      const int end = intersection.tgt2src_offsets[elm + 1];

      double N2 = 0.0, D2 = 0.0, C = 0.0, QD = 0.0;

      for (int i = start; i < end; ++i) {
        const int current_src_elm = intersection.tgt2src_indices[i];
        auto src_elm_vert_coords =
          get_vert_coords_of_elem(src_coords, src_faces2nodes, current_src_elm);
        r3d::Polytope<2> poly;
        r3d::intersect_simplices(poly, tgt_elm_vert_coords,
                                 src_elm_vert_coords);
        auto nverts = poly.nverts;
        int order[r3d::MaxVerts<2>::value];
        auto m = get_polygon_cycle_ccw(poly, order);
        auto poly_area = r3d::measure(poly);

        double sum_area = 0.0;
        for (int j = 1; j < nverts - 1; ++j) {
          // build triangle from poly.verts[order[0]], poly.verts[order[j]],
          // poly.verts[order[j+1]]
          auto& p0 = poly.verts[order[0]].pos;
          auto& p1 = poly.verts[order[j]].pos;
          auto& p2 = poly.verts[order[j + 1]].pos;

          Omega_h::Few<Omega_h::Vector<2>, 3> tri_coords;

          tri_coords[0] = {p0[0], p0[1]};
          tri_coords[1] = {p1[0], p1[1]};
          tri_coords[2] = {p2[0], p2[1]};

          Omega_h::Few<Omega_h::Vector<2>, 2> basis;

          basis[0] = tri_coords[1] - tri_coords[0];
          basis[1] = tri_coords[2] - tri_coords[0];

          Omega_h::Real area =
            Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
          sum_area += area;

          const double EPS_AREA = abs_tol + rel_tol * ploy_area;
          if (area <= EPS_AREA)
            continue; // drops duplicates and colinear/degenerates

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
        }
      }

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

#endif
