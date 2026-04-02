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
 */
#ifndef PCMS_TRANSFER_LOAD_VECTOR_INTEGRATOR_HPP
#define PCMS_TRANSFER_LOAD_VECTOR_INTEGRATOR_HPP

#include <Omega_h_shape.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <vector>
#include <pcms/transfer/mesh_intersection.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms
{
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
    auto ip_vec = MeshField::getIntegrationPoints<MeshField::Triangle>(order);
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
  const IntersectionResults& intersection, const Omega_h::Reals& source_values);
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
                                     const Omega_h::Reals& source_values);
} // namespace pcms

#endif // PCMS_TRANSFER_LOAD_VECTOR_INTEGRATOR_HPP
