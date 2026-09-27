#ifndef PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP
#define PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP

#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/utility/assert.h"
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <MeshField_Config.hpp>
#include <Omega_h_shape.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms::detail
{

// Shared checks for a scalar Cartesian Lagrange space on a simplex mesh
// (triangles in 2D, tetrahedra in 3D), independent of order. Order is validated
// separately by the callers below.
inline void CheckOmegaHScalarSimplexLayout(
  CoordinateSystem coordinate_system,
  const std::shared_ptr<const OmegaHLagrangeLayout>& layout,
  const char* context, const char* role)
{
  if (layout == nullptr) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must use OmegaHLagrangeLayout");
  }
  if (layout->GetNumComponents() != 1) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must have exactly one component");
  }
  if (coordinate_system != CoordinateSystem::Cartesian) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must use Cartesian coordinates");
  }
  const Omega_h::Mesh& mesh = layout->GetMesh();
  if (mesh.dim() != 2 && mesh.dim() != 3) {
    throw pcms_error(std::string(context) + ": " + role +
                     " mesh must be 2D or 3D");
  }
  if (mesh.family() != OMEGA_H_SIMPLEX) {
    throw pcms_error(std::string(context) + ": " + role +
                     " mesh must be a simplex (triangle/tetrahedron) mesh");
  }
}

// Strict order-1 check (used where only P1 is supported, e.g. Monte-Carlo RHS).
inline void CheckOmegaHScalarP1Layout(
  CoordinateSystem coordinate_system,
  const std::shared_ptr<const OmegaHLagrangeLayout>& layout,
  const char* context, const char* role)
{
  CheckOmegaHScalarSimplexLayout(coordinate_system, layout, context, role);
  if (layout->GetOrder() != 1) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must be order-1");
  }
}

// Conservative-projection check: any supported Lagrange order (P0 or P1). The
// intersection integrator handles source and target orders independently, so
// this replaces the strict P1 requirement on those paths.
inline void CheckOmegaHScalarLagrangeLayout(
  CoordinateSystem coordinate_system,
  const std::shared_ptr<const OmegaHLagrangeLayout>& layout,
  const char* context, const char* role)
{
  CheckOmegaHScalarSimplexLayout(coordinate_system, layout, context, role);
  const int order = layout->GetOrder();
  if (order != 0 && order != 1) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must be order-0 or order-1");
  }
}

// Map barycentric coordinates on a simplex (Dim+1 barycentric components) to
// the global Cartesian point, given the simplex's Dim+1 vertex coordinates.
template <int Dim>
[[nodiscard]] OMEGA_H_INLINE Omega_h::Vector<Dim> GlobalFromBarycentric(
  const Omega_h::Vector<Dim + 1>& barycentric_coord,
  const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& verts_coord)
{
  Omega_h::Vector<Dim> real_coords;
  for (int d = 0; d < Dim; ++d) {
    real_coords[d] = 0.0;
  }
  for (int i = 0; i < Dim + 1; ++i) {
    for (int d = 0; d < Dim; ++d) {
      real_coords[d] += barycentric_coord[i] * verts_coord[i][d];
    }
  }
  return real_coords;
}

// Walk each face of a clipped r3d polyhedron exactly once and fan it into
// triangles, invoking op(v0, v1, v2) for each triangle (v0 is the face's anchor
// vertex, so a face with k vertices yields k-2 triangles). Every vertex of an
// r3d clipped Polytope<3> has exactly three face-neighbors (pnbrs); marking
// each directed edge as it is consumed guarantees every face is emitted once.
// This mirrors the edge-marking traversal buried inside r3d::reduce, which r3d
// does not expose for reuse, so the traversal is reproduced here once and
// shared. Vertices arrive as r3d::Vector<3> (indexable [0..2]).
template <typename TriangleOp>
OMEGA_H_INLINE void ForEachPolytopeFaceTriangle(const r3d::Polytope<3>& poly,
                                                TriangleOp&& op)
{
  // emarks[v][p] == 1 once the directed edge (v, pnbr p) has been consumed.
  int emarks[r3d::Polytope<3>::max_verts][3] = {{}};
  for (int vstart = 0; vstart < poly.nverts; ++vstart) {
    for (int pstart = 0; pstart < 3; ++pstart) {
      if (emarks[vstart][pstart]) {
        continue;
      }
      int pnext = pstart;
      int vcur = vstart;
      emarks[vcur][pnext] = 1;
      int vnext = poly.verts[vcur].pnbrs[pnext];
      const auto face_v0 = poly.verts[vcur].pos;

      // Move to the second edge of this face.
      int np = 0;
      for (np = 0; np < 3; ++np) {
        if (poly.verts[vnext].pnbrs[np] == vcur) {
          break;
        }
      }
      vcur = vnext;
      pnext = (np + 1) % 3;
      emarks[vcur][pnext] = 1;
      vnext = poly.verts[vcur].pnbrs[pnext];

      // Fan the face into triangles anchored at face_v0.
      while (vnext != vstart) {
        op(face_v0, poly.verts[vnext].pos, poly.verts[vcur].pos);

        // Advance around the face.
        for (np = 0; np < 3; ++np) {
          if (poly.verts[vnext].pnbrs[np] == vcur) {
            break;
          }
        }
        vcur = vnext;
        pnext = (np + 1) % 3;
        emarks[vcur][pnext] = 1;
        vnext = poly.verts[vcur].pnbrs[pnext];
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Star decomposition of a clipped r3d polytope into quadrature simplices.
// ---------------------------------------------------------------------------

// Boundary (Dim-1)-simplices of a clipped r3d polytope as its vertex graph
// encodes them; op(facet) receives the facet's Dim vertices in walk order.
//
// This is the one dimension-specific step of the decomposition, because r3d's
// graph has no dimension-generic boundary iterator (r3d::reduce is itself
// written once per dimension). In 2D, pnbrs[0] is the next vertex around the
// single cycle and every edge is a facet. In 3D, the three pnbrs are the
// cyclically ordered neighbors and faces are recovered by the edge-marking walk
// in ForEachPolytopeFaceTriangle, which also fans them into triangles.
template <typename FacetOp>
OMEGA_H_INLINE void ForEachPolytopeBoundarySimplex(const r3d::Polytope<2>& poly,
                                                   FacetOp&& op)
{
  for (int v = 0; v < poly.nverts; ++v) {
    r3d::Few<r3d::Vector<2>, 2> facet;
    facet[0] = poly.verts[v].pos;
    facet[1] = poly.verts[poly.verts[v].pnbrs[0]].pos;
    op(facet);
  }
}

template <typename FacetOp>
OMEGA_H_INLINE void ForEachPolytopeBoundarySimplex(const r3d::Polytope<3>& poly,
                                                   FacetOp&& op)
{
  ForEachPolytopeFaceTriangle(poly, [&](const r3d::Vector<3>& a,
                                        const r3d::Vector<3>& b,
                                        const r3d::Vector<3>& c) {
    r3d::Few<r3d::Vector<3>, 3> facet;
    facet[0] = a;
    facet[1] = b;
    facet[2] = c;
    op(facet);
  });
}

// Signed measure of the simplex (apex, facet[0], ..., facet[Dim-1]) in the
// order the boundary walk emits it.
template <int Dim>
[[nodiscard]] OMEGA_H_INLINE Omega_h::Real StarSimplexSignedMeasure(
  const Omega_h::Vector<Dim>& apex,
  const r3d::Few<r3d::Vector<Dim>, Dim>& facet)
{
  Omega_h::Few<Omega_h::Vector<Dim>, Dim> basis;
  for (int i = 0; i < Dim; ++i) {
    for (int d = 0; d < Dim; ++d) {
      basis[i][d] = facet[i][d] - apex[d];
    }
  }
  return Omega_h::simplex_size_from_basis(basis);
}

// Star-decompose a clipped r3d polytope into simplices, invoking
// op(simplex, measure) for every piece with |measure| > eps. The pieces'
// measures sum to the polytope's measure.
//
// The apex is the polytope's first vertex: every boundary facet (from
// ForEachPolytopeBoundarySimplex) is lifted to it, and the facets incident to
// that vertex give zero-measure pieces that the epsilon filter drops. Any apex
// works, because by the divergence theorem the signed sum over a closed
// boundary is the enclosed measure regardless of where the pieces meet; a
// vertex apex just yields the fewest pieces (one for a clean simplex, against
// Dim+1 from the centroid), and each piece is one quadrature rule downstream.
//
// The measures are *signed*, and that matters: when the clipping planes pass
// through the polytope's own vertices (a source element sharing a face plane
// with the target, always the case on the same mesh), r3d sees signed distances
// of order 1e-13 with mixed signs and splices new vertices at O(1) fractions
// along the edges between them. The measure it reports is still exact because
// r3d integrates with signed pieces, but the boundary graph is folded: some
// facets come back with reversed orientation and cancel against their mirror.
// Taking |measure| per piece breaks that cancellation and over-counts by up to
// tens of percent of the element, so the sign is kept and the overall
// orientation is fixed from the total.
//
// The walk's handedness is not fixed by r3d (it follows the order the clip left
// the neighbor lists in), so the sign that makes the pieces add up to a
// positive measure is read off the total. Folded pairs cancel in that sum, so
// it is the true measure to roundoff and its sign is well defined.
template <int Dim, typename SimplexOp>
OMEGA_H_INLINE void ForEachPolytopeStarSimplex(const r3d::Polytope<Dim>& poly,
                                               const double eps, SimplexOp&& op)
{
  Omega_h::Vector<Dim> apex;
  for (int d = 0; d < Dim; ++d) {
    apex[d] = poly.verts[0].pos[d];
  }

  Omega_h::Real signed_total = 0.0;
  ForEachPolytopeBoundarySimplex(
    poly, [&](const r3d::Few<r3d::Vector<Dim>, Dim>& facet) {
      signed_total += StarSimplexSignedMeasure<Dim>(apex, facet);
    });
  const Omega_h::Real orientation = (signed_total < 0.0) ? -1.0 : 1.0;

  ForEachPolytopeBoundarySimplex(
    poly, [&](const r3d::Few<r3d::Vector<Dim>, Dim>& facet) {
      const Omega_h::Real measure =
        orientation * StarSimplexSignedMeasure<Dim>(apex, facet);
      if (Kokkos::fabs(measure) > eps) {
        Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1> simplex;
        simplex[0] = apex;
        for (int i = 0; i < Dim; ++i) {
          for (int d = 0; d < Dim; ++d) {
            simplex[i + 1][d] = facet[i][d];
          }
        }
        op(simplex, measure);
      }
    });
}

// Driver over the sub-simplices (triangles in 2D, tets in 3D) that tile each
// target element's intersection with the source mesh: clip the target against
// every intersecting source element and star-decompose the result. Invokes
// op(sub_simplex_coords, src_elm, measure) for every non-degenerate piece. The
// measure is signed (see ForEachPolytopeStarSimplex); the pieces of one
// intersection always sum to its measure, so callers must use the measure as a
// weight rather than a size.
template <int Dim, typename SimplexOp>
OMEGA_H_INLINE void ForEachIntersectionSubsimplex(
  const int elm, const IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_elems2nodes, const Omega_h::LOs& src_elems2nodes,
  SimplexOp&& op)
{
  auto tgt_elm_vert_coords =
    get_vert_coords_of_elem<Dim>(tgt_coords, tgt_elems2nodes, elm);
  const int start = intersection.tgt2src_offsets[elm];
  const int end = intersection.tgt2src_offsets[elm + 1];

  for (int i = start; i < end; ++i) {
    const int current_src_elm = intersection.tgt2src_indices[i];
    auto src_elm_vert_coords = get_vert_coords_of_elem<Dim>(
      src_coords, src_elems2nodes, current_src_elm);
    r3d::Polytope<Dim> poly;
    r3d::intersect_simplices(poly, tgt_elm_vert_coords, src_elm_vert_coords);
    if (poly.nverts < Dim + 1) {
      continue;
    }
    const double poly_measure = Kokkos::fabs(r3d::measure(poly));
    const double eps =
      PCMS_INTERSECTION_ABS_TOL + PCMS_INTERSECTION_REL_TOL * poly_measure;

    ForEachPolytopeStarSimplex<Dim>(
      poly, eps,
      [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& simplex,
          Omega_h::Real measure) { op(simplex, current_src_elm, measure); });
  }
}

// Maps spatial dimension to the MeshField simplex topology whose reference
// quadrature rules we use (triangle in 2D, tetrahedron in 3D).
template <int Dim>
struct SimplexTopology;
template <>
struct SimplexTopology<2>
{
  static constexpr MeshField::Mesh_Topology value = MeshField::Triangle;
};
template <>
struct SimplexTopology<3>
{
  static constexpr MeshField::Mesh_Topology value = MeshField::Tetrahedron;
};

// Barycentric integration points and weights for a reference simplex (triangle
// in 2D, tetrahedron in 3D), taken from MeshField's predefined quadrature rules
// and staged on device for use in element integration kernels.
//
// MeshField::getIntegrationPoints returns a host std::vector, which cannot be
// dereferenced inside a device kernel, so the (tiny) rule is copied into device
// Kokkos views once at construction. Each barycentric point has Dim+1
// components.
//
// The quadrature order is a runtime argument because the required polynomial
// accuracy depends on the source and target element orders (degree =
// source_order + target_order), which are only known at construction.
template <int Dim>
struct IntegrationData
{
  Kokkos::View<Omega_h::Real * [Dim + 1]>
    bary_coords;                        // barycentric coordinates
  Kokkos::View<Omega_h::Real*> weights; // quadrature weights

  explicit IntegrationData(int order)
  {
    auto ip_vec =
      MeshField::getIntegrationPoints<SimplexTopology<Dim>::value>(order);
    const std::size_t num_ip = ip_vec.size();

    bary_coords =
      Kokkos::View<Omega_h::Real * [Dim + 1]>("bary_coords", num_ip);
    weights = Kokkos::View<Omega_h::Real*>("weights", num_ip);

    auto bary_coords_host = Kokkos::create_mirror_view(bary_coords);
    auto weights_host = Kokkos::create_mirror_view(weights);
    for (std::size_t i = 0; i < num_ip; ++i) {
      // MeshField returns points in reduced parametric coordinates: only the
      // first Dim barycentric components are stored, with the last implied by
      // the partition of unity. Expand to the full Dim+1 barycentric form
      // consumed by GlobalFromBarycentric.
      Omega_h::Real last = 1.0;
      for (int d = 0; d < Dim; ++d) {
        const Omega_h::Real xi = ip_vec[i].param[d];
        bary_coords_host(i, d) = xi;
        last -= xi;
      }
      bary_coords_host(i, Dim) = last;
      weights_host(i) = ip_vec[i].weight;
    }
    Kokkos::deep_copy(bary_coords, bary_coords_host);
    Kokkos::deep_copy(weights, weights_host);
  }

  int size() const { return bary_coords.extent(0); }
};

// Target Lagrange basis on a simplex (triangle in 2D, tetrahedron in 3D),
// parameterized by spatial dimension and element order, for the
// conservative-projection RHS assembly. Order 0 is a single element-constant
// DOF; order 1 is the Dim+1 vertex (barycentric) DOFs. Higher orders slot in as
// additional specializations, kept in lock-step with element_dispatch.h.
//
// Each specialization provides, for a target element `elm` with local vertex
// ids `verts` and vertex coordinates `tgt_verts`:
//   ndof            number of local target DOFs
//   Index(...)      active PETSc row for local dof k
//   Values(pt, ...) basis values at the (global) integration point pt
template <int Dim, int Order>
struct TargetSimplexBasis;

template <int Dim>
struct TargetSimplexBasis<Dim, 0>
{
  static constexpr int ndof = 1;

  template <typename Permutation>
  KOKKOS_INLINE_FUNCTION static LO Index(
    const Permutation& permutation, int elm,
    const Omega_h::Few<Omega_h::LO, Dim + 1>&, int /*k*/)
  {
    return permutation(elm);
  }

  KOKKOS_INLINE_FUNCTION static void Values(
    const Omega_h::Vector<Dim>&,
    const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&, Omega_h::Real out[ndof])
  {
    out[0] = 1.0;
  }
};

template <int Dim>
struct TargetSimplexBasis<Dim, 1>
{
  static constexpr int ndof = Dim + 1;

  template <typename Permutation>
  KOKKOS_INLINE_FUNCTION static LO Index(
    const Permutation& permutation, int /*elm*/,
    const Omega_h::Few<Omega_h::LO, Dim + 1>& verts, int k)
  {
    return permutation(verts[k]);
  }

  // P1 basis functions are the barycentric coordinates of the target element
  // evaluated at the (global) integration point.
  KOKKOS_INLINE_FUNCTION static void Values(
    const Omega_h::Vector<Dim>& pt,
    const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& tgt_verts,
    Omega_h::Real out[ndof])
  {
    const auto bary = Omega_h::barycentric_from_global<Dim, Dim>(pt, tgt_verts);
    for (int i = 0; i < ndof; ++i) {
      out[i] = bary[i];
    }
  }
};

} // namespace pcms::detail

#endif // PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP
