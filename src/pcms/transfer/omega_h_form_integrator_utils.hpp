#ifndef PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP
#define PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP

#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/backends/omega_h/mesh_intersection.hpp"
#include "pcms/backends/omega_h/simplex_quadrature.hpp"
#include "pcms/utility/assert.hpp"
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <MeshField_Config.hpp>
#include <Omega_h_shape.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms::transfer::detail
{

using omega_h::detail::IntegrationData;

// The layout and Omega_h geometry of a scalar Lagrange space on a simplex mesh
// (triangles in 2D, tetrahedra in 3D) in Cartesian coordinates.
struct OmegaHScalarSpace
{
  std::shared_ptr<const omega_h::OmegaHLagrangeLayout> layout;
  std::shared_ptr<const omega_h::OmegaHGeometry> geometry;
};

// Checks that `space` is a scalar Lagrange space of an order in [min_order,
// max_order] on the Omega_h geometry of a 2D or 3D simplex mesh, in the
// Cartesian coordinate system of the mesh dimension, and returns its layout and
// geometry. The integrators take their integrals in the mesh coordinates
// without a metric, which is only the physical measure in Cartesian
// coordinates; other systems are refused until the geometry's metric enters
// the measures.
inline OmegaHScalarSpace ResolveOmegaHScalarSpace(const FunctionSpace& space,
                                                  int min_order, int max_order,
                                                  const char* context,
                                                  const char* role)
{
  const std::string prefix = std::string(context) + ": " + role + " ";
  auto layout = std::dynamic_pointer_cast<const omega_h::OmegaHLagrangeLayout>(
    space.GetLayout());
  if (layout == nullptr) {
    throw pcms_error(prefix + "space must use OmegaHLagrangeLayout");
  }
  if (layout->GetNumComponents() != 1) {
    throw pcms_error(prefix + "space must have exactly one component");
  }
  auto geometry = std::dynamic_pointer_cast<const omega_h::OmegaHGeometry>(
    space.GetGeometry());
  if (geometry == nullptr || &geometry->GetMesh() != &layout->GetMesh()) {
    throw pcms_error(prefix +
                     "space must be on the Omega_h geometry of its layout's "
                     "mesh");
  }
  const Omega_h::Mesh& mesh = geometry->GetMesh();
  if (mesh.dim() != 2 && mesh.dim() != 3) {
    throw pcms_error(prefix + "mesh must be 2D or 3D");
  }
  if (!pcms::detail::SameCoordinateSystem(space.GetCoordinateSystem(),
                                          csys::CreateCartesian(mesh.dim()))) {
    throw pcms_error(prefix +
                     "space must use the Cartesian coordinate system of its "
                     "mesh dimension; map other systems explicitly");
  }
  if (mesh.family() != OMEGA_H_SIMPLEX) {
    throw pcms_error(prefix +
                     "mesh must be a simplex (triangle/tetrahedron) mesh");
  }
  const int order = layout->GetOrder();
  if (order < min_order || order > max_order) {
    throw pcms_error(
      prefix + "space must be order-" + std::to_string(min_order) +
      (max_order > min_order ? " or order-" + std::to_string(max_order)
                             : std::string()));
  }
  return {std::move(layout), std::move(geometry)};
}

// Measure of the element of each of `points`, which lie in the elements
// `owners` names, on an affine simplex geometry: |det J| / Dim! from the
// geometry's cell Jacobians.
inline Kokkos::View<Real*, DeviceMemorySpace> ElementMeasuresAtPoints(
  const Geometry& geometry, CoordinateView<DeviceMemorySpace> points,
  const OwnerHints& owners, const char* context)
{
  const auto* omega_h_geometry =
    dynamic_cast<const omega_h::OmegaHGeometry*>(&geometry);
  if (omega_h_geometry == nullptr) {
    throw pcms_error(std::string(context) +
                     ": element measures need an Omega_h geometry");
  }
  const LO n = static_cast<LO>(points.GetValues().extent(0));
  PCMS_ALWAYS_ASSERT(owners.entities.extent(0) == static_cast<size_t>(n));
  const auto element_measures =
    Omega_h::measure_elements_real(&omega_h_geometry->GetMesh());
  const auto element_ids = owners.entities;
  Kokkos::View<Real*, DeviceMemorySpace> measures(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "element_measures"), n);
  Kokkos::parallel_for(
    "ElementMeasuresAtPoints",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      measures(i) = Kokkos::fabs(element_measures[element_ids(i)]);
    });
  return measures;
}

// Map barycentric coordinates on a simplex (Dim+1 barycentric components) to
// the global point, given the simplex's Dim+1 vertex coordinates. The
// intersection uses it on the sub-simplices r3d produces.
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
// op(simplex, measure) for every piece with |measure| > eps.
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

// loop over the sub-simplices (triangles in 2D, tets in 3D) that partition each
// target element's intersection with the source mesh
template <int Dim, typename SimplexOp>
OMEGA_H_INLINE void ForEachIntersectionSubsimplex(
  const int elm, const omega_h::IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_elems2nodes, const Omega_h::LOs& src_elems2nodes,
  SimplexOp&& op)
{
  auto tgt_elm_vert_coords =
    omega_h::get_vert_coords_of_elem<Dim>(tgt_coords, tgt_elems2nodes, elm);
  const int start = intersection.tgt2src_offsets[elm];
  const int end = intersection.tgt2src_offsets[elm + 1];

  for (int i = start; i < end; ++i) {
    const int current_src_elm = intersection.tgt2src_indices[i];
    auto src_elm_vert_coords = omega_h::get_vert_coords_of_elem<Dim>(
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

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP
