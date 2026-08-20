#ifndef PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP
#define PCMS_TRANSFER_OMEGA_H_FORM_INTEGRATOR_UTILS_HPP

#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"
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
  std::shared_ptr<const CoordinateSystem> coordinate_system,
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
  const Omega_h::Mesh& mesh = layout->GetMesh();
  if (mesh.dim() != 2 && mesh.dim() != 3) {
    throw pcms_error(std::string(context) + ": " + role +
                     " mesh must be 2D or 3D");
  }
  if (!SameCoordinateSystem(coordinate_system,
                            csys::Cartesian::Create(mesh.dim()))) {
    throw pcms_error(std::string(context) + ": " + role +
                     " space must use the Cartesian coordinate system of its "
                     "mesh dimension; map other systems explicitly");
  }
  if (mesh.family() != OMEGA_H_SIMPLEX) {
    throw pcms_error(std::string(context) + ": " + role +
                     " mesh must be a simplex (triangle/tetrahedron) mesh");
  }
}

// Strict order-1 check (used where only P1 is supported, e.g. Monte-Carlo RHS).
inline void CheckOmegaHScalarP1Layout(
  std::shared_ptr<const CoordinateSystem> coordinate_system,
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
  std::shared_ptr<const CoordinateSystem> coordinate_system,
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

[[nodiscard]] OMEGA_H_INLINE int RemoveDuplicateVerticesAndFixLinks(
  r3d::Polytope<2>& poly, const double tol = 1e-12)
{
  const int old_n = poly.nverts;
  int new_n = 0;

  for (int i = 0; i < old_n; ++i) {
    const auto& pi = poly.verts[i].pos;
    bool dup = false;
    for (int j = 0; j < new_n; ++j) {
      const auto& pj = poly.verts[j].pos;
      if (Kokkos::fabs(pi[0] - pj[0]) < tol &&
          Kokkos::fabs(pi[1] - pj[1]) < tol) {
        dup = true;
        break;
      }
    }
    if (!dup) {
      poly.verts[new_n] = poly.verts[i];
      ++new_n;
    }
  }
  poly.nverts = new_n;

  if (new_n >= 3) {
    double cx = 0.0;
    double cy = 0.0;
    for (int i = 0; i < new_n; ++i) {
      cx += poly.verts[i].pos[0];
      cy += poly.verts[i].pos[1];
    }
    cx /= new_n;
    cy /= new_n;

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

    for (int i = 0; i < new_n; ++i) {
      const int prev = (i - 1 + new_n) % new_n;
      const int next = (i + 1) % new_n;
      poly.verts[i].pnbrs[0] = prev;
      poly.verts[i].pnbrs[1] = next;
    }

    double area = r3d::measure(poly);
    if (area < 0.0) {
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
    for (int i = 0; i < new_n; ++i) {
      poly.verts[i].pnbrs[0] = -1;
      poly.verts[i].pnbrs[1] = -1;
    }
  }

  return new_n;
}

// 2D: fan the clipped intersection polygon into triangles anchored at vertex 0,
// invoking op(sub_triangle, src_elm, area) for each non-degenerate piece.
template <typename SimplexOp>
OMEGA_H_INLINE void ForEachIntersectionSubtriangleImpl(
  const int elm, const IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_elems2nodes, const Omega_h::LOs& src_elems2nodes,
  SimplexOp&& op)
{
  auto tgt_elm_vert_coords =
    get_vert_coords_of_elem<2>(tgt_coords, tgt_elems2nodes, elm);
  const int start = intersection.tgt2src_offsets[elm];
  const int end = intersection.tgt2src_offsets[elm + 1];

  for (int i = start; i < end; ++i) {
    const int current_src_elm = intersection.tgt2src_indices[i];
    auto src_elm_vert_coords =
      get_vert_coords_of_elem<2>(src_coords, src_elems2nodes, current_src_elm);
    r3d::Polytope<2> poly;
    r3d::intersect_simplices(poly, tgt_elm_vert_coords, src_elm_vert_coords);
    auto nverts = RemoveDuplicateVerticesAndFixLinks(poly, 1e-12);
    auto poly_area = r3d::measure(poly);

    for (int j = 1; j < nverts - 1; ++j) {
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

      const double eps_area =
        PCMS_INTERSECTION_ABS_TOL + PCMS_INTERSECTION_REL_TOL * poly_area;
      if (area <= eps_area) {
        continue;
      }

      op(tri_coords, current_src_elm, area);
    }
  }
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

// 3D: star-decompose the clipped intersection polyhedron into tetrahedra from
// its centroid. The centroid lies strictly inside the convex intersection, so
// lifting each boundary-face triangle (enumerated by
// ForEachPolytopeFaceTriangle) to the centroid tiles the polyhedron without
// overlap; op(sub_tet, src_elm, volume) fires for each non-degenerate piece.
template <typename SimplexOp>
OMEGA_H_INLINE void ForEachIntersectionSubtetImpl(
  const int elm, const IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_elems2nodes, const Omega_h::LOs& src_elems2nodes,
  SimplexOp&& op)
{
  auto tgt_elm_vert_coords =
    get_vert_coords_of_elem<3>(tgt_coords, tgt_elems2nodes, elm);
  const int start = intersection.tgt2src_offsets[elm];
  const int end = intersection.tgt2src_offsets[elm + 1];

  for (int i = start; i < end; ++i) {
    const int current_src_elm = intersection.tgt2src_indices[i];
    auto src_elm_vert_coords =
      get_vert_coords_of_elem<3>(src_coords, src_elems2nodes, current_src_elm);
    r3d::Polytope<3> poly;
    r3d::intersect_simplices(poly, tgt_elm_vert_coords, src_elm_vert_coords);
    if (poly.nverts < 4) {
      continue;
    }
    const double poly_vol = Kokkos::fabs(r3d::measure(poly));
    const double eps_vol =
      PCMS_INTERSECTION_ABS_TOL + PCMS_INTERSECTION_REL_TOL * poly_vol;

    // Centroid of the (convex) intersection polyhedron: interior apex.
    Omega_h::Vector<3> apex = {0.0, 0.0, 0.0};
    for (int v = 0; v < poly.nverts; ++v) {
      apex[0] += poly.verts[v].pos[0];
      apex[1] += poly.verts[v].pos[1];
      apex[2] += poly.verts[v].pos[2];
    }
    apex[0] /= poly.nverts;
    apex[1] /= poly.nverts;
    apex[2] /= poly.nverts;

    // Lift each boundary-face triangle to the interior centroid to form a tet.
    ForEachPolytopeFaceTriangle(poly, [&](const r3d::Vector<3>& a,
                                          const r3d::Vector<3>& b,
                                          const r3d::Vector<3>& c) {
      Omega_h::Few<Omega_h::Vector<3>, 4> tet_coords;
      tet_coords[0] = apex;
      tet_coords[1] = {a[0], a[1], a[2]};
      tet_coords[2] = {b[0], b[1], b[2]};
      tet_coords[3] = {c[0], c[1], c[2]};

      Omega_h::Few<Omega_h::Vector<3>, 3> basis;
      basis[0] = tet_coords[1] - tet_coords[0];
      basis[1] = tet_coords[2] - tet_coords[0];
      basis[2] = tet_coords[3] - tet_coords[0];

      const Omega_h::Real vol =
        Kokkos::fabs(Omega_h::tet_volume_from_basis(basis));
      if (vol > eps_vol) {
        op(tet_coords, current_src_elm, vol);
      }
    });
  }
}

// Dimension-generic driver over the sub-simplices (triangles in 2D, tets in 3D)
// that tile each target element's intersection with the source mesh. Invokes
// op(sub_simplex_coords, src_elm, measure) for every non-degenerate piece.
template <int Dim, typename SimplexOp>
OMEGA_H_INLINE void ForEachIntersectionSubsimplex(
  const int elm, const IntersectionResults& intersection,
  const Omega_h::Reals& tgt_coords, const Omega_h::Reals& src_coords,
  const Omega_h::LOs& tgt_elems2nodes, const Omega_h::LOs& src_elems2nodes,
  SimplexOp&& op)
{
  if constexpr (Dim == 3) {
    ForEachIntersectionSubtetImpl(elm, intersection, tgt_coords, src_coords,
                                  tgt_elems2nodes, src_elems2nodes, op);
  } else {
    ForEachIntersectionSubtriangleImpl(elm, intersection, tgt_coords,
                                       src_coords, tgt_elems2nodes,
                                       src_elems2nodes, op);
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
