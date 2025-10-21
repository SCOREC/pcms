#ifndef PCMS_INTERPOLATOR_LOAD_VECTOR_INTEGRATOR_HPP
#define PCMS_INTERPOLATOR_LOAD_VECTOR_INTEGRATOR_HPP

#include <Omega_h_shape.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <vector>
#include "mesh_intersection.hpp"
#include <Kokkos_MathematicalFunctions.hpp>

// computes the load vector for each element
// This routine is only applicable linear elements
// It calculates the rhs of the conservative field transfer formula from source
// to the target mesh

static const double EPS_AREA = 1e-16;
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

template <int order>
struct IntegrationData
{
  Kokkos::View<MeshField::Vector3*> bary_coords;
  Kokkos::View<Omega_h::Real*> weights;

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

  int size() const { return bary_coords.extent(0); }
};

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
      //	printf("number of intersections : %d\n", end-start);
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
        //			printf("No of verts : %d\n", nverts);
        // decompose if more than 3 vertices
        //			printf(" area from polytope : %f\n",
        //r3d::measure(poly));
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

          if (area <= EPS_AREA)
            continue; // drops duplicates and colinear/degenerates
                      //				printf(" area = %f\n", area);

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
        // printf(" sum_area = %.16e\n", sum_area);
        // printf(" poly_area = %.16e\n", poly_area);
        // OMEGA_H_CHECK(nearly_equal(sum_area, poly_area));
      }

      for (int j = 0; j < 3; ++j) {
        elmLoadVector(elm * 3 + j) = part_integration[j];
      }
    });

  return elmLoadVector;
}

struct Errors
{
  double proj_err;
  double cons_err;
};

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

          const double EPS_AREA = Kokkos::fmax(abs_tol, rel_tol * area);
          if (area <= EPS_AREA)
            continue; // drops duplicates and colinear/degenerates
                      //				printf(" area = %f\n", area);

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
