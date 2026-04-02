#include <pcms/field_transfer/mesh_intersection.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_vtk.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Mesh intersection test with source and target", "[intersection]")
{

  Omega_h::Library lib;

  Omega_h::Reals coords_target({
    0.0, 0.0, // v0
    1.0, 0.0, // v1
    1.0, 1.0, // v2
    0.0, 1.0  // v3
  });

  // Target Mesh with two triangles
  // Two triangles, CCW
  // T0: (v0,v1,v3) = (0,1,3)
  // T1: (v1,v2,v3) = (1,2,3)
  Omega_h::LOs ev2v_target({0, 1, 3, 1, 2, 3});

  Omega_h::Mesh tgt_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&tgt_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_target, coords_target);

  // 2D coordinates (x,y) : 4 vertices
  Omega_h::Reals coords_source({
    0.0, 0.0, // v0
    0.5, 0.0, // v1
    1.0, 0.0, // v2
    0.0, 0.5, // v3
    0.5, 0.5, // v4
    1.0, 0.5, // v5
    0.0, 1.0, // v6
    0.5, 1.0, // v7
    1.0, 1.0  // v8
  });

  // TARGET triangulations
  Omega_h::LOs ev2v_source(
    {0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 3, 4, 7, 3, 7, 6, 4, 5, 8, 4, 8, 7});
  Omega_h::Mesh src_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&src_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_source, coords_source);

  REQUIRE(src_mesh.dim() == 2);
  REQUIRE(tgt_mesh.dim() == 2);

  int num_tgt_elems = tgt_mesh.nelems();

  REQUIRE(src_mesh.dim() == 2);
  REQUIRE(tgt_mesh.dim() == 2);

  const int nsrc = src_mesh.nelems();
  const int ntgt = tgt_mesh.nelems();

  REQUIRE(nsrc > 0);
  REQUIRE(ntgt > 0);

  const auto src_coords = src_mesh.coords();
  const auto tgt_coords = tgt_mesh.coords();

  const auto src_faces2verts =
    src_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto tgt_faces2verts =
    tgt_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  auto tgt_elm_area = Omega_h::measure_elements_real(&tgt_mesh);

  // ============================
  // Compute intersection mapping
  // ============================
  auto intersection = intersectTargets(src_mesh, tgt_mesh);

  auto offsets = Omega_h::HostRead(intersection.tgt2src_offsets);
  auto indices = Omega_h::HostRead(intersection.tgt2src_indices);

  SECTION("Offsets and indices sizes")
  {
    REQUIRE(offsets.size() == static_cast<std::size_t>(ntgt + 1));
    REQUIRE(indices.size() == static_cast<std::size_t>(offsets[ntgt]));
  }

  SECTION("All intersection elements indices valid")
  {
    for (int i = 0; i < indices.size(); ++i) {
      REQUIRE(indices[i] >= 0);
      REQUIRE(indices[i] < nsrc);
    }
  }

  SECTION("Each target element should intersect at least one source element")
  {
    for (int t = 0; t < ntgt; ++t) {
      REQUIRE(offsets[t + 1] - offsets[t] > 0);
    }
  }

  SECTION("Total intersection area matches target mesh area")
  {
    Omega_h::Write<Omega_h::Real> total_intersected_area(ntgt, 0.0);

    Omega_h::parallel_for(
      ntgt, OMEGA_H_LAMBDA(int t) {
        auto tgt_vert_coords =
          get_vert_coords_of_elem(tgt_coords, tgt_faces2verts, t);
        int start = intersection.tgt2src_offsets[t];
        int end = intersection.tgt2src_offsets[t + 1];

        Omega_h::Real sum_area = 0.0;

        for (int i = start; i < end; ++i) {
          int sid = intersection.tgt2src_indices[i];
          auto src_vert_coords =
            get_vert_coords_of_elem(src_coords, src_faces2verts, sid);

          r3d::Polytope<2> poly;
          r3d::intersect_simplices(poly, tgt_vert_coords, src_vert_coords);

          sum_area += r3d::measure(poly);
        }
        total_intersected_area[t] = sum_area;
        printf("element = %d, num source intersections = %d, sum area = %f\n",
               t, end - start, sum_area);
      });

    auto expected = Omega_h::HostRead(tgt_elm_area);
    auto computed = Omega_h::HostRead(Omega_h::read(total_intersected_area));

    double tol = 1e-6;

    for (int t = 0; t < ntgt; ++t) {
      CAPTURE(t, expected[t], computed[t]);
      CHECK_THAT(expected[t], Catch::Matchers::WithinAbs(computed[t], tol));
    }
  }

  SECTION("Simple r3d intersection sanity test")
  {

    r3d::Few<r3d::Vector<2>, 3> A = {r3d::Vector<2>{{0.2, 0.2}},
                                     r3d::Vector<2>{{0.7, 0.2}},
                                     r3d::Vector<2>{{0.2, 0.7}}};

    r3d::Few<r3d::Vector<2>, 3> B = {r3d::Vector<2>{{0.0, 0.0}},
                                     r3d::Vector<2>{{1.0, 0.0}},
                                     r3d::Vector<2>{{0.0, 1.0}}};

    r3d::Polytope<2> P1;
    r3d::intersect_simplices(P1, A, B);
    double vol1 = r3d::measure(P1);
    REQUIRE(P1.nverts == 3);
    REQUIRE(Omega_h::are_close(vol1, 0.5 * 0.5 * 0.5));
  }
}
