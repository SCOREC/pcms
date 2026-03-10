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

  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();

  /*
     We make two 2D meshes on the unit square.
     Target mesh is coarse   : 4 × 4 = 32 triangles
     Source mesh is refined  : 8 × 8 = 128 triangles
     This ensures non-trivial intersection patterns.
  */

  int tgt_n = 4; // target elements per direction
  int src_n = 8; // source elements per direction (finer)

  Omega_h::Mesh tgt_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX,
                       /*x-length=*/1, /*y-length=*/1, /*z-length=*/0,
                       /*nx=*/tgt_n, /*ny=*/tgt_n, /*nz=*/0.0,
                       /*symmetric=*/false);

  Omega_h::Mesh src_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, src_n, src_n, 0, false);

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

    r3d::Polytope<2> P2;
    r3d::intersect_simplices(P2, B, A);
    double vol2 = r3d::measure(P2);
    REQUIRE(P2.nverts == 3);
    REQUIRE(Omega_h::are_close(vol2, 0.5 * 0.5 * 0.5));
  }
}
