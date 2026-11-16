#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/utility/mesh_geometry.h>

TEST_CASE("pcms::get_entity_centroids supports simplex mesh entities")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  SECTION("2D edge centroids")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 4, 4, 0, false);

    auto centroids = pcms::get_entity_centroids(mesh, Omega_h::EDGE);
    REQUIRE(centroids.size() == mesh.nedges() * mesh.dim());
  }

  SECTION("2D face centroids")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 4, 4, 0, false);

    auto centroids = pcms::get_entity_centroids(mesh, Omega_h::FACE);
    REQUIRE(centroids.size() == mesh.nfaces() * mesh.dim());
  }

  SECTION("3D edge centroids")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 3, 3, 3, false);

    auto centroids = pcms::get_entity_centroids(mesh, Omega_h::EDGE);
    REQUIRE(centroids.size() == mesh.nedges() * mesh.dim());
  }

  SECTION("3D element centroids")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 3, 3, 3, false);

    auto centroids = pcms::get_entity_centroids(mesh, Omega_h::REGION);
    REQUIRE(centroids.size() == mesh.nelems() * mesh.dim());
  }

  SECTION("Vertex coordinates passthrough")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 2, 2, 0, false);

    auto centroids = pcms::get_entity_centroids(mesh, Omega_h::VERT);
    auto coords = mesh.coords();
    REQUIRE(centroids.size() == coords.size());
  }
}

TEST_CASE("pcms::distance_squared handles 2D and 3D")
{
  SECTION("2D distance")
  {
    const Omega_h::Real p1[3] = {1.0, 2.0, 99.0};
    const Omega_h::Real p2[3] = {4.0, 6.0, -5.0};
    REQUIRE(pcms::distance_squared(p1, p2, 2) == Catch::Approx(25.0));
  }

  SECTION("3D distance")
  {
    const Omega_h::Real p1[3] = {1.0, 2.0, 3.0};
    const Omega_h::Real p2[3] = {4.0, 6.0, 8.0};
    REQUIRE(pcms::distance_squared(p1, p2, 3) == Catch::Approx(50.0));
  }
}
