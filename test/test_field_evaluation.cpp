#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include "field_test_utils.h"
#include <array>
#include <cmath>
#include <vector>

using pcms::Real;

// Non-linear test function used to exercise quadratic interpolation.
KOKKOS_INLINE_FUNCTION static Real sin_f(Real x, Real y)
{
  return std::sin(20 * x * y) / 2 + 0.5;
}

// Standard interior test points shared across all evaluation tests.
static const std::vector<Real> kEvalCoords = {
  0.7681, 0.886,  0.5337, 0.5205,   0.8088, 0.1513, 0.13,
  0.43,   0.5484, 0.8263, 0.006119, 0.8642, 0.5889, 0.5622,
  0.9268, 0.1749, 0.2615, 0.1468,   0.9793, 0.9612,
};

TEST_CASE("evaluate linear 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory->CreateFunction<Real>();

  pcms::test::SetField(
    field.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  pcms::test::CheckEvaluation(
    factory, field, pcms::test::StandardEvalCoords2D(),
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
}

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("evaluate quadratic 2d meshfields_field")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 2, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  // Quadratic DOF holders span vertices and edge midpoints; the layout's DOF
  // coordinates cover both, so SetField samples sin_f at every holder.
  auto field = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field, OMEGA_H_LAMBDA(Real x, Real y) { return sin_f(x, y); });

  pcms::test::CheckEvaluation(
    factory, field, kEvalCoords,
    OMEGA_H_LAMBDA(Real x, Real y) { return std::sin(20 * x * y) / 2 + 0.5; },
    1.0e-2);
}
#endif

#ifdef PCMS_ENABLE_MESHFIELDS
// Recomputes the expected DOF-holder coordinates from mesh geometry exactly as
// the layout's isoparametric map does for affine simplex elements, accumulated
// per entity dimension in the same block order the layout stores:
//   dim 0 (vertex) holders -> the vertex's own coordinate,
//   dim 1 (edge)   holders -> the midpoint of the edge's two endpoint vertices.
// Returns a flat node-major vector matching GetDOFHolderCoordinates() layout.
static std::vector<Real> ComputeManualHolderCoords(
  Omega_h::Mesh& mesh, const std::array<int, 4>& nodes_per_dim)
{
  const int dim = mesh.dim();
  auto mesh_coords_h = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  auto edge2vtx_h = Omega_h::HostRead<Omega_h::LO>(mesh.ask_verts_of(1));

  std::vector<Real> expected;
  for (int e_dim = 0; e_dim <= dim; ++e_dim) {
    if (nodes_per_dim[static_cast<size_t>(e_dim)] != 1)
      continue;
    const int num_ents = mesh.nents(e_dim);
    for (int ent = 0; ent < num_ents; ++ent) {
      if (e_dim == 0) {
        for (int d = 0; d < dim; ++d)
          expected.push_back(mesh_coords_h[ent * dim + d]);
      } else { // e_dim == 1: edge midpoint
        const int a = edge2vtx_h[2 * ent + 0];
        const int b = edge2vtx_h[2 * ent + 1];
        for (int d = 0; d < dim; ++d)
          expected.push_back(
            0.5 * (mesh_coords_h[a * dim + d] + mesh_coords_h[b * dim + d]));
      }
    }
  }
  return expected;
}

// Builds a MeshFields layout of the given order on a small known box and checks
// that GetDOFHolderCoordinates() matches the manual reference computation.
static void CheckMeshFieldsDofHolderCoordsAgainstManual(Omega_h::Library& lib,
                                                        int order)
{
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 2, 3, 0, false);

  std::array<int, 4> nodes_per_dim{};
  if (order == 1)
    nodes_per_dim = {1, 0, 0, 0};
  else if (order == 2)
    nodes_per_dim = {1, 1, 0, 0};

  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);
  auto layout = factory->GetLayout();

  // Sanity-check the block layout: vertices then edges.
  if (order == 2) {
    const auto offsets = layout->GetEntOffsets();
    const int nverts = mesh.nents(0);
    REQUIRE(offsets[0] == 0);
    REQUIRE(offsets[1] == static_cast<size_t>(nverts));
  }

  const std::vector<Real> expected =
    ComputeManualHolderCoords(mesh, nodes_per_dim);

  auto coords = layout->GetDOFHolderCoordinates().GetValues();
  const int dim = mesh.dim();
  REQUIRE(coords.extent(1) == static_cast<size_t>(dim));
  REQUIRE(coords.extent(0) == expected.size() / static_cast<size_t>(dim));

  auto coords_h = pcms::test::CopyCoordinatesToHost(
    coords, static_cast<int>(coords.extent(0)), dim);

  size_t idx = 0;
  for (size_t r = 0; r < coords_h.extent(0); ++r) {
    for (int d = 0; d < dim; ++d) {
      INFO("row " << r << " dim " << d);
      REQUIRE(coords_h(r, d) == Catch::Approx(expected[idx]));
      ++idx;
    }
  }
}

TEST_CASE("MeshFields: linear dof-holder coordinates match manual reference")
{
  auto lib = Omega_h::Library{};
  CheckMeshFieldsDofHolderCoordsAgainstManual(lib, 1);
}

TEST_CASE("MeshFields: quadratic dof-holder coordinates match manual reference")
{
  auto lib = Omega_h::Library{};
  CheckMeshFieldsDofHolderCoordsAgainstManual(lib, 2);
}
#endif

TEST_CASE("evaluate quadratic 2d omega_h_field throws")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);

  REQUIRE_THROWS_AS(pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 2, 1, pcms::CoordinateSystem::Cartesian, "global",
                      pcms::LagrangeFunctionSpace::Backend::OmegaH),
                    pcms::pcms_error);
}
